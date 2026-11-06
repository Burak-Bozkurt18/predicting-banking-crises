# =============================================================================
# Early-warning models: estimation, forecasting, evaluation
#
# Design
#   * ROC-AUC is the main metric (hyperparameter selection, DeLong comparison).
#   * PR-AUC is reported as a supplementary metric with country-cluster
#     bootstrap confidence intervals.
#   * Hyperparameters are re-tuned every `retune_every` years using only
#     information available at the retuning date.
#   * Every analysis (main + robustness checks) is a "scenario" that runs
#     through the SAME functions. Expensive steps (tuning, forecasting,
#     bootstrap) are cached on disk. The cache key contains a hash of the
#     data, the configuration AND the code of the model definition, so:
#       - changing one model's workflow/grid recomputes only that model,
#       - adding a robustness check does not recompute anything else,
#       - rerunning an unchanged script is (almost) instant.
# =============================================================================

# 1 Load Packages =============================================================

library(tidyverse)
library(tidymodels)
library(modelsummary)
library(kableExtra)
library(pROC)
library(foreach)
library(doParallel)
library(themis)
library(patchwork)

# 2 Configuration =============================================================

last_known_crisis_year <- 2025
first_year             <- 1970
first_test_year        <- 2005
retune_every           <- 5     # re-tune hyperparameters every 5 years
n_assess               <- 10    # number of assessment years per tuning round
n_boot                 <- 1000  # bootstrap replicates for PR-AUC intervals

outcomes   <- paste0("precrisis", 1:4)
models     <- c("Logit", "RF", "SVM", "MLP")
ml_models  <- setdiff(models, "Logit")
grid_size  <- c(RF = 30, SVM = 30, MLP = 15)

horizon_of    <- function(outcome) as.integer(sub("precrisis", "", outcome))
max_test_year <- function(outcome) last_known_crisis_year - horizon_of(outcome)
k_label       <- function(outcome) paste0("K = ", horizon_of(outcome))

# Retuning dates: 2005, 2010, 2015, 2020
retune_years   <- seq(first_test_year, last_known_crisis_year - 1, by = retune_every)
retune_year_of <- function(y) max(retune_years[retune_years <= y])

cache_dir <- "cache"
table_dir <- "output/tables"
fig_dir   <- "output/figures"
walk(c(cache_dir, table_dir, fig_dir), dir.create, showWarnings = FALSE)

# 3 Helpers ===================================================================

# Evaluate `expr` only if `file` does not exist; otherwise read it from disk
cached <- function(file, expr) {
  if (file.exists(file)) {
    message("cache hit: ", basename(file))
    return(read_rds(file))
  }
  res <- expr
  write_rds(res, file, compress = "gz")
  res
}

cache_file <- function(type, model, outcome, spec, key) {
  file.path(cache_dir,
            sprintf("%s_%s_%s_%s_%s.rds", type, model, outcome, spec,
                    substr(key, 1, 10)))
}

write_table <- function(x, name) {
  writeLines(as.character(x), file.path(table_dir, paste0(name, ".tex")))
  invisible(x)
}

stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.1 ~ "*",
            TRUE ~ "")
}

# 4 Data and scenarios ========================================================

panel <- read_rds("data/final/panel.rds")

# One cleaned panel for everything (logit and ML see identical data)
panel_base <- panel |>
  filter(crisis != 1) |>
  mutate(
    advanced = factor(advanced),
    region   = factor(region),
    across(all_of(outcomes), factor)
  )

stopifnot(identical(levels(panel_base$precrisis1), c("0", "1")))

panel_logit <- panel_base   # used for the in-sample logits below

specs <- list(
  baseline = c("advanced", "cgdppriv_growth"),
  broad    = c("advanced", "cgdppriv_growth", "govcgdp_growth",
               "rgdpgrowth", "inflation"),
  full     = c("advanced", "cgdppriv_growth", "govcgdp_growth",
               "rgdpgrowth", "inflation", "ltd", "bcagdp",
               "bmgdpgrowth", "bmtr", "sprr")
)

full_vars_rc2 <- c("advanced", "tlpriv_rgrowth", "govcgdp_growth",
                   "rgdpgrowth", "inflation", "ltd", "bcagdp",
                   "bm_rgrowth", "bmtr", "sprr")

full_vars_rc3 <- c("region", "cgdppriv_growth", "govcgdp_growth",
                   "rgdpgrowth", "inflation", "ltd", "bcagdp",
                   "bmgdpgrowth", "bmtr", "sprr")

# A scenario = a set of specifications (+ optional sample restriction).
# Adding a robustness check means adding one entry here.
scenarios <- list(
  main = list(specs = specs, complete_on = NULL),
  # RC1: baseline and broad on the same sample as the full specification
  rc1  = list(specs = specs[c("baseline", "broad")], complete_on = specs$full),
  # RC2: real growth variables
  rc2  = list(specs = list(full = full_vars_rc2), complete_on = NULL),
  # RC3: region instead of advanced
  rc3  = list(specs = list(full = full_vars_rc3), complete_on = NULL)
)
# RC4 (excluding COVID years) needs no new estimation: it re-evaluates the
# main predictions up to 2019 (see section 11).

# Estimation data for one scenario x outcome x specification.
# Output columns: iso3c, year, truth (factor "0"/"1"), predictors
get_panel <- function(scn, outcome, spec) {
  vars <- scn$specs[[spec]]
  panel_base |>
    select(iso3c, year, all_of(outcome),
           all_of(unique(c(vars, scn$complete_on)))) |>
    drop_na() |>
    select(iso3c, year, all_of(outcome), all_of(vars)) |>
    rename(truth = all_of(outcome))
}

# 5 Theme =====================================================================

theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "serif", size = 12),
      axis.text = element_text(size = 9, colour = "black"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      panel.grid = element_line(linetype = "dashed")
    )
)

# 6 Sample composition table ==================================================

fmt_range <- function(x) {
  if (length(x) == 0) return("--")
  r <- range(x)
  if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], "--", r[2])
}

create_sample_summary <- function(outcome, spec) {
  h         <- horizon_of(outcome)
  train_end <- first_test_year - 1 - h
  test_end  <- max_test_year(outcome)
  
  s <- get_panel(scenarios$main, outcome, spec) |>
    mutate(y = truth == "1")
  
  train   <- filter(s, year >= first_year, year <= train_end)
  gap <- filter(s, year > train_end, year < first_test_year)
  test    <- filter(s, year >= first_test_year, year <= test_end)
  
  if (nrow(train) + nrow(gap) + nrow(test) != nrow(s)) {
    warning("Train + gap + test != full sample for ", outcome)
  }
  
  tibble(
    observations   = nrow(s),
    precrisisK     = sum(s$y),
    countries      = n_distinct(s$iso3c),
    advanced_share = (mean(as.integer(s$advanced)) - 1) * 100,
    train_years    = fmt_range(train$year),
    training_obs   = nrow(train),
    training_share = mean(train$y) * 100,
    gap_obs    = nrow(gap),
    test_years     = fmt_range(test$year),
    testing_obs    = nrow(test),
    testing_share  = mean(test$y) * 100
  )
}

tab_df <- expand_grid(specification = names(specs), outcome = outcomes) |>
  mutate(summary = map2(outcome, specification, create_sample_summary),
         specification = factor(specification, levels = names(specs)),
         Outcome = k_label(outcome)) |>
  unnest(summary) |>
  arrange(specification, outcome)

panel_index <- tab_df |>
  count(specification) |>
  mutate(label = str_to_title(as.character(specification))) %>%
  { setNames(.$n, .$label) }

tab_sample <- tab_df |>
  transmute(
    Outcome,
    obs       = prettyNum(observations, big.mark = ","),
    pre_obs   = prettyNum(precrisisK, big.mark = ","),
    countries = as.character(countries),
    advanced  = formatC(advanced_share, format = "f", digits = 1),
    train_yrs = train_years,
    train_obs = prettyNum(training_obs, big.mark = ","),
    train_pre = formatC(training_share, format = "f", digits = 1),
    gap_obs   = prettyNum(gap_obs, big.mark = ","),
    test_yrs  = test_years,
    test_obs  = prettyNum(testing_obs, big.mark = ","),
    test_pre  = formatC(testing_share, format = "f", digits = 1)
  ) |>
  kbl(
    format = "latex", booktabs = TRUE, align = paste0("l", strrep("r", 11)),
    escape = FALSE, linesep = "",
    col.names = c("Outcome", "Obs.", "PreCrisis obs.", "Countries",
                  "Advanced (\\%)", "Years", "Obs.", "PreCrisis (\\%)",
                  "Obs.", "Years", "Obs.", "PreCrisis (\\%)"),
    caption = "Sample composition by model specification and crisis horizon",
    label = "sample_summary"
  ) |>
  pack_rows(index = panel_index, bold = TRUE, italic = FALSE,
            latex_gap_space = "0.4em") |>
  add_header_above(c(" " = 1, "Full sample" = 4, "Training" = 3,
                     "Gap" = 1, "Testing" = 3)) |>
  kable_styling(latex_options = c("hold_position", "scale_down"),
                position = "center")

tab_sample

# 7 In-sample logits (full sample, descriptive) ===============================

gof_custom <- function(model) {
  mf        <- model.frame(model)
  rows      <- as.integer(rownames(mf))
  n_countries <- length(unique(panel_logit$iso3c[rows]))
  y         <- model.response(mf)
  ll_model  <- as.numeric(logLik(model))
  ll_null   <- as.numeric(logLik(glm(y ~ 1, family = binomial())))
  data.frame(`Num. Countries` = n_countries,
             `McFadden R²` = round(1 - ll_model / ll_null, 3))
}

models_insample <- map(
  outcomes,
  \(o) map(specs, \(rhs) glm(reformulate(rhs, response = o),
                             data = panel_logit, family = binomial()))
) |> setNames(outcomes)

walk2(
  list(c("precrisis1", "precrisis2"), c("precrisis3", "precrisis4")),
  c("Early-warning logits, 1-2 year horizons",
    "Early-warning logits, 3-4 year horizons"),
  \(os, ttl) {
    modelsummary(
      models_insample[os], title = ttl, gof_function = gof_custom,
      gof_omit = "RMSE|Log.Lik.|Std.Errors|AIC", coef_omit = "(Intercept)",
      stars = TRUE, vcov = ~iso3c, shape = "cbind", output = "latex"
    ) |> print()
  }
)

# 8 Model definitions =========================================================
# Everything that defines a model's architecture lives here. Editing one of
# these functions invalidates the cache of that model only.

make_rf_workflow <- function(data) {
  recipe(truth ~ ., data = data) |>
    step_rm(iso3c, year) |>
    workflow() |>
    add_model(
      rand_forest(trees = 500, mtry = tune(), min_n = tune()) |>
        set_engine("ranger", importance = "permutation") |>
        set_mode("classification")
    )
}

make_rf_grid <- function(data, size) {
  n_predictors <- ncol(data) - 3   # minus truth, iso3c, year
  grid_random(
    mtry(range = c(2L, n_predictors)),
    min_n(range = c(2L, 30L)),
    size = size
  )
}

make_svm_workflow <- function(data) {
  recipe(truth ~ ., data = data) |>
    step_rm(iso3c, year) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_predictors()) |>
    step_upsample(truth, over_ratio = 0.5) |>
    workflow() |>
    add_model(
      svm_rbf(cost = tune(), rbf_sigma = tune()) |>
        set_engine("kernlab") |>
        set_mode("classification")
    )
}

make_svm_grid <- function(data, size) {
  grid_space_filling(cost(), rbf_sigma(), size = size)
}

make_mlp_workflow <- function(data) {
  recipe(truth ~ ., data = data) |>
    step_rm(iso3c, year) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_predictors()) |>
    step_upsample(truth, over_ratio = 0.5) |>
    workflow() |>
    add_model(
      mlp(hidden_units = tune(), penalty = tune(),
          epochs = 50, activation = "relu") |>
        set_engine("brulee", stop_iter = 5) |>
        set_mode("classification")
    )
}

make_mlp_grid <- function(data, size) {
  grid_space_filling(hidden_units(), penalty(), size = size)
}

workflow_fns <- list(RF = make_rf_workflow,  SVM = make_svm_workflow,
                     MLP = make_mlp_workflow)
grid_fns     <- list(RF = make_rf_grid,      SVM = make_svm_grid,
                     MLP = make_mlp_grid)

# 9 Pipeline functions ========================================================

## 9.1 Tuning (ROC-AUC, rolling re-tuning) ====================================

# Folds for a retuning at the beginning of `retune_year`.
# Labels of year t become known in t + horizon, so assessment years end at
# retune_year - 1 - horizon, and the analysis window of each fold ends
# `horizon` + 1 years before its assessment year (gap).
make_cv_periods <- function(horizon, retune_year) {
  assess_end <- retune_year - 1 - horizon
  tibble(
    id              = paste0("Fold", seq_len(n_assess)),
    analysis_start  = first_year,
    assessment_year = (assess_end - n_assess + 1):assess_end
  ) |>
    mutate(analysis_end = assessment_year - 1 - horizon)
}

make_folds <- function(data, cv_periods) {
  splits <- pmap(
    cv_periods,
    \(id, analysis_start, assessment_year, analysis_end) {
      make_splits(
        x = list(
          analysis   = which(data$year >= analysis_start &
                               data$year <= analysis_end),
          assessment = which(data$year == assessment_year)
        ),
        data = data
      )
    }
  )
  manual_rset(splits = splits, ids = cv_periods$id)
}

tune_once <- function(data, model, horizon, retune_year) {
  message("  tuning ", model, ", retune year ", retune_year)
  set.seed(222)
  folds <- make_folds(data, make_cv_periods(horizon, retune_year))
  wf    <- workflow_fns[[model]](data)
  grid  <- grid_fns[[model]](data, grid_size[[model]])
  
  res <- tune_grid(
    wf, resamples = folds, grid = grid,
    metrics = metric_set(roc_auc),
    control = control_grid(save_pred = FALSE, verbose = FALSE)
  )
  best <- show_best(res, metric = "roc_auc", n = 1)
  
  tibble(
    retune_year = retune_year,
    params      = list(select_best(res, metric = "roc_auc")),
    cv_roc_auc  = best$mean,
    cv_n_folds  = best$n
  )
}

tune_model <- function(data, model, horizon, tune_years) {
  map(tune_years, \(R) tune_once(data, model, horizon, R)) |> list_rbind()
}

## 9.2 Forecasting ============================================================

# Expanding window, gap of `horizon` years, parameters of the most recent
# retuning date.
forecast_ml <- function(data, tuned, model, horizon, forecast_years) {
  set.seed(222)
  wf_template <- workflow_fns[[model]](data)
  
  map(forecast_years, \(y) {
    message("Forecasting ", model, ", year ", y)
    params <- tuned$params[[which(tuned$retune_year == retune_year_of(y))]]
    train  <- filter(data, year >= first_year, year <= y - 1 - horizon)
    test   <- filter(data, year == y)
    if (nrow(test) == 0) return(NULL)
    
    fit <- finalize_workflow(wf_template, params) |> fit(data = train)
    test |>
      transmute(iso3c, year, truth,
                .pred_1 = predict(fit, new_data = test, type = "prob")$.pred_1)
  }) |> list_rbind()
}

forecast_logit <- function(data, spec_vars, horizon, forecast_years) {
  formula <- reformulate(spec_vars, response = "truth")
  
  map(forecast_years, \(y) {
    train <- filter(data, year >= first_year, year <= y - 1 - horizon)
    test  <- filter(data, year == y)
    if (nrow(test) == 0) return(NULL)
    
    fit <- glm(formula, data = train, family = binomial())
    test |>
      transmute(iso3c, year, truth,
                .pred_1 = predict(fit, newdata = test, type = "response"))
  }) |> list_rbind()
}

## 9.3 Cache keys =============================================================
# The key changes when the data, configuration or model code changes.
# Helper functions not listed in `code` (e.g. `get_panel`) are not tracked;
# delete the cache folder if you change those.

tune_key <- function(data, model, horizon, max_year) {
  rlang::hash(list(
    data  = data, model = model, horizon = horizon, max_year = max_year,
    cfg   = list(retune_years, n_assess, first_year, grid_size[[model]]),
    code  = map(list(workflow_fns[[model]], grid_fns[[model]], make_cv_periods,
                     make_folds, tune_once, tune_model), deparse)
  ))
}

## 9.4 Predictions for one scenario cell ======================================

get_preds <- function(scn_name, model, outcome, spec) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", scn_name, " | ", model,
          " | ", outcome, " | ", spec)
  scn   <- scenarios[[scn_name]]
  data  <- get_panel(scn, outcome, spec)
  h     <- horizon_of(outcome)
  years <- first_test_year:max_test_year(outcome)
  
  if (model == "Logit") {
    return(forecast_logit(data, scn$specs[[spec]], h, years))
  }
  
  tk <- tune_key(data, model, h, max(years))
  tuned <- cached(
    cache_file("tune", model, outcome, spec, tk),
    tune_model(data, model, h, retune_years[retune_years <= max(years)])
  )
  
  fk <- rlang::hash(list(tk, years, deparse(forecast_ml)))
  cached(
    cache_file("oos", model, outcome, spec, fk),
    forecast_ml(data, tuned, model, h, years)
  )
}

# All models x outcomes x specifications of one scenario.
# Long tibble: scenario, model, outcome, specification, preds (list-column)
run_scenario <- function(scn_name, which_models = models) {
  scn <- scenarios[[scn_name]]
  expand_grid(model = which_models, outcome = outcomes,
              specification = names(scn$specs)) |>
    mutate(
      scenario = scn_name,
      preds = pmap(list(model, outcome, specification),
                   \(m, o, s) get_preds(scn_name, m, o, s)),
      .before = 1
    )
}

# 10 Run scenarios ============================================================
# Tip: for quick iterations run only run_scenario("main") or a single model,
# e.g. run_scenario("main", which_models = c("Logit", "RF")).

cl <- makePSOCKcluster(parallel::detectCores() - 1)
registerDoParallel(cl)

oos_all <- map(names(scenarios), run_scenario) |> list_rbind()

stopCluster(cl)
registerDoSEQ()

# 11 Evaluation ===============================================================

restrict_years <- function(oos, max_year) {
  oos |> mutate(preds = map(preds, \(d) filter(d, !is.na(.pred_1),
                                               year <= max_year)))
}



## 11.2 DeLong test: Logit vs. each ML model ==================================
# Paired DeLong test on the pooled out-of-sample predictions. Both models are
# scored on exactly the same rows (joined on iso3c and year).
# NOTE: DeLong assumes independent observations, which the panel violates
# (overlapping pre-crisis windows, common crisis waves). Treat p-values as
# optimistic.

delong_one <- function(p_logit, p_ml) {
  tryCatch({
    d <- inner_join(
      select(p_logit, iso3c, year, truth, p_logit = .pred_1),
      select(p_ml, iso3c, year, truth_ml = truth, p_ml = .pred_1),
      by = c("iso3c", "year")
    ) |> filter(!is.na(p_logit), !is.na(p_ml))
    stopifnot(all(d$truth == d$truth_ml))
    
    roc_l <- pROC::roc(d$truth, d$p_logit, levels = c("0", "1"),
                       direction = "<", quiet = TRUE)
    roc_m <- pROC::roc(d$truth, d$p_ml, levels = c("0", "1"),
                       direction = "<", quiet = TRUE)
    test  <- pROC::roc.test(roc_l, roc_m, method = "delong", paired = TRUE,
                            alternative = "two.sided")
    
    auc_l <- as.numeric(pROC::auc(roc_l))
    auc_m <- as.numeric(pROC::auc(roc_m))
    diff  <- auc_m - auc_l
    z     <- unname(test$statistic)
    se    <- if (is.finite(z) && z != 0) abs(diff / z) else NA_real_
    
    tibble(n_matched = nrow(d), auc_logit = auc_l, auc_ml = auc_m,
           diff = diff, se = se,
           ci_lower = diff - 1.96 * se, ci_upper = diff + 1.96 * se,
           z = z, p_value = test$p.value)
  }, error = function(e) {
    tibble(n_matched = NA_integer_, auc_logit = NA_real_, auc_ml = NA_real_,
           diff = NA_real_, se = NA_real_, ci_lower = NA_real_,
           ci_upper = NA_real_, z = NA_real_, p_value = NA_real_)
  })
}

# Holm adjustment within each scenario x horizon x specification
# (3 comparisons: Logit vs. RF, SVM, MLP). Change `group_by` to widen or
# narrow the family of tests.
delong_table <- function(oos, max_year = Inf) {
  oos <- restrict_years(oos, max_year)
  logit <- oos |>
    filter(model == "Logit") |>
    select(scenario, outcome, specification, p_logit = preds)
  
  oos |>
    filter(model != "Logit") |>
    inner_join(logit, by = c("scenario", "outcome", "specification")) |>
    mutate(res = map2(p_logit, preds, delong_one)) |>
    select(-preds, -p_logit) |>
    unnest(res) |>
    group_by(scenario, outcome, specification) |>
    mutate(p_holm = p.adjust(p_value, method = "holm")) |>
    ungroup()
}

## 11.3 ROC-AUC with standard DeLong confidence intervals =====================

delong_roc_ci <- function(preds, level = 0.95) {
  preds <- filter(preds, !is.na(.pred_1))
  roc_obj <- pROC::roc(preds$truth, preds$.pred_1,
                       levels = c("0", "1"), direction = "<", quiet = TRUE)
  ci <- pROC::ci.auc(roc_obj, conf.level = level, method = "delong")
  tibble(auc = as.numeric(pROC::auc(roc_obj)),
         ci_lower = as.numeric(ci[1]), ci_upper = as.numeric(ci[3]))
}

roc_ci_table <- function(oos, max_year = Inf) {
  restrict_years(oos, max_year) |>
    mutate(result = map(preds, delong_roc_ci)) |>
    select(-preds) |>
    unnest(result)
}


## 11.4 PR-AUC with country-cluster bootstrap intervals =======================

boot_pr_auc <- function(preds, n_boot = 1000, seed = 123, level = 0.95) {
  preds <- filter(preds, !is.na(.pred_1))
  stat <- function(d) {
    if (n_distinct(d$truth) < 2) return(NA_real_)
    suppressWarnings(tryCatch(
      pr_auc_vec(d$truth, d$.pred_1, event_level = "second"),
      error = \(e) NA_real_))
  }
  est <- stat(preds)
  by_cntry <- unname(split(preds, preds$iso3c))
  n_cntry <- length(by_cntry)
  set.seed(seed)
  boots <- map_dbl(seq_len(n_boot), \(i) {
    stat(list_rbind(by_cntry[sample.int(n_cntry, n_cntry, replace = TRUE)]))
  })
  boots <- boots[is.finite(boots)]
  if (length(boots) < 2) stop("Insufficient valid bootstrap replicates")
  tibble(auc = est, se = sd(boots),
         ci_lower = unname(quantile(boots, (1 - level) / 2)),
         ci_upper = unname(quantile(boots, 1 - (1 - level) / 2)),
         n_valid = length(boots), prevalence = mean(preds$truth == "1"))
}

pr_boot_table <- function(oos, max_year = Inf) {
  restrict_years(oos, max_year) |>
    mutate(n_pos = map_int(preds, \(p) sum(p$truth == "1", na.rm = TRUE)),
           n_pos_countries = map_int(preds, \(p) n_distinct(p$iso3c[p$truth == "1"])),
           boot = map(preds, \(p) {
             key <- rlang::hash(list(p, n_boot, deparse(boot_pr_auc)))
             cached(file.path(cache_dir, paste0("prboot_", substr(key, 1, 10), ".rds")),
                    boot_pr_auc(p, n_boot = n_boot))
           })) |>
    select(-preds) |>
    unnest(boot)
}

# 12 Tables ===================================================================

latex_auc_ci_table <- function(ab, metric, caption, label) {

  wide <- ab |>
    mutate(
      horizon = k_label(outcome),
      spec = str_to_title(specification),
      cell = sprintf("%.3f [%.3f, %.3f]",
                     auc, ci_lower, ci_upper)
    ) |>
    select(horizon, spec, model, cell,
           any_of(c("prevalence", "n_pos", "n_pos_countries"))) |>
    pivot_wider(names_from = model, values_from = cell)

  if (metric == "PR") {
    wide <- wide |>
      select(horizon, spec, all_of(models),
             prevalence, n_pos, n_pos_countries)
  } else {
    wide <- wide |>
      select(horizon, spec, all_of(models))
  }

  if (metric == "PR") {
    col_names <- c("Horizon", "Specification", models,
                   "Prevalence", "Pre-crisis obs.", "Countries")
    alignment <- "llccccrrr"
  } else {
    col_names <- c("Horizon", "Specification", models)
    alignment <- "llcccc"
  }

  wide |>
    kbl(
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      linesep = "",
      align = alignment,
      col.names = col_names,
      caption = caption,
      label = label
    ) |>
    collapse_rows(
      columns = 1,
      latex_hline = "major",
      valign = "top"
    )
}


latex_delong_table <- function(dl, caption, label) {
  dl |>
    mutate(
      horizon = k_label(outcome),
      spec    = str_to_title(specification),
      model   = factor(model, levels = ml_models),
      diff_s  = sprintf("%.3f%s", diff, stars(p_holm)),
      ci      = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
      p_holm_s = sprintf("%.3f", p_holm)
    ) |>
    arrange(outcome, specification, model) |>
    select(horizon, spec, model, auc_logit, auc_ml, diff_s, ci, p_holm_s) |>
    kbl(format = "latex", booktabs = TRUE, digits = 3, linesep = "",
        align = "lllrrrrr",
        col.names = c("Horizon", "Specification", "Model", "AUC Logit",
                      "AUC ML", "Difference", "95\\% CI", "p (Holm)"),
        escape = FALSE, caption = caption, label = label) |>
    collapse_rows(columns = 1:2, latex_hline = "major", valign = "top") |>
    footnote(
      general = paste(
        "Difference = AUC(ML) - AUC(Logit). Paired DeLong test, two-sided,",
        "on pooled out-of-sample predictions; the CI is the difference",
        "$\\pm$ 1.96 DeLong standard errors. Holm adjustment within each",
        "horizon and specification (three comparisons). Stars refer to the",
        "Holm-adjusted p-value: * p<0.1, ** p<0.05, *** p<0.01."
      ),
      general_title = "Notes: ", threeparttable = TRUE, escape = FALSE
    )
}

# Main results: standard DeLong ROC intervals; cluster-bootstrap PR intervals.
main_oos <- filter(oos_all, scenario == "main")
roc_ci_main  <- roc_ci_table(main_oos)
pr_boot_main <- pr_boot_table(main_oos)

write_table(latex_auc_ci_table(roc_ci_main, "ROC",
                               "Out-of-sample ROC-AUC with standard DeLong confidence intervals", "roc_ci_main"),
            "roc_ci_main")
write_table(latex_auc_ci_table(pr_boot_main, "PR",
                               "Out-of-sample PR-AUC with country-cluster bootstrap intervals", "pr_boot_main"),
            "pr_boot_main")

# 13 Robustness checks ========================================================
# All checks go through the same functions. RC4 re-evaluates main predictions.

checks <- tribble(
  ~check,                        ~scn,   ~max_year,
  "Main",                        "main", Inf,
  "RC1: reduced sample",         "rc1",  Inf,
  "RC2: real growth vars",       "rc2",  Inf,
  "RC3: region",                 "rc3",  Inf,
  "RC4: excluding COVID years",  "main", 2019
)

delong_checks <- pmap(checks, \(check, scn, max_year) {
  delong_table(filter(oos_all, scenario == scn), max_year) |>
    mutate(check = check)
}) |> list_rbind()

walk(which(checks$check != "Main"), \(i) {
  ck <- checks$check[i]
  scn <- checks$scn[i]
  max_year <- checks$max_year[i]
  tag <- tolower(sub(":.*", "", ck))
  oos <- filter(oos_all, scenario == scn)
  roc_ci <- roc_ci_table(oos, max_year)
  write_table(latex_auc_ci_table(roc_ci, "ROC",
                                 paste0("Out-of-sample ROC-AUC (", ck, ")"),
                                 paste0("roc_ci_", tag)), paste0("roc_ci_", tag))
  pr_boot <- pr_boot_table(oos, max_year)
  write_table(latex_auc_ci_table(pr_boot, "PR",
                                 paste0("Out-of-sample PR-AUC (", ck, ")"),
                                 paste0("pr_boot_", tag)), paste0("pr_boot_", tag))
  write_table(
    latex_delong_table(filter(delong_checks, check == ck),
                       paste0("Logit vs. machine learning, DeLong test (", ck, ")"),
                       paste0("delong_", tag)),
    paste0("delong_", tag)
  )
})

# 14 ROC and PR curves ========================================================
# Pooled out-of-sample predictions, one panel per horizon.

model_colours <- c(Logit = "black", RF = "#1b9e77", SVM = "#d95f02",
                   MLP = "#7570b3")

plot_roc_pr <- function(oos, scn_name = "main", spec_name = "full",
                        max_year = Inf) {
  d <- oos |>
    filter(scenario == scn_name, specification == spec_name) |>
    restrict_years(max_year) |>
    mutate(horizon = k_label(outcome),
           model   = factor(model, levels = models))
  
  roc_df <- d |>
    mutate(curve = map(preds, \(p) roc_curve(p, truth, .pred_1,
                                             event_level = "second"))) |>
    select(model, horizon, curve) |>
    unnest(curve) |>
    mutate(fpr = 1 - specificity)
  
  pr_df <- d |>
    mutate(curve = map(preds, \(p) pr_curve(p, truth, .pred_1,
                                            event_level = "second"))) |>
    select(model, horizon, curve) |>
    unnest(curve) |>
    filter(is.finite(precision))
  
  prev <- d |>
    mutate(prev = map_dbl(preds, \(p) mean(p$truth == "1"))) |>
    group_by(horizon) |>
    summarise(prev = first(prev), .groups = "drop")
  
  p_roc <- ggplot(roc_df, aes(fpr, sensitivity, colour = model, linetype = model)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_path(linewidth = 0.6) +
    facet_wrap(~horizon, nrow = 1) +
    scale_colour_manual(values = model_colours) +
    coord_equal() +
    labs(x = "False positive rate", y = "True positive rate",
         colour = NULL, linetype = NULL)
  
  p_pr <- ggplot(pr_df, aes(recall, precision, colour = model, linetype = model)) +
    geom_hline(data = prev, aes(yintercept = prev),
               colour = "grey50", linewidth = 0.4) +
    geom_path(linewidth = 0.6) +
    facet_wrap(~horizon, nrow = 1) +
    scale_colour_manual(values = model_colours) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Recall", y = "Precision", colour = NULL, linetype = NULL)
  
  (p_roc / p_pr) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# Main scenario, one figure per specification
walk(names(specs), \(s) {
  ggsave(file.path(fig_dir, paste0("roc_pr_main_", s, ".png")),
         plot_roc_pr(oos_all, "main", s), width = 11, height = 7)
})

# Robustness check 1
walk(c("baseline", "broad"), \(s) {
  ggsave(file.path(fig_dir, paste0("roc_pr_rc1_", s, ".png")),
         plot_roc_pr(oos_all, "rc1", s), width = 11, height = 7)
})

# Robustness checks 2 and 3
walk(c("rc2", "rc3"), \(r) {
  ggsave(file.path(fig_dir, paste0("roc_pr_", r, "_full.png")),
         plot_roc_pr(oos_all, r, "full"), width = 11, height = 7)
})

# Robustness check 4 (no COVID years)
walk(names(specs), \(s) {
  ggsave(file.path(fig_dir, paste0("roc_pr_rc4_", s, ".png")),
         plot_roc_pr(oos_all, "main", s, max_year = 2019), width = 11, height = 7)
})
