# Predicting Systemic Banking Crises: Logit vs. Machine Learning Early-Warning Models

Replication code for my master's thesis *"Predicting Banking Crises: A Comparative Study of Logit Models and Machine Learning Approaches"* at the Institute of Statistics and Econometrics, Faculty of Business, Economics and Social Sciences, Kiel University (2026). Supervisor: Prof. Dr. Kai Carstensen.

Please note: The master's thesis is still in progress. It will be completed by December 1, 2026. 

## Overview

The thesis asks whether machine learning models predict systemic banking crises better out of sample than the standard logit early-warning model. Using a country-year panel from 1970 to 2025, it compares:

| Model | Implementation |
|-------|----------------|
| Logit | `glm` (benchmark) |
| Random Forest (RF) | `ranger` via tidymodels |
| Support Vector Machine (SVM, RBF kernel) | `kernlab` via tidymodels |
| Multilayer Perceptron (MLP) | `brulee` via tidymodels |

**Outcome.** Binary pre-crisis indicators `precrisis1` to `precrisis4`: a systemic banking crisis starts within the next K = 1, 2, 3 or 4 years. Crisis dates come from the Laeven & Valencia Systemic Banking Crises Database (2026 update). Crisis years themselves are excluded from estimation.

**Predictors.** Three specifications of increasing size (`baseline`, `broad`, `full`), covering credit, government debt, GDP growth, inflation, loans-to-deposits, current account, broad money, reserves, stock returns and an advanced-economy dummy.

**Evaluation design**
- Pseudo-real-time expanding-window forecasts, with test years from 2005 on.
- A gap of K years between the training window and the forecast year, so that only labels known at the forecast date are used.
- Hyperparameters are re-tuned every 5 years (2005, 2010, 2015, 2020) using only information available at that date.
- ROC-AUC is the main metric. PR-AUC is reported as a supplement with country-cluster bootstrap confidence intervals.
- Logit and ML models are compared with paired DeLong tests and Holm-adjusted p-values.
- Robustness checks:
  - RC1: reduced common sample
  - RC2: real growth variables
  - RC3: region instead of advanced-economy dummy
  - RC4: excluding COVID years

## Repository structure

```
.
├── scripts/
│   ├── 01a_clean_data.R             # Read raw files, convert to long format, add ISO codes
│   ├── 01b_construct_indicators.R   # Combine sources into harmonised indicators
│   ├── 01c_create_panel.R           # Build the final country-year panel and outcomes
│   ├── 02_descriptive_analysis.R    # Crisis counts, event-study plots, summary statistics
│   └── 03_estimation_and_prediction.R  # Tuning, forecasting, evaluation, tables, figures
├── data/
│   ├── raw/        # Original downloads (see "Data sources")
│   ├── interim/    # Cleaned datasets and constructed indicators (.rds)
│   └── final/      # panel.rds, the estimation dataset
├── cache/          # Cached tuning results, out-of-sample predictions, bootstrap results
├── output/
│   ├── figures/    # ROC and PR curves (.png)
│   └── tables/     # LaTeX tables (.tex)
├── renv/, renv.lock, .Rprofile      # Reproducible package environment
└── Predicting_Banking_Crises.Rproj
```

## Data sources

| Source | Content |
|--------|---------|
| IMF Systemic Banking Crises Database (Laeven & Valencia, 2026) | Crisis start and end years |
| BIS | Credit, credit-to-GDP, CPI, property prices |
| IMF: Global Debt Database, World Economic Outlook, Regional Economic Outlook for Africa, Public Finances in Modern History, Monetary and Financial Statistics, National Economic Accounts | Debt, GDP, inflation, current account, monetary aggregates, interest rates, share prices |
| OECD | Share prices, interest rates, property prices, exchange rates |
| Eurostat | Short-term interest rates |
| World Bank: WDI, Global Financial Development | Credit, broad money, reserves, loans-to-deposits, stock returns |
| Jordà-Schularick-Taylor Macrohistory Database (R6) | Long-run interest rates, money, loans-to-deposits |

Where several sources cover the same variable, `01b_construct_indicators.R` chooses the source with the most observations per country as the primary series. It then splices in secondary sources that correlate with it sufficiently (min. correlation 0.8 over at least 8 overlapping years) after rescaling. A few documented manual corrections for implausible values are applied there as well.

> **Note on raw data:** The raw files are included in `data/raw/`.

## Reproducing the results

**Requirements:** R 4.6.0 and the packages pinned in `renv.lock`. The ML tuning is computationally heavy and uses all but one CPU core.

```r
# 1. Restore the package environment
renv::restore()

# 2. Build the dataset (run from the project root)
source("scripts/01a_clean_data.R")
source("scripts/01b_construct_indicators.R")
source("scripts/01c_create_panel.R")

# 3. Descriptive analysis
source("scripts/02_descriptive_analysis.R")

# 4. Estimation, forecasting and evaluation
source("scripts/03_estimation_and_prediction.R")
```

Open `Predicting_Banking_Crises.Rproj` first so that relative paths resolve correctly.

**Caching.** Tuning, forecasting and bootstrap steps in script 03 are stored in `cache/`. Cache keys hash the data, configuration and model code, so re-running an unchanged script is nearly instant, and changing one model's workflow only recomputes that model. If you change helper functions that are not tracked in the key (e.g. `get_panel`), delete `cache/` to force a full recomputation.

## Outputs

- `output/tables/*.tex`: ROC-AUC and PR-AUC tables with confidence intervals, and DeLong test tables for the main analysis and each robustness check.
- `output/figures/*.png`: pooled out-of-sample ROC and precision-recall curves for each horizon and specification.

## Citation

If you use this code, please cite:

> Burak Bozkurt (2026). *Predicting Banking Crises: A Comparative Study of Logit Models and Machine Learning Approaches*. Master's thesis, Kiel University.

## License

The code in this repository is released under the [MIT License](LICENSE).
The datasets in `data/raw/` remain subject to the terms of their original
providers (BIS, IMF, OECD, Eurostat, World Bank, Jordà-Schularick-Taylor)
and are not covered by this license.
