# ml4time2event

Machine learning for time-to-event analysis. Provides tools for predicting survival outcomes and competing risks using statistical and machine learning methods.

## Installation

```r
# Install from GitHub with all runtime dependencies
install.packages("pak")  # or remotes
pak::pak("dakdemir-nmdp/ml4time2event")

# Alternatively with renv (recreates the development library)
install.packages("renv")
renv::restore()
```

## Reproducible Development Environment

This repo uses [`renv`](https://rstudio.github.io/renv) to lock dependency versions. To recreate the exact package set used during development:

```r
install.packages("renv")
renv::restore()
```

Packages install into the project-local `renv/library` directory (configured via `.Renviron`) so the environment stays isolated from your global R library. Some dependencies, such as `igraph`, may require external toolchains (e.g., `gfortran`) when binaries are unavailable on your platform; install those system requirements before running `renv::restore()`.

When you change dependencies, run `renv::snapshot()` to update `renv.lock` so teammates get the same versions.

## Pipeline Quickstart

```r
library(ml4time2event)
library(dplyr)

# Survival pipeline -------------------------------------------------------
lung_df <- get_lung_survival_data()

surv_pipeline <- ml4t2e_fit_pipeline(
  data = lung_df,
  analysis_type = "survival",
  timevar = "time",
  eventvar = "status",
  models = c("glmnet", "coxph"),
  include_rf = FALSE,
  prediction_times = seq(0, 1000, length.out = 50)
)

surv_preds <- predict(
  surv_pipeline,
  newdata = lung_df[1:5, ],
  new_times = seq(0, 730, length.out = 25)
)

surv_preds$predictions$models_used

# Persist the trained pipeline (optional)
tmp_path <- tempfile("lung_pipeline_", fileext = ".rds")
ml4t2e_save_pipeline(surv_pipeline, tmp_path)
restored_pipeline <- ml4t2e_load_pipeline(tmp_path)

# Competing-risks pipeline -----------------------------------------------
bmt_df <- get_bmt_competing_risks_data()

cr_pipeline <- ml4t2e_fit_pipeline(
  data = bmt_df,
  analysis_type = "competing_risks",
  timevar = "ftime",
  eventvar = "status",
  models = c("FG", "cox"),
  include_rf = FALSE,
  prediction_times = seq(0, 150, length.out = 40)
)

cr_preds <- predict(
  cr_pipeline,
  newdata = bmt_df[1:4, ],
  ensemble_method = "average"
)

cr_preds$predictions$models_used

# SHAP-based explainability -----------------------------------------------
# Explain predictions using SHAP values
shap_result <- ml4t2e_calculate_shap(
  pipeline = surv_pipeline,
  data = lung_df[1:50, ],
  time_horizon = 365,  # 1-year expected time lost
  nsim = 100
)

# Variable importance plot
ml4t2e_shap_importance(shap_result)

# Dependence plot showing feature effects
ml4t2e_shap_dependence(shap_result, feature = "age")

# Explain individual prediction
ml4t2e_shap_waterfall(shap_result, obs_id = 1)
```

The pipelines combine preprocessing, model fitting, ensemble construction, prediction,
persistence, and explainability in a single object that can be reloaded and used for production scoring.

## Features

- **Pipeline API**: `ml4t2e_fit_pipeline()` creates end-to-end survival or competing-risks workflows with preprocessing, modeling, and persistence.
- **Survival Models**: Cox, Random Forest, XGBoost, GAM, BART, DeepSurv, GLMNet, GBM, RuleFit
- **Competing Risks**: Fine-Gray, cause-specific Cox, and ML models for competing risks
- **Ensemble Methods**: Averaging, weighted averaging, super learner stacking
- **Metrics**: C-index, Brier score, integrated metrics, expected time lost
- **Explainability**: SHAP-based variable importance and dependence analysis for interpretable predictions
- **Comprehensive Testing**: 1900+ passing tests

## Supported Models

- **Survival (RunSurvModels `models` argument)**: `glmnet`, `coxph`, `rulefit`, `xgboost`, `gam`, `gbm`, `ExpSurvReg`, `WeibSurvReg`, `bart`, `deepsurv` (plus two Random Forest baselines when `include_rf = TRUE`).
- **Competing Risks (RunCRModels `models` argument)**: `FG`, `rulefit`, `bart`, `cox`, `xgboost`, `gam`, `survreg` (plus two Random Forest baselines when `include_rf = TRUE`).

## Vignettes

- `vignettes/comprehensive_survival_analysis.R` - Complete survival analysis workflow
- `vignettes/comprehensive_competing_risks_analysis.R` - Competing risks analysis
- `vignettes/shap_explainability.Rmd` - SHAP-based explainability and interpretation

## License

GPL (>= 2)

## Author

Deniz Akdemir (dakdemir@nmdp.org)
