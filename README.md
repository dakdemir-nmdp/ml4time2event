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

## Quickstart (3–5 lines)

```r
library(ml4time2event)
library(dplyr)

# Survival ---------------------------------------------------------------
lung <- get_lung_survival_data() %>% mutate(status = if_else(status == 2, 1L, 0L))
surv_task <- ml4t2e_task_surv(lung, time = "time", event = "status")
surv_fit  <- ml4t2e_fit(surv_task, models = c("cox", "random_forest"), ensemble = "auto")
surv_preds <- predict(surv_fit, type = "survival")
ml4t2e_evaluate(surv_fit)
autoplot(surv_fit)

# Competing risks --------------------------------------------------------
bmt <- get_bmt_competing_risks_data() %>%
  mutate(
    cause = if_else(status == 0L, NA_integer_, status),
    status = if_else(status == 0L, 0L, 1L)
  )
cr_task <- ml4t2e_task_cr(bmt, time = "ftime", status = "status", cause = "cause")
cr_fit  <- ml4t2e_fit(cr_task, models = c("cox", "fine_gray"), ensemble = "auto")
cr_preds <- predict(cr_fit, type = "cif")
ml4t2e_evaluate(cr_fit)
autoplot(cr_fit)
```

The `ml4t2e_task_*()` constructors standardise raw datasets. `ml4t2e_fit()` drives the registry-backed engines (now returned as tidy `t2e_fit` objects), while `predict()`, `ml4t2e_evaluate()`, and `autoplot()` all work with the same tidy outputs for both survival and competing-risk workflows.

### Pipelines, recipes, and resampling

```r
library(rsample)
library(recipes)

pipe <- ml4t2e_pipeline(
  outcome  = list(type = "survival", time = "time", event = "status"),
  models   = c("cox", "gbm"),
  ensemble = "auto",
  recipe   = recipe(~ ., data = lung) %>%
    step_impute_median(all_numeric_predictors()),
  resampling = vfold_cv(lung, v = 3)
)

pipe$fit(lung)
pipe$summary()
pipe_preds <- pipe$predict(lung[1:5, ], times = seq(0, 365, length.out = 25))
autoplot(pipe_preds)
```

The pipeline object encapsulates preprocessing (via `recipes`), modelling (`ml4t2e_fit()`), optional resampling, and prediction/evaluation helpers. The fitted pipeline stores the trained model, the prepped recipe, training metrics, and (when requested) cross-validation summaries. Use `$predict()` with new data to bake the recipe and score in one step; `$evaluate()` applies the same processing when you provide hold-out data.

## Features

- Task constructors `ml4t2e_task_surv()` / `ml4t2e_task_cr()` create tidy `t2e_task` objects with metadata and defaults.
- Registry-backed engine discovery via `ml4t2e_list_models()` with package requirements and descriptive tags.
- Unified training and inference: `ml4t2e_fit()` → `predict()` → `ml4t2e_evaluate()` share tidy `t2e_*` outputs for survival and competing-risk paths.
- Lightweight ensembling: enable `ensemble = "auto"` or `"simple"` for averaged predictions recorded alongside per-model scores.
- Pipeline API: `ml4t2e_pipeline()` wraps recipes, resampling, fitting, and persistence in a single R6 object.
- Compatibility helper `ml4t2e_fit_pipeline()` maps the legacy signature onto the new pipeline engine for quick migrations.
- Explainability helpers: `ml4t2e_explain()` + `autoplot()` produce SHAP-style importance, dependence, and curve diagnostics.

## Supported Models

- Survival engines (via `ml4t2e_list_models("survival")`): `cox`, `glmnet`, `random_forest`, `gbm`, `bart`, `gam`, `survreg`, `rulefit`, `shallownn`, `ttah`, `xgboost`.
- Competing-risk engines (via `ml4t2e_list_models("competing_risk")`): `cox`, `fine_gray`, `cr_random_forest`, `cr_bart`, `cr_gam`, `cr_rulefit`, `cr_survreg`, `cr_shallownn`, `cr_ttah`, `cr_xgboost`.


## Vignettes

- `vignettes/comprehensive_survival_analysis.R` - Complete survival analysis workflow
- `vignettes/comprehensive_competing_risks_analysis.R` - Competing risks analysis
- `vignettes/shap_explainability.Rmd` - SHAP-based explainability and interpretation

## License

GPL (>= 2)

## Author

Deniz Akdemir (dakdemir@nmdp.org)
