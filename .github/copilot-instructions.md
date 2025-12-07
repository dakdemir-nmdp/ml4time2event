# ml4time2event Copilot Instructions

These instructions describe how to work effectively inside the **ml4time2event** codebase using GitHub Copilot (Chat/Agents). They provide the context required to extend, test, and finalize the R package for survival and competing-risk modelling.

---

## Project Overview

ml4time2event is an R package that unifies statistical and machine-learning approaches for time-to-event modelling. It wraps task constructors, engine registries, ensembling, explainability helpers, and pipelines for both survival and competing-risk outcomes.

Primary goals:
- Construct tidy survival/competing-risk tasks from raw data.
- Train a registry-backed set of engines (Cox, glmnet, GBM, random forest, BART, xgboost, etc.).
- Provide lightweight ensembling, evaluation metrics, and explainability (SHAP, plots).
- Offer reproducible pipelines that integrate recipes, resampling, and prediction APIs.

The package is designed for programmatic use from R scripts, vignettes, and tidy workflows (rsample/recipes/tidymodels). It uses `renv` for dependency pinning.

---

## Architecture Overview

### Directory Highlights
- `R/`: primary source files grouped by concern.
  - `core_*`: core infrastructure (`core_tasks.R`, `core_fit.R`, `core_predictions.R`, `core_registry.R`, pipelines, R6 base classes).
  - `surv_*.R` & `cr_*.R`: model-specific wrappers for survival vs competing-risk engines plus metric utilities.
  - `ensemble_*`: ensembling helpers and persistence.
  - `data_*.R`, `vignette_data_utils.R`: dataset access/cleaning utilities.
  - `shap_*`, `plotting_functions.R`: explainability + autoplot logic.
  - `recipes_steps.R`: custom recipe steps that integrate with `recipes::step_*`.
- `tests/`: `testthat` suite and fixtures.
- `vignettes/`: longer tutorials (comprehensive survival/competing-risk workflows, SHAP explainability).
- `data/` + `data-raw/`: packaged datasets and scripts to build them.
- `renv/` + `renv.lock`: reproducible dependency environment.

### Key Components
- **Task constructors** (`ml4t2e_task_surv`, `ml4t2e_task_cr`, `ml4t2e_task_from_df`): validate inputs, store metadata, and prepare defaults (time/status/cause columns, factor levels, stratification info).
- **Model registry** (`ml4t2e_list_models`, `ml4t2e_register_model`): describes supported engines, required packages, default hyperparameters, and whether they operate on survival vs competing-risk data.
- **Fitting API** (`ml4t2e_fit`, `ml4t2e_fit_pipeline`): orchestrates preprocessing, model training, ensembling, and evaluation across tasks.
- **Pipeline API** (`ml4t2e_pipeline`, `pipeline_simple.R`, `pipeline_wrappers.R`): R6 objects that combine recipes, resampling, fitting, prediction, evaluation, persistence.
- **Explainability** (`ml4t2e_explain`, `shap_analysis.R`, `shap_plotting.R`, `autoplot.t2e_*`): compute SHAP-like importances, dependence plots, survival/CIF curves.
- **Metrics** (`surv_metrics.R`, `cr_metrics.R`): concordance, Brier/IBS, time-dependent metrics, CIF calibration.

---

## Core Usage Workflow

1. **Prepare Data**
   - Start with a tibble/data.frame containing predictors plus outcome columns.
   - Survival: `time`, `status` (0 = censored, 1 = event).
   - Competing risk: `ftime`, `status` (0 = censored, 1 = event), `cause` (factor/int cause for events).

2. **Construct Tasks**
   ```r
   surv_task <- ml4t2e_task_surv(lung, time = "time", event = "status")
   cr_task   <- ml4t2e_task_cr(bmt, time = "ftime", status = "status", cause = "cause")
   ```
   Tasks store metadata (censoring, metrics, resampling defaults).

3. **Fit Models / Pipelines**
   ```r
   surv_fit <- ml4t2e_fit(
     surv_task,
     models   = c("cox", "random_forest", "gbm"),
     ensemble = "auto",
     resampling = rsample::vfold_cv(surv_task$data, v = 3)
   )

   pipe <- ml4t2e_pipeline(
     outcome    = list(type = "survival", time = "time", event = "status"),
     models     = c("cox", "glmnet"),
     recipe     = recipes::recipe(~ ., data = lung) %>% recipes::step_impute_median(all_numeric_predictors()),
     ensemble   = "simple",
     resampling = rsample::vfold_cv(lung, v = 3)
   )
   pipe$fit(lung)
   ```

4. **Predict & Evaluate**
   ```r
   surv_preds <- predict(surv_fit, type = "survival")
   ml4t2e_evaluate(surv_fit)   # returns tibble of metrics
   autoplot(surv_fit)          # autoplot methods for survival curves/CIFs

   pipe$predict(new_data)      # uses stored recipe + model
   pipe$evaluate(holdout_data)
   pipe$save("artifacts/")
   ```

5. **Explainability**
   ```r
   expl <- ml4t2e_explain(surv_fit, new_data = lung[1:50, ], approach = "shap")
   autoplot(expl, type = "importance")
   ```

---

## Testing & Validation

- Tests live in `tests/testthat/`. Use `testthat` 3e syntax (`test_that`, `expect_*`). Run via `devtools::test()` or `testthat::test_local()`.
- Continuous integration should run `devtools::check()` (R CMD check) plus targeted scripts (`pkgdown`, `lintr`) as needed.
- When adding engines or metrics:
  - Create fixtures under `tests/fixtures` or inline data frames.
  - Add deterministic seeds via `set.seed()` or `withr::with_seed()`.
  - Keep runtime low (< a few seconds) by reducing resampling folds or dataset sizes inside tests.
- Example manual workflow:
  ```r
  renv::restore()
  devtools::document()
  devtools::test()
  devtools::check()
  pkgdown::build_site()
  ```

---

## Critical Patterns & Conventions

- **Roxygen2**: exported functions/classes must have roxygen blocks; regenerate docs via `devtools::document()` and commit resulting `man/` + `NAMESPACE`.
- **Naming**: exported helpers use the `ml4t2e_` prefix; internal helpers use camelCase or descriptive verbs within relevant files.
- **Dependencies**: prefer CRAN packages. Use `Suggests` for optional heavy engines; guard usage with `rlang::check_installed()`.
- **Error handling**: use `rlang::abort()`/`cli::cli_abort()` with informative messages; avoid silent failures.
- **Parallelism**: `future`, `foreach`, and base parallelization need deterministic seeds; document or enforce safe defaults in pipelines.
- **Data objects**: packaged datasets live in `data/` as `.rda`; update via scripts under `data-raw/`.
- **Environments**: renv locks versions; after adding packages call `renv::snapshot()`.

---

## Performance Considerations

- Heavy engines (GBM, random forest, BART, xgboost, SHAP explainers) can be costly. Provide options to reduce iterations, tree depth, or use smaller folds for quick tests/examples.
- Reuse recipes and intermediate matrices when possible; pipeline objects cache them.
- When benchmarking, capture `Sys.time()` or use `bench::mark()`; document dataset size and hardware.

---

## Documentation & Release

- README covers installation, quickstart (survival + competing risk), renv restore instructions, and how to run tests/checks.
- Vignettes provide end-to-end tutorials. Keep them reproducible; wrap expensive chunks in `if (interactive())` or `\donttest{}` as appropriate.
- Releases require:
  - Clean `devtools::check()` (0 ERROR/WARNING; NOTES justified).
  - Updated `DESCRIPTION` version + `NEWS.md`.
  - Synced docs (`devtools::document()`).
  - Optional pkgdown build & GitHub release notes.

---

## Quality Assurance & Verification

- **Model Results**: Review results from models to ensure they make sense, behave as expected, and are correct. Check that predictions align with domain expectations (e.g., survival probabilities decrease over time).
- **UI Consistency**: Ensure the user interface (API signatures, object structures, output formats) is consistent across the package.
- **Documentation**: Verify that documentation is correct and complete. Ensure examples run and descriptions match the code behavior.

---

## Working Effectively with Copilot

- Provide context: file path, function, desired behavior, data shapes (time/status/cause), resampling setup.
- Mention whether you need survival or competing-risk behavior; specify metrics (C-index, Brier, CIF) and acceptable runtime.
- For tests, specify `testthat` expectations and datasets. For docs, supply target audience and examples to include.
- When editing pipelines/registries, note which engines should remain stable for backwards compatibility.

These guidelines should give Copilot enough context to help you clean up, extend, and finalize ml4time2event without regressing its survival/competing-risk workflows.
