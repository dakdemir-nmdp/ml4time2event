# Object-Oriented Refactor Plan for `ml4time2event`

## Objectives
- Introduce a coherent object-oriented layer (R6-based) for survival and competing risk (CR) models, ensemblers, and pipelines while keeping user-facing APIs stable.
- Restructure the `R/` directory into logical subfolders to clarify responsibilities and ease navigation.
- Reduce the number of exported functions so that only the primary user workflows remain public; keep helpers and implementation details internal.

## Clean‑Slate Policy (No Backward Compatibility)
- Treat this as a new, first-version API with no compatibility shims.
- Legacy function names and list-based returns are retired in favor of a simple, object-driven, tidy-output interface.
- NAMESPACE will export only the new high-level verbs and S3 methods described below.

## User Experience Principles
- Make the common path trivial: a 3–5 line “fit → predict → plot” workflow for both survival and CR.
- Predictable, type-stable returns: predictions and metrics always return tidy frames with consistent columns.
- Sensible defaults with few knobs: good defaults for time grids, metrics, and plots; advanced tuning available but tucked away.
- Progressive disclosure: export high-level verbs first; power users can opt into objects and registries.
- Clear, actionable messages: friendly errors and guidance toward the easiest next step.
- No backward compatibility: legacy function names are retired; users adopt the new API.

## Current State Observations
- All model implementations are independent functions returning lists, leading to duplicated validation, inconsistent outputs, and ad-hoc method dispatch.
- Ensembling, pipelines, and SHAP utilities interact with model outputs via loosely defined list structures, making it fragile to extend or substitute components.
- `R/` is a flat directory containing ~40 files with mixed responsibilities (model fitting, preprocessing, plotting, utilities), and `NAMESPACE` exports nearly every function.
- Documentation (vignettes, reference docs) assumes the current functional API, so compatibility and messaging must be managed carefully.

## Architectural Direction
- **Object system**: Adopt the `R6` package (already lightweight; add to `Imports`) to provide encapsulated model objects with clear lifecycles. Expose S3 generics (`predict`, `print`, `summary`, `autoplot`) that delegate to underlying R6 instances to preserve base R ergonomics.
- **Core abstractions**:
  - `TimeToEventModel` (abstract R6) with required methods: `initialize(spec)`, `fit(task)`, `predict_survival(newdata, times)`, `predict_risk(newdata, times)`, `model_info()`, and optional `tune()`.
  - `SurvivalModel` and `CompetingRiskModel` subclasses, each implementing engine-specific logic (Cox, BART, GBM, etc.) and sharing validation helpers.
  - `ModelRegistry` module to register model classes with metadata (tags, supported outcomes, required packages) used by pipelines and vignettes.
  - `EnsemblerBase` (R6) with concrete `SurvivalEnsembler` and `CREnsemble`. Handles member model management, weight estimation, and prediction orchestration.
  - `Pipeline` object encapsulating data preparation (`recipes`), model fitting, ensembling, and persistence; integrates with the high-level verbs.
- **Clean API only**: No legacy wrappers or bridging layers. All users interact via the new task/fit/predict/plot verbs and S3 methods.

## Minimal, Friendly Public API
- High-level verbs (exported):
  - `ml4t2e_task_surv(data, time, event, features = NULL, id = NULL, time_units = NULL)`
  - `ml4t2e_task_cr(data, time, status, cause, features = NULL, id = NULL, time_units = NULL)`
  - `ml4t2e_fit(task, models = c("cox", "gbm"), ensemble = c(FALSE, TRUE, "auto", "stack", "simple"), recipe = minimal_data_recipe(), resampling = NULL, seed = NULL, controls = list())`
  - `predict(object, newdata = NULL, times = NULL, type = c("survival", "risk", "cif"), include = c("ensemble", "all", "best", "model:cox", "model:gbm"), ...)` where `object` is a `t2e_fit`.
  - `ml4t2e_evaluate(fit_or_preds, task = NULL, metrics = c("ibs", "c_index"), include = c("ensemble", "all"))`
  - `autoplot(x, what = c("curves", "calibration", "roc", "shap", "importance", "dependence", "waterfall"), ...)` where `x` is a `t2e_pred`, `t2e_fit`, or explanation tibble.
  - `ml4t2e_explain(fit, newdata = NULL, method = c("auto", "shap", "permutation"), include = c("ensemble", "best", "all", "model:cox"), times = NULL, top_n = 20)`
  - `ml4t2e_save(fit, path)` and `ml4t2e_load(path)` for persistence.
  - `ml4t2e_list_models(outcome = c("survival", "cr"))` to discover available engines and required packages.

- Type-stable predictions (tibble return):
  - Survival: columns `id`, `time`, `surv`, `model`, `ensemble` (logical), `set` ("train"/"test").
  - CR: columns `id`, `time`, `cause`, `cif`, `model`, `ensemble`, `set`.
  - Risk/LP: columns `id`, `risk` or `lp`, optional `time`.

- Print and summary:
  - `print(fit)` shows models used, data summary, times grid, metrics observed on train.
  - `summary(fit)` returns a list/tibble with per-model metrics and weights (if ensemble).

- Minimal quickstart (survival):
```
task <- ml4t2e_task_surv(lung, time = "time", event = "status")
fit  <- ml4t2e_fit(task, models = c("cox", "gbm"), ensemble = "auto")
pred <- predict(fit, type = "survival")
autoplot(pred, what = "curves")
```

- Minimal quickstart (CR):
```
task <- ml4t2e_task_cr(bmt, time = "ftime", status = "fstatus", cause = "cause")
fit  <- ml4t2e_fit(task, models = c("fg", "rf"), ensemble = TRUE)
pred <- predict(fit, type = "cif")
autoplot(pred, what = "curves")
```

## Directory Restructuring Proposal
```
R/
  core/                # registries, shared constants, zzz setup
  utils/               # general + math utilities, validation helpers
  data/                # data ingestion & recipes
  tasks/               # task definitions (survival vs CR datasets)
  models/
    survival/          # Cox, BART, GBM, glmnet, etc.
    competing_risk/    # Fine-Gray, CR-Cox, etc.
  ensemble/            # ensembling strategies, stacking, weighting
  pipelines/           # pipeline orchestration & persistence
  metrics/             # performance metrics (surv + CR)
  plotting/            # visualization helpers
  shap/                # SHAP analysis & plotting
```
- Convert existing `.R` files into these subfolders, grouping helpers with their objects. Update `DESCRIPTION` (`Collate` not required when using roxygen2) and ensure `devtools::load_all()` still sources recursively.

### File Mapping (Current → Target)
- Models (survival):
  - `R/surv_cox.R` → `R/models/survival/cox.R`
  - `R/surv_bart.R` → `R/models/survival/bart.R`
  - `R/surv_gbm.R` → `R/models/survival/gbm.R`
  - `R/surv_glmnet.R` → `R/models/survival/glmnet.R`
  - `R/surv_random_forest.R` → `R/models/survival/rfsrc.R`
  - `R/surv_rulefit.R` → `R/models/survival/rulefit.R`
  - `R/surv_gam.R` → `R/models/survival/gam.R`
  - `R/surv_reg.R` → `R/models/survival/survreg.R`
  - `R/surv_xgboost.R` → `R/models/survival/xgboost.R`
- Models (CR):
  - `R/cr_cox.R` → `R/models/competing_risk/cs_cox.R`
  - `R/cr_fine_gray.R` → `R/models/competing_risk/fine_gray.R`
  - `R/cr_random_forest.R` → `R/models/competing_risk/rfsrc.R`
  - `R/cr_xgboost.R` → `R/models/competing_risk/xgboost.R`
  - `R/cr_gam.R`, `R/cr_survreg.R`, `R/cr_bart.R`, `R/cr_rulefit.R`, `R/cr_shallownn.R` → `R/models/competing_risk/...`
- Ensembling:
  - `R/surv_ensemble.R`, `R/cr_ensemble.R`, `R/ensemble_methods.R`, `R/ensemble_persistence.R` → `R/ensemble/`
- Pipelines:
  - `R/pipeline_wrappers.R`, `R/outcome_prediction.R`, `R/run_simulation.R` → `R/pipelines/`
- Metrics:
  - `R/surv_metrics.R`, `R/cr_metrics.R`, `R/model_evaluation.R` → `R/metrics/`
- Data & recipes:
  - `R/data_io.R`, `R/data_cleaning.R`, `R/data_transform.R`, `R/recipes_steps.R` → `R/data/`
- Plotting & SHAP:
  - `R/plotting_functions.R`, `R/shap_plotting.R` → `R/plotting/`
  - `R/shap_analysis.R` → `R/shap/analysis.R`
- Utilities:
  - `R/general_utils.R`, `R/math_utils.R`, `R/ttah_utils.R` → `R/utils/`
  - `R/survival_target_builder.R`, `R/vignette_data_utils.R`, `R/data_summary.R` → `R/utils/` or `R/tasks/` depending on usage.

## Workstreams & Key Steps
1. **Discovery & Preparation**
   - Inventory current functions, their dependencies, and exported API (e.g., via `PACKAGE_FUNCTIONS_REFERENCE.md`, `NAMESPACE`).
   - Tag critical user workflows (model training, prediction, pipelines, ensembling) that must remain stable.
   - Set up regression test harness (extend `tests/` as needed) to capture current behavior before refactor.
   - Decide on optional dependencies policy per engine (e.g., xgboost, randomForestSRC) and provide graceful errors with install hints.
   - Define target classes and S3 surface: `t2e_task`, `t2e_fit`, `t2e_pred`, `t2e_ensemble` and their S3 methods.
   - Define explanation artifact shape (`t2e_explain`) and persistence format (`.rds` with `api_version`).

2. **Foundation & Utilities**
   - Introduce `R6` to `Imports` and add a `core/r6_base.R` file defining abstract classes and mixins (validation, serialization).
   - Centralize common validation, preprocessing, and profiling utilities into `utils/`, replacing duplicated logic inside model functions.
   - Define a lightweight `Task` object (`SurvivalTask`, `CRTask`) capturing data, formula, and metadata to standardize inputs.
   - Add friendly messaging utilities (consider `cli`, `glue`) for errors and informative messages.
   - Establish a `times` policy: default to pretty sequence over observed event times (e.g., 100 quantiles), allow explicit override, and store `time_units`.
   - Define ensemble weight estimation defaults (stacking on IBS for survival, iIBS or log-loss equivalent for CR) with pluggable objective.

3. **Model Refactor (Survival)**
   - For each survival engine, create a corresponding R6 class under `models/survival/` implementing the `SurvivalModel` interface.
   - Ensure all classes return consistent prediction objects (`SurvivalPrediction` S3 class) with survival curves, risk scores, and meta info.
   - No legacy wrappers; implement the new classes and S3 surfaces only.
   - Implement `predict_survival()`, `predict_risk()`, and optional `predict_time(q = 0.5)` for median survival time.
   - Support engine-specific `controls` routed via the registry; validate and document common knobs.

4. **Model Refactor (Competing Risk)**
   - Mirror the survival approach in `models/competing_risk/`, aligning method names (`predict_cif`, `predict_event_time`), and shared utilities (cause-specific handling).
   - Update CR metrics and plotting functions to consume standardized prediction objects.
   - Normalize cause labels; store a `cause_map` between original labels and integer codes in the task.
   - Define `event_of_interest` handling and filtering for evaluation and plotting APIs.

5. **Ensembler & Pipeline Overhaul**
   - Redesign `SurvEnsemble` and `CREnsemble` as `EnsemblerBase` subclasses; consolidate weighting logic and persistence in `ensemble/`.
   - Refactor `pipeline_wrappers.R` into `pipelines/` with an R6 `Pipeline` class orchestrating tasks, models, ensemblers, and SHAP.
   - Ensure serialization works with new objects via `saveRDS` + class versioning; optionally provide minimal save/load helpers.
   - Provide simple cross-validation and resampling helpers; expose `ml4t2e_fit(..., resampling = vfolds(5))` with consistent summaries.
   - Ensembling UX: `ensemble = FALSE | TRUE | "auto" | "stack" | "simple"`; advanced options under `controls$ensemble`.

6. **API & Namespace Curation**
   - Define the public API surface (high-level tasks and fit functions, discovery, and S3 methods). Keep it minimal and consistent.
   - Remove `@export` from internal helpers; add `@keywords internal` where appropriate. Regenerate `NAMESPACE`.
   - Export only: `ml4t2e_task_surv`, `ml4t2e_task_cr`, `ml4t2e_fit`, `ml4t2e_list_models`, `ml4t2e_evaluate`, `ml4t2e_explain`, `ml4t2e_save`, `ml4t2e_load`, plus S3 methods for `predict.t2e_fit`, `autoplot.t2e_pred`, `autoplot.t2e_fit`, `autoplot.t2e_explain`, `print.t2e_fit`, `summary.t2e_fit`, `print.t2e_ensemble`, `summary.t2e_ensemble`.
   - Keep internal: individual engine classes/fitters, registries, low-level validators, metrics implementations, SHAP internals.

7. **Documentation & Vignettes**
   - Update roxygen documentation for new classes and high-level verbs; ensure `@examples` use the object workflow.
   - Revise vignettes (e.g., `vignettes/comprehensive_survival_analysis.Rmd`) to demonstrate the new objects and the simple end-to-end flow.
   - Refresh package reference docs (`pkgdown`, if used) and README to describe the new architecture and quickstart.
   - Add a concise “Quickstart” vignette centered on the high-level verbs and common patterns.
   - Include a “Choosing Models” vignette listing engines, trade-offs, and required packages returned by `ml4t2e_list_models()`.

8. **Testing & Validation**
   - Expand unit tests to cover each model class, ensemble, and pipeline scenario with snapshot-based prediction checks.
   - Add integration tests for pipelines and SHAP.
   - Run CRAN checks (`devtools::check()`) iteratively to catch S3/S4 method issues and namespace regressions.
   - Add golden-tests for printing and plotting outputs (e.g., vdiffr for ggplot) to preserve UX.
   - Confirm type stability of returns with `checkmate` assertions in tests.

9. **Polish & Release**
   - Finalize object and prediction class documentation; ensure examples run quickly.
   - Prepare README and a “Quickstart” vignette focusing on the 3–5 line workflow.
   - Tag a clean initial release (e.g., `0.1.0`) with the new API and structure.

## Export Strategy Guidelines
- Export only high-level tasks, fit, model discovery, and S3 methods:
  - Functions: `ml4t2e_task_surv`, `ml4t2e_task_cr`, `ml4t2e_fit`, `ml4t2e_list_models`, `ml4t2e_evaluate`.
  - S3: `predict.t2e_fit`, `autoplot.t2e_pred`, `autoplot.t2e_fit`, `print.t2e_fit`, `summary.t2e_fit`, `print.t2e_ensemble`, `summary.t2e_ensemble`.
- Keep internal: engine classes, registries, validators, metrics implementations, SHAP internals, plotting helpers; annotate with `@keywords internal`.
- Regenerate `NAMESPACE` via roxygen; avoid exporting engine-specific or low-level utilities.
- Add contributor guidelines in `design/architecture.md` describing how to add a new engine with minimal surface area exposure.

## Risks & Mitigations
- **Surface creep**: Resist adding engine-specific exports; keep extensions behind the registry and `controls`.
- **Complexity creep**: Keep R6 classes thin; share behavior through mixins/utility functions; document class responsibilities clearly.
- **Learning curve for contributors**: Provide developer notes (`design/architecture.md`) describing object lifecycle, registry usage, and folder conventions.
- **Package build issues**: Validate that subfolder structure is compatible with `roxygen2` and R CMD build; add checks in CI.
- **Optional dependency churn**: At first call to an engine, verify presence and version; provide informative install message with `install.packages()` hint.
- **Performance regressions**: Benchmark key paths (Cox, GBM, RF, FG); store reference timings and memory footprints.

## Suggested Sequence & Milestones
1. Week 1: Discovery, architecture skeleton, utilities consolidation.
2. Weeks 2-3: Refactor survival engines (R6 + S3); adjust tests.
3. Week 4: Complete CR model refactor and ensemble overhaul.
4. Week 5: Pipeline refactor, serialization, documentation updates.
5. Week 6: Namespace pruning, vignette rewrites, release preparation.

## Acceptance Criteria per Milestone
- Skeleton landed: `R/core/`, `R/models/survival/cox.R` example class, `ml4t2e_task_surv`, `ml4t2e_fit` prototype; passes checks.
- Survival parity: Cox + GBM + glmnet engines refactored; `predict()` and `ml4t2e_evaluate()` stable with snapshots.
- CR parity: Fine-Gray + RF engines refactored; consistent CIF predictions and plots; ensemble weights reproducible.
- Pipeline parity: End-to-end quickstarts in vignettes run and render without warnings.
- Namespace curated: Only high-level verbs exported; internal helpers documented as internal; `NAMESPACE` regenerated.

## Additional UX Details
- Defaults:
  - `times`: 100 equally spaced points between min observed and 95th percentile of event/censoring times.
  - Metrics: survival uses `integratedBrier` and `integratedC`; CR uses `integratedBrierCR` and `integratedConcordanceCR`.
  - Plots: `autoplot(..., what = "curves")` shows mean curves with 95% ribbons when resampling is used.
- Messages:
  - Use `cli::cli_abort`/`cli::cli_warn` with guidance and examples; no deprecation layer.
- Reproducibility:
  - `ml4t2e_fit(..., seed = 42)` sets seeds across engines where applicable; store RNG kind in object metadata.
- Tidy integration:
  - All predictions return tibbles; plotting accepts those tibbles directly; `autoplot` builds ggplot2 objects.
- Optional advanced knobs:
  - Expose `controls = list()` in `ml4t2e_fit` to pass engine-specific tuning in a structured way, but keep defaults performant.

## Single vs Multi‑Model Flow
- Single model: `ml4t2e_fit(task, models = "cox", ensemble = FALSE)`
  - Predict: `predict(fit, type = "survival")` returns only the Cox predictions.
  - Evaluate: `ml4t2e_evaluate(fit)` yields a one-row metrics table.
- Multiple models: `ml4t2e_fit(task, models = c("cox","gbm","glmnet"), ensemble = "auto")`
  - Predict: `predict(fit, include = "all")` returns stacked predictions (one row per id-time-model) plus an `ensemble = TRUE/FALSE` flag.
  - Filtered predict: `predict(fit, include = c("ensemble","model:cox"))` for only ensemble and a specific model.
  - Evaluate: `ml4t2e_evaluate(fit, include = "all")` returns per-model and ensemble metrics.

## Ensembling UX
- `ensemble` argument controls whether and how ensembling is done:
  - `FALSE`: no ensemble; keep per-model fits.
  - `TRUE`/`"auto"`: stacking with cross-validated weights optimized for the default metric family.
  - `"stack"`: explicit stacking; `controls$ensemble` can set objective and regularization.
  - `"simple"`: simple average or median; choose via `controls$ensemble$rule = "mean"|"median"`.
- Weights and method recorded in `summary(fit)`; reproducible via stored `resampling` and `seed`.
- Predictions always include model identifiers; ensemble predictions identified with `model = "ensemble"` and `ensemble = TRUE`.

## Predictions & Evaluations
- `predict()` returns tidy tibbles; supports `include` to filter models and ensemble.
- `ml4t2e_evaluate()` accepts either a `t2e_fit` or a predictions tibble; returns a tidy metrics table with columns: `model`, `metric`, `value`, `split` (train/test/cv), and optional `cause`.
- For CR tasks, specify `event_of_interest` in `controls` or as an argument to evaluation/plotting when needed.

## Summaries
- `print(fit)`: compact overview: models, time grid summary, resampling, ensemble method, basic metrics.
- `summary(fit)`: detailed tibble: per-model metrics, weights (if ensemble), runtime, package versions, and parameter snapshots (top-level only).

## XAI Components
- `ml4t2e_explain()` provides a unified interface:
  - `method = "auto"` chooses SHAP if supported, else permutation importance.
  - Returns a tidy explanation tibble tagged with class `t2e_explain` and columns like `id`, `feature`, `value`, `importance`, optional `time`.
  - `include` filters which models to explain (e.g., just `ensemble` or a specific engine).
- `autoplot(t2e_explain, what = "importance"|"dependence"|"waterfall")` plots common explanation views.

## Persistence
- `ml4t2e_save(fit, path)`: saves via `saveRDS` with fields `api_version`, `time_grid`, `time_units`, `registry_versions`, `controls`, and engine-specific states.
- `ml4t2e_load(path)`: loads and validates `api_version`, warns on newer/older engines with actionable guidance.
- Optional `ml4t2e_save_preds(pred, path)` for caching large prediction tibbles (kept internal unless needed).

## Object Interface Details (Concise)
- `TimeToEventModel` (R6): `initialize(spec)`, `fit(task)`, `predict_survival(newdata, times)`, `predict_risk(newdata, times)`, `predict_time(q)`, `model_info()`.
- `SurvivalModel` and `CompetingRiskModel` extend base and can share mixins for preprocessing and validation.
- `Task` objects (`t2e_task_surv`, `t2e_task_cr`) store: roles (id, time, event/status/cause), features, `time_units`, `cause_map`, and `data_hash`.
- `t2e_fit`: holds `task`, list of fitted engines, optional `t2e_ensemble`, `time_grid`, `seed`, `api_version`.
- S3: `predict.t2e_fit()` returns a `t2e_pred` tibble; `autoplot.t2e_pred()` builds ggplot2 objects; `summary.t2e_fit()`/`print.t2e_fit()` show concise model/metric info.

## Naming & Style Conventions
- Files: lower_snake_case; classes: UpperCamelCase; functions: lower_snake_case.
- Internal-only helpers: leading dot (e.g., `.check_task()`), `@keywords internal`, no `@export`.
- Package hooks in `R/core/zzz.R` (e.g., register engines on load).
