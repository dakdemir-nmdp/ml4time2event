# Copilot Instructions for ml4time2event

## Overview
This is an R package for machine learning in time-to-event analysis, providing survival and competing risks models. Key components include model training (`SurvModel_*`, `CRModel_*`), prediction (`Predict_*`), data preprocessing (imputation, recipes), and evaluation (concordance, Brier score).

## Architecture
- **Models**: Separate files for each algorithm (e.g., `R/surv_cox.R` for Cox models, `R/cr_bart.R` for competing risks BART). Models return consistent lists: `list(model, times, varprof)`.
- **Data Flow**: Data loading (`readData`) → Cleaning/Imputation (`ImputeMissinRecordsData`, recipes) → Splitting (`t2edata_split`) → Training → Prediction → Evaluation (`timedepConcordance`, `integratedBrier`).
- **Utilities**: Variable profiling (`VariableProfile`), metrics (`BrierScore`).
- **Structure**: Functional organization in `R/`; examples in `examples/`; tests in `tests/testthat/`.

## Critical Development Workflows
- **R Path**: Use `/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/R` for terminal commands
- **Development**: `devtools::load_all()` to load package, `devtools::test()` for tests, `lintr::lint_package()` for style
- **Testing**: `devtools::test(filter = 'test-name')` for specific tests; tests organized by functionality: `test-surv-models/` (survival models), `test-cr-models/` (competing risks), `test-data/` (data processing)
- **Model Training**: Call `SurvModel_Cox(data, expvars, timevar, eventvar, ...)` with `...` for flexible parameter passing
- **Prediction**: Use `Predict_SurvModel_Cox(modelout, newdata)`; returns `list(Probs, Times)`
- **Evaluation**: Pass predictions to `timedepConcordance(predsurv, predsurvtimes, obstimes, obsevents)`

## Package-Specific Patterns
- **Function Naming**: All model functions follow `SurvModel_*` or `CRModel_*` pattern, predictions use `Predict_*` prefix
- **Return Structures**: Models always return `list(model, times, varprof)`, predictions return `list(Probs, Times)`
- **Variable Handling**: Use `VariableProfile()` for consistent factor/numeric profiling; store factor levels for prediction consistency
- **Data I/O**: `readData()` supports csv/xlsx/sas7bdat with automatic lowercase column conversion
- **Survival Objects**: Always use `survival::Surv(timevar, eventvar)` for time-to-event outcomes
- **Error Handling**: Validate inputs early with `stop()` messages; use `tryCatch()` for optional package dependencies

## Integration Points & Dependencies
- **Core Dependencies**: survival, glmnet, xgboost, recipes, partykit - check DESCRIPTION for full list
- **Optional Models**: BART, randomForestSRC, mgcv - wrapped in `requireNamespace()` checks
- **Data Processing**: Heavy use of dplyr pipes, recipes for preprocessing, missRanger for imputation
- **Testing Framework**: testthat with nested directory structure (`test-surv-models/`, `test-cr-models/`)

## Code Style (see CLAUDE.md)
- **Naming**: CamelCase for exported functions (`SurvModel_Cox`), snake_case for files/internal vars
- **Documentation**: Roxygen2 with `@title`, `@param`, `@return`; use `@importFrom` for all external functions
- **Formatting**: 2-space indent, `<-` for assignment, max 80 chars per line
- **Function Structure**: Validate inputs first, return named lists consistently

## Examples
- **Basic Workflow**: `examples/01_basic_survival.R` shows complete Cox model pipeline with PBC data
- **Competing Risks**: Use `CRModel_FineGray` for subdistribution hazards, `PredictCRModels` for predictions
- **Data Recipes**: `minimal_data_recipe`, `prep_data_recipe` for feature engineering and imputation

Reference `CLAUDE.md` for detailed style guidelines. Focus on modular, testable functions with consistent API patterns.