# ml4time2event Package Functions Reference

## Data Loading and Preparation

### Data Loading
- **`get_lung_survival_data()`** - Load prepared NCCTG lung cancer survival dataset
- **`get_bmt_competing_risks_data()`** - Load prepared BMT competing risks dataset
- **`readData(file)`** - Read data from CSV or other formats

### Data Cleaning
- **`ImputeMissinRecordsData(data, dontuse, num.trees)`** - Impute missing values using random forest

### Data Splitting
- **`t2edata_split(data, prop)`** - Split data into train/test sets for time-to-event analysis

### Data Transformation
- **`ZeroOneScalerData(data)`** - Scale numeric variables to [0,1] range
- **`ZeroOneScalerApplierData(data, mins, maxs)`** - Apply pre-calculated scaling
- **`UndoZeroOneScalerApplierData(data, mins, maxs)`** - Reverse scaling transformation
- **`NumVarstCatsData(data, numgroups, cuts)`** - Convert numeric to categorical variables

### Data Summary and Profiling
- **`VariableProfile(data, variables)`** - Generate variable profiles and summaries
- **`SummaryTableDataSummmary(data, UseVars)`** - Create gtsummary tables
- **`pairwiserelationshipsDataSummmary(data)`** - Calculate pairwise distance correlations
- **`gethighcorvarsDataSummmary(pmat, corcutoff)`** - Identify highly correlated variable pairs
- **`OneAgainstRestCorDataSummmary(data)`** - Calculate each variable's correlation against all others

## Model Training

### Survival Models
- **`RunSurvModels(datatrain, ExpVars, timevar, eventvar, models, ...)`** - Train ensemble of survival models
  - Available models: "coxph", "glmnet", "rulefit", "xgboost", "gam", "gbm", "ExpSurvReg", "bart", "shallownn", "ttah", "RF"
  - Returns: SurvEnsemble object containing all fitted models

### Competing Risks Models
- **`RunCRModels(datatrain, ExpVars, timevar, eventvar, models, ...)`** - Train ensemble of competing risks models
  - Available models: "cox", "FG", "xgboost", "gam", "bart", "rulefit", "survreg", "shallownn", "ttah", "RF"
  - Returns: CREnsemble object containing all fitted models

### Pipeline API (High-Level)
- **`ml4t2e_fit_pipeline(data, analysis_type, timevar, eventvar, models, ...)`** - Fit complete preprocessing + modeling pipeline
  - analysis_type: "survival" or "competing_risks"
  - Returns: ml4t2e_pipeline object

## Prediction

### Survival Predictions
- **`PredictSurvModels(models, newdata, new_times, ensemble_method)`** - Generate survival predictions
  - Returns: List with ModelPredictions, NewProbs (ensemble), NewTimes, models_used

### Competing Risks Predictions
- **`PredictCRModels(models, newdata, new_times, ensemble_method)`** - Generate CIF predictions
  - Returns: List with ModelPredictions, NewProbs (ensemble), NewTimes, models_used

### Outcome Predictions
- **`PredictAllPossibleOutcomesSurvOrCifs(data, modelslist, modeltypes, times)`** - Predict from mixed model types
- **`calculate_expected_time_lost(times, event_probs, upper_limit, lower_limit)`** - Calculate expected time lost (ETL) by integrating event probabilities

## Model Evaluation

### Unified Evaluation Functions
- **`EvaluateSurvModels(models, data, timevar, eventvar, eval_times, prediction_times)`** - Comprehensive survival model evaluation
  - Returns: Leaderboard with integrated_brier, integrated_c, brier_scores, c_index, rank
- **`EvaluateCRModels(models, data, timevar, eventvar, eval_times, prediction_times, cause)`** - Comprehensive competing risks evaluation

### Survival Metrics
- **`BrierScore(predsurv, pred_times, obstimes, obsevents, eval_times)`** - Brier score at specific times
- **`integratedBrier(predsurv, pred_times, obstimes, obsevents, eval_times)`** - Integrated Brier score
- **`timedepConcordance(predsurv, pred_times, obstimes, obsevents)`** - Time-dependent concordance
- **`integratedC(predsurv, pred_times, obstimes, obsevents)`** - Integrated concordance index

### Competing Risks Metrics
- **`BrierScoreCR(SurvObj, Predictions, eval_times, cause, pred_times)`** - CIF Brier score
- **`integratedBrierCR(SurvObj, Predictions, eval_times, cause, pred_times)`** - Integrated CIF Brier score
- **`timedepConcordanceCR(SurvObj, Predictions, time, cause, pred_times)`** - Time-dependent concordance for CR
- **`integratedConcordanceCR(SurvObj, Predictions, eval_times, cause, pred_times)`** - Integrated concordance for CR

## Visualization

### Survival/CIF Curve Plotting
- **`plot_survival_curves(predictions, model_names, patients_to_plot, colors, highlight_ensemble, ...)`** - Plot survival curves
  - Supports single or multiple models
  - Automatic ensemble highlighting
  - Faceting by patient
- **`plot_cif_curves(predictions, model_names, patients_to_plot, colors, highlight_ensemble, ...)`** - Plot CIF curves
  - Similar interface to plot_survival_curves

### Performance Visualization
- **`plot_model_performance(performance_df, metric, highlight_ensemble, flip_coords)`** - Create performance comparison plots
  - Metrics: "concordance", "brier", "prediction_error", "expected_time_lost"
  - Supports horizontal bar charts

## Model Explainability (SHAP)

### SHAP Calculation
- **`ml4t2e_calculate_shap(pipeline, data, time_horizon, background, nsim)`** - Calculate SHAP values
  - Returns: ml4t2e_shap object with shap_values, baseline, predictions, feature_values
- **`ml4t2e_shap_predict_fn(pipeline, time_horizon, ensemble_method)`** - Create prediction function for SHAP

### SHAP Visualization
- **`ml4t2e_shap_importance(shap_result, max_features, plot_type)`** - Variable importance plot
  - plot_type: "beeswarm" or "bar"
- **`ml4t2e_shap_dependence(shap_result, feature, color_by)`** - Feature dependence plot
- **`ml4t2e_shap_waterfall(shap_result, obs_id, max_features)`** - Individual prediction explanation

## Model Persistence

### Save/Load Functions
- **`SaveEnsemble(ensemble_model, file, compress)`** - Save ensemble to disk
- **`LoadEnsemble(file)`** - Load ensemble from disk
- **`SurvEnsemble(ensemble_list)`** - Create SurvEnsemble S3 class
- **`CREnsemble(ensemble_list)`** - Create CREnsemble S3 class
- **`is.SurvEnsemble(x)`** - Check if object is SurvEnsemble
- **`is.CREnsemble(x)`** - Check if object is CREnsemble

## Helper Utilities

### Mathematical Utilities
- **`Integrator(times, scores, minmax, scale)`** - Numerical integration helper

### Model Utilities
- **`prepare_newdata_for_model(newdata, train_data, ExpVars)`** - Prepare prediction data
- **`baseline_cumhaz_stepfun(fit)`** - Extract baseline cumulative hazard
- **`cifMatInterpolaltor(probsMat, times, new_times)`** - Interpolate CIF matrices

## Individual Model Functions

### Survival Models (Individual)
- **`SurvModel_Cox(data, ExpVars, timevar, eventvar, ...)`** - Cox proportional hazards
- **`SurvModel_glmnet(data, ExpVars, timevar, eventvar, ...)`** - Elastic net regularized Cox
- **`SurvModel_RF(data, ExpVars, timevar, eventvar, ...)`** - Random forest
- **`SurvModel_xgboost(data, ExpVars, timevar, eventvar, ...)`** - XGBoost
- **`SurvModel_GAM(data, ExpVars, timevar, eventvar, ...)`** - Generalized additive model
- **`SurvModel_gbm(data, ExpVars, timevar, eventvar, ...)`** - Gradient boosting
- **`SurvModel_BART(data, ExpVars, timevar, eventvar, ...)`** - Bayesian additive regression trees
- **`SurvModel_ShallowNN(data, ExpVars, timevar, eventvar, ...)`** - Shallow neural network
- **`SurvModel_TTAH(data, ExpVars, timevar, eventvar, ...)`** - Time-to-any-hazard model
- **`SurvModel_rulefit(data, ExpVars, timevar, eventvar, ...)`** - RuleFit
- **`SurvModel_SurvReg(data, ExpVars, timevar, eventvar, dist)`** - Parametric survival regression

### Competing Risks Models (Individual)
- **`CRModel_Cox(data, ExpVars, timevar, eventvar, ...)`** - Cause-specific Cox
- **`CRModel_FineGray(data, ExpVars, timevar, eventvar, ...)`** - Fine-Gray subdistribution hazard
- **`CRModel_RF(data, ExpVars, timevar, eventvar, ...)`** - Random forest for CR
- **`CRModel_xgboost(data, ExpVars, timevar, eventvar, ...)`** - XGBoost for CR
- **`CRModel_GAM(data, ExpVars, timevar, eventvar, ...)`** - GAM for CR
- **`CRModel_BART(data, ExpVars, timevar, eventvar, ...)`** - BART for CR
- **`CRModel_ShallowNN(data, ExpVars, timevar, eventvar, ...)`** - Shallow NN for CR
- **`CRModel_TTAH(data, ExpVars, timevar, eventvar, ...)`** - TTAH for CR
- **`CRModel_rulefit(data, ExpVars, timevar, eventvar, ...)`** - RuleFit for CR
- **`CRModel_SurvReg(data, ExpVars, timevar, eventvar, ...)`** - Parametric CR model

### Prediction Functions (Individual Models)
Each model has corresponding prediction functions:
- **`Predict_SurvModel_*`** - Predict from individual survival models
- **`Predict_CRModel_*`** - Predict from individual competing risks models

## Key Function Categories Summary

1. **Data Pipeline**: Load → Clean → Transform → Profile → Split
2. **Modeling Pipeline**: Train → Predict → Evaluate → Visualize → Save
3. **Explainability Pipeline**: Fit → Calculate SHAP → Visualize Importance/Dependence/Waterfall
4. **Evaluation Pipeline**: Predict → Metrics (Brier, Concordance, ETL) → Leaderboard → Plot

## Typical Workflow

### Survival Analysis
```r
# 1. Load and prepare data
data <- get_lung_survival_data()
split <- t2edata_split(data, prop = 0.7)

# 2. Train ensemble
models <- RunSurvModels(
  datatrain = split$Train,
  ExpVars = c("age", "sex", "ph.karno"),
  timevar = "time",
  eventvar = "status",
  models = c("coxph", "glmnet", "xgboost")
)

# 3. Generate predictions
preds <- PredictSurvModels(
  models = models,
  newdata = split$Test,
  new_times = NULL,  # Use model's default grid
  ensemble_method = "average"
)

# 4. Evaluate
leaderboard <- EvaluateSurvModels(
  models = models,
  data = split$Test,
  timevar = "time",
  eventvar = "status"
)

# 5. Visualize
plot_survival_curves(
  predictions = list(
    Cox = list(Probs = preds$ModelPredictions$CPH_Model, Times = preds$NewTimes),
    Ensemble = list(Probs = preds$NewProbs, Times = preds$NewTimes)
  ),
  patients_to_plot = 1:3
)

# 6. Save
SaveEnsemble(models, file = "my_ensemble.rds")
```

### Pipeline API (Simplified)
```r
# All-in-one pipeline
pipeline <- ml4t2e_fit_pipeline(
  data = data,
  analysis_type = "survival",
  timevar = "time",
  eventvar = "status",
  models = c("coxph", "glmnet")
)

# SHAP explainability
shap_result <- ml4t2e_calculate_shap(
  pipeline = pipeline,
  data = data[1:50, ],
  time_horizon = 365
)

ml4t2e_shap_importance(shap_result)
```
