# Test script for ml4time2event Vignette Code

# Load necessary libraries
# Ensure these are installed: install.packages(c("survival", "dplyr", "recipes", "rsample", "ggplot2", "riskRegression"))
library(ml4time2event)
library(survival)
library(dplyr)
library(recipes)
library(rsample)
library(ggplot2)
library(riskRegression) # Needed for pbc data and some evaluation metrics

# --- Survival Analysis Example using 'colon' dataset ---
message("--- Starting Survival Analysis Example (colon dataset) ---")

# 1. Load and Preprocess Data
library(survival) # Load library *before* checking/using
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' needed for this example to run. Please install it.", call. = FALSE)
}
utils::data(colon, package = "survival") # Use utils::data
colon_df <- colon %>%
  filter(etype == 2) %>% # Focus on death event for standard survival
  mutate(
    status = case_when(
      status == 1 ~ 1, # Death
      TRUE ~ 0        # Censored
    ),
    sex = factor(sex, levels = c(0, 1), labels = c("female", "male")),
    obstruct = factor(obstruct),
    perfor = factor(perfor),
    adhere = factor(adhere),
    nodes = as.numeric(nodes), # Ensure nodes is numeric
    differ = factor(differ),
    extent = factor(extent),
    surg = factor(surg),
    node4 = factor(node4)
  ) %>%
  select(-etype, -id) %>% # Remove unnecessary columns
  na.omit() # Keep complete cases for simplicity in this example

# Ensure time is positive
colon_df <- colon_df %>% filter(time > 0)

# Define outcome and features
outcome_vars_surv <- c("time", "status")
feature_vars_surv <- setdiff(names(colon_df), outcome_vars_surv)

# 2. Split Data
set.seed(123)
# Use strata argument with the status variable, remove outcome_vars from call
# Assume t2edata_split returns a list: [[1]] train, [[2]] test
split_list_surv <- t2edata_split(colon_df, prop = 0.7, strata = "status")
train_data_surv <- split_list_surv[[1]] # Access train data from list
test_data_surv <- split_list_surv[[2]]  # Access test data from list

# 3. Create Data Recipe
# Using a minimal recipe for demonstration
# Add expvar argument
surv_recipe <- t2emodel_data_recipe_init(traindata = train_data_surv, timevar = "time", eventvar = "status", expvar = feature_vars_surv, idvars = NULL) %>%
  step_impute_median(all_numeric_predictors()) %>% # Impute numeric NAs
  step_impute_mode(all_nominal_predictors()) %>% # Impute categorical NAs
  step_dummy(all_nominal_predictors()) %>% # Create dummy variables
  step_nzv(all_predictors()) # Remove near-zero variance predictors

# Prep and bake
prepped_recipe_surv <- prep(surv_recipe, training = train_data_surv)
baked_train_surv <- bake(prepped_recipe_surv, new_data = train_data_surv)
baked_test_surv <- bake(prepped_recipe_surv, new_data = test_data_surv)

# 4. Train Survival Models
# Define models to run
surv_model_list <- c(
  "Cox", "RF", "glmnet", "xgboost", "BART", "GAM", "rulefit", "gbm", "SurvReg"
)

# Define evaluation times (e.g., quantiles of event times in test set)
eval_times_surv <- quantile(baked_test_surv$time[baked_test_surv$status == 1], probs = c(0.25, 0.5, 0.75))

# Run models
# Note: BART can be slow, consider removing if time is an issue
# Note: GAM might require specific formula handling depending on features
# Note: Rulefit can also be computationally intensive
# For demonstration, let's run a subset
surv_model_subset <- c("Cox", "RF", "glmnet", "xgboost") # "BART", "GAM", "rulefit", "gbm", "SurvReg"

# Use datatrain, timevar, eventvar arguments
trained_surv_models <- RunSurvModels(
  datatrain = baked_train_surv,
  timevar = "time",
  eventvar = "status",
  ExpVars = feature_vars_surv, # Add ExpVars argument
  models_to_run = surv_model_subset,
  eval_times = eval_times_surv,
  n_cores = 1 # Set to > 1 for parallel processing if desired
  # Add other relevant parameters like hyperparameters if needed
)

message("Survival models trained.")

# 5. Predict Survival Probabilities
# Use timevar, eventvar arguments
surv_predictions <- PredictSurvModels(
  trained_models = trained_surv_models,
  new_data = baked_test_surv,
  timevar = "time",
  eventvar = "status",
  predict_horizon = eval_times_surv
)

message("Survival predictions generated.")
# print(head(surv_predictions)) # Inspect predictions

# 6. Evaluate Survival Models
# Time-dependent Concordance Index
# Use timevar, eventvar arguments
surv_cindex <- timedepConcordance(
  predictions = surv_predictions,
  test_data = baked_test_surv,
  timevar = "time",
  eventvar = "status",
  predict_horizon = eval_times_surv
)
message("Survival C-index calculated.")
print("Survival C-index:")
print(surv_cindex)

# Integrated Brier Score (requires riskRegression package or similar logic)
# Need predicted survival probabilities at various times
# Let's predict at more time points for IBS
fine_eval_times_surv <- seq(min(eval_times_surv), max(eval_times_surv), length.out = 50)
# Use timevar, eventvar arguments
surv_predictions_fine <- PredictSurvModels(
  trained_models = trained_surv_models,
  new_data = baked_test_surv,
  timevar = "time",
  eventvar = "status",
  predict_horizon = fine_eval_times_surv
)

# Calculate IBS (using ml4time2event's integratedBrier if available and adapted, or riskRegression)
# Assuming ml4time2event::integratedBrier exists and works similarly
# ibs_scores_surv <- integratedBrier(
#   predictions = surv_predictions_fine, # Needs matrix/df format: rows=subjects, cols=times
#   test_data = baked_test_surv,
#   outcome_vars = outcome_vars_surv,
#   predict_horizon = fine_eval_times_surv
# )
# message("Survival IBS calculated (placeholder).")
# print(ibs_scores_surv)

# Using riskRegression::Score for IBS as an alternative/verification
# Prepare data for Score function
surv_formula <- Surv(time, status) ~ 1
pred_list_rr_surv <- lapply(surv_model_subset, function(model_name) {
    pred_df <- surv_predictions_fine %>% filter(model == model_name)
    # Reshape: rows = subjects, cols = times
    pred_matrix <- tidyr::pivot_wider(pred_df, names_from = predict_horizon, values_from = surv_prob) %>%
                   select(-model, -any_of(c("time", "status", names(baked_test_surv)))) %>% # Keep only prediction columns
                   as.matrix()
    # Ensure order matches test data
    # This assumes predictions are returned in the order of new_data
    return(1 - pred_matrix) # Score expects risk (1 - survival prob)
})
names(pred_list_rr_surv) <- surv_model_subset

# Calculate IBS using riskRegression
ibs_calc_surv <- tryCatch({
  riskRegression::Score(
    object = pred_list_rr_surv,
    formula = surv_formula,
    data = baked_test_surv,
    metrics = "ibs",
    times = fine_eval_times_surv,
    summary = "ibs",
    conf.int = FALSE # Faster calculation
  )$Brier$score
}, error = function(e) {
  message("Error calculating IBS with riskRegression: ", e$message)
  return(NULL)
})

message("Survival IBS calculated using riskRegression.")
print("Survival IBS (lower is better):")
print(ibs_calc_surv)


# 7. Ensemble Prediction (Example: Averaging)
# Predict survival probabilities from all models at specific horizons
# Use the predictions already generated (surv_predictions)
ensemble_surv_preds <- PredictAllPossibleOutcomesSurvOrCifs(
  predictions = surv_predictions,
  ensemble_method = "average", # or "weighted_average" with weights
  outcome_type = "survival"
)

message("Survival ensemble predictions generated.")
# print(head(ensemble_surv_preds))

# Evaluate ensemble (e.g., C-index) - requires formatting like other predictions
# Add 'model' column and combine
ensemble_surv_preds$model <- "Ensemble_Avg"
all_surv_preds_for_eval <- bind_rows(surv_predictions, ensemble_surv_preds)

# Use timevar, eventvar arguments
ensemble_cindex_surv <- timedepConcordance(
  predictions = all_surv_preds_for_eval,
  test_data = baked_test_surv,
  timevar = "time",
  eventvar = "status",
  predict_horizon = eval_times_surv
)
message("Survival Ensemble C-index calculated.")
print("Survival C-index (including Ensemble):")
print(ensemble_cindex_surv)


# 8. ML Pipeline Example (Recipe + Model Training + Prediction)
# This essentially combines steps 3, 4, 5
# Update function signature and internal calls to use timevar/eventvar
run_survival_pipeline <- function(train_df, test_df, timevar, eventvar, feature_vars, model_type, eval_times) {
  message(paste("Running pipeline for:", model_type))
  # 1. Recipe
  # Note: Using standard recipe here, not t2emodel_data_recipe_init
  recipe_obj <- recipe(train_df) %>%
    update_role(all_of(timevar), all_of(eventvar), new_role = "outcome") %>%
    update_role(all_of(feature_vars), new_role = "predictor") %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>% # Use one_hot for compatibility
    step_nzv(all_predictors())

  prepped_recipe <- prep(recipe_obj, training = train_df)
  baked_train <- bake(prepped_recipe, new_data = train_df)
  baked_test <- bake(prepped_recipe, new_data = test_df)

  # 2. Train Model
  # Use datatrain, timevar, eventvar, ExpVars arguments
  trained_model_list <- RunSurvModels(
    datatrain = baked_train,
    timevar = timevar,
    eventvar = eventvar,
    ExpVars = setdiff(names(baked_train), c(timevar, eventvar)), # Derive ExpVars within pipeline
    models_to_run = c(model_type), # Run only the specified model
    eval_times = eval_times,
    n_cores = 1
  )

  # 3. Predict
  # Use timevar, eventvar arguments
  predictions <- PredictSurvModels(
    trained_models = trained_model_list,
    new_data = baked_test,
    timevar = timevar,
    eventvar = eventvar,
    predict_horizon = eval_times
  )

  # 4. Evaluate (optional, example C-index)
  # Use timevar, eventvar arguments
  c_index <- timedepConcordance(
    predictions = predictions,
    test_data = baked_test,
    timevar = timevar,
    eventvar = eventvar,
    predict_horizon = eval_times
  )
  message(paste("Pipeline C-index for", model_type, ":"))
  print(c_index)

  return(list(model = trained_model_list, predictions = predictions, c_index = c_index))
}

# Run pipeline for one model, passing timevar/eventvar
pipeline_rf_surv <- run_survival_pipeline(
  train_df = train_data_surv,
  test_df = test_data_surv,
  timevar = "time",
  eventvar = "status",
  feature_vars = feature_vars_surv,
  model_type = "RF",
  eval_times = eval_times_surv
)
message("Survival ML Pipeline example completed for RF.")


# --- Competing Risks Analysis Example using 'pbc' dataset ---
message("\n--- Starting Competing Risks Example (pbc dataset) ---")

# 1. Load and Preprocess Data
library(survival) # Load library *before* checking/using
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' needed for this example to run. Please install it.", call. = FALSE)
}
utils::data(pbc, package = "survival") # Use utils::data
# For competing risks, need time, event status (0=censored, 1=event1, 2=event2, ...), features
# Events: 1=death, 2=transplant
pbc_df <- pbc %>%
  mutate(
    cr_status = case_when(
      status == 2 ~ 1, # Death is event 1
      status == 1 ~ 2, # Transplant is event 2
      status == 0 ~ 0  # Censored
    ),
    # Convert relevant columns to factors/numeric
    sex = factor(sex),
    ascites = factor(ascites),
    hepato = factor(hepato),
    spiders = factor(spiders),
    edema = factor(edema),
    stage = factor(stage)
  ) %>%
  select(-id, -status) %>% # Remove original status and id
  na.omit() # Keep complete cases for simplicity

# Ensure time is positive
pbc_df <- pbc_df %>% filter(time > 0)

# Define outcome and features
outcome_vars_cr <- c("time", "cr_status")
feature_vars_cr <- setdiff(names(pbc_df), outcome_vars_cr)

# 2. Split Data
set.seed(456)
# Use strata argument with the status variable, remove outcome_vars from call
# Assume t2edata_split returns a list: [[1]] train, [[2]] test
split_list_cr <- t2edata_split(pbc_df, prop = 0.7, strata = "cr_status")
train_data_cr <- split_list_cr[[1]] # Access train data from list
test_data_cr <- split_list_cr[[2]]  # Access test data from list

# 3. Create Data Recipe
# Add expvar argument
cr_recipe <- t2emodel_data_recipe_init(traindata = train_data_cr, timevar = "time", eventvar = "cr_status", expvar = feature_vars_cr, idvars = NULL) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors())

# Prep and bake
prepped_recipe_cr <- prep(cr_recipe, training = train_data_cr)
baked_train_cr <- bake(prepped_recipe_cr, new_data = train_data_cr)
baked_test_cr <- bake(prepped_recipe_cr, new_data = test_data_cr)

# 4. Train Competing Risks Models
# Define models to run
cr_model_list <- c("FineGray", "Cox", "RF", "BART", "rulefit") # Add others as available

# Define evaluation times
eval_times_cr <- quantile(baked_test_cr$time[baked_test_cr$cr_status != 0], probs = c(0.25, 0.5, 0.75))

# Run models (subset for speed)
cr_model_subset <- c("FineGray", "Cox", "RF") # "BART", "rulefit"

# Use datatrain, timevar, eventvar arguments
trained_cr_models <- RunCRModels(
  datatrain = baked_train_cr,
  timevar = "time",
  eventvar = "cr_status",
  ExpVars = feature_vars_cr, # Add ExpVars argument
  models_to_run = cr_model_subset,
  eval_times = eval_times_cr,
  event_of_interest = 1, # Example: Focus on event 1 (death) for evaluation like C-index
  n_cores = 1
)

message("Competing risks models trained.")

# 5. Predict Cumulative Incidence Functions (CIFs)
# Use timevar, eventvar arguments
cr_predictions <- PredictCRModels(
  trained_models = trained_cr_models,
  new_data = baked_test_cr,
  timevar = "time",
  eventvar = "cr_status",
  predict_horizon = eval_times_cr,
  event_of_interest = "all" # Predict CIF for all events
)

message("Competing risks predictions generated.")
# print(head(cr_predictions)) # Inspect predictions (will have cif_event_1, cif_event_2, etc.)

# 6. Evaluate Competing Risks Models
# Time-dependent Concordance Index for a specific event
# Use timevar, eventvar arguments
cr_cindex_event1 <- timedepConcordanceCR(
  predictions = cr_predictions,
  test_data = baked_test_cr,
  timevar = "time",
  eventvar = "cr_status",
  predict_horizon = eval_times_cr,
  event_of_interest = 1 # Evaluate C-index for event 1 (death)
)
message("Competing Risks C-index (Event 1) calculated.")
print("Competing Risks C-index (Event 1):")
print(cr_cindex_event1)

# Integrated Brier Score for CR (requires careful handling of multiple events)
# Using riskRegression::Score as an example
# Prepare data for Score function
cr_formula <- Hist(time, cr_status) ~ 1 # Using Hist for competing risks
pred_list_rr_cr <- lapply(cr_model_subset, function(model_name) {
    pred_df <- cr_predictions %>% filter(model == model_name, event_type == 1) # Focus on CIF for event 1
    # Reshape: rows = subjects, cols = times
    pred_matrix <- tidyr::pivot_wider(pred_df, names_from = predict_horizon, values_from = cif_pred) %>%
                   select(-model, -event_type, -any_of(c("time", "cr_status", names(baked_test_cr)))) %>%
                   as.matrix()
    return(pred_matrix) # Score expects risk (CIF)
})
names(pred_list_rr_cr) <- cr_model_subset

# Calculate IBS using riskRegression for event 1
ibs_calc_cr <- tryCatch({
  riskRegression::Score(
    object = pred_list_rr_cr,
    formula = cr_formula,
    data = baked_test_cr,
    metrics = "ibs",
    times = eval_times_cr, # Use eval times for CR IBS
    summary = "ibs",
    cause = 1, # Specify the cause for IBS
    conf.int = FALSE
  )$Brier$score
}, error = function(e) {
  message("Error calculating CR IBS with riskRegression: ", e$message)
  return(NULL)
})

message("Competing Risks IBS (Event 1) calculated using riskRegression.")
print("Competing Risks IBS for Event 1 (lower is better):")
print(ibs_calc_cr)


# 7. Ensemble Prediction (Example: Averaging CIFs)
ensemble_cr_preds <- PredictAllPossibleOutcomesSurvOrCifs(
  predictions = cr_predictions,
  ensemble_method = "average",
  outcome_type = "competing_risks"
)

message("Competing risks ensemble predictions generated.")
# print(head(ensemble_cr_preds))

# Evaluate ensemble (e.g., C-index for event 1)
ensemble_cr_preds$model <- "Ensemble_Avg"
all_cr_preds_for_eval <- bind_rows(cr_predictions, ensemble_cr_preds)

# Use timevar, eventvar arguments
ensemble_cindex_cr1 <- timedepConcordanceCR(
  predictions = all_cr_preds_for_eval,
  test_data = baked_test_cr,
  timevar = "time",
  eventvar = "cr_status",
  predict_horizon = eval_times_cr,
  event_of_interest = 1
)
message("Competing Risks Ensemble C-index (Event 1) calculated.")
print("Competing Risks C-index for Event 1 (including Ensemble):")
print(ensemble_cindex_cr1)


# 8. ML Pipeline Example (Recipe + CR Model Training + Prediction)
# Update function signature and internal calls to use timevar/eventvar
run_cr_pipeline <- function(train_df, test_df, timevar, eventvar, feature_vars, model_type, eval_times, event_interest) {
  message(paste("Running CR pipeline for:", model_type))
  # 1. Recipe
  # Note: Using standard recipe here, not t2emodel_data_recipe_init
  recipe_obj <- recipe(train_df) %>%
    update_role(all_of(timevar), all_of(eventvar), new_role = "outcome") %>%
    update_role(all_of(feature_vars), new_role = "predictor") %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_nzv(all_predictors())

  prepped_recipe <- prep(recipe_obj, training = train_df)
  baked_train <- bake(prepped_recipe, new_data = train_df)
  baked_test <- bake(prepped_recipe, new_data = test_df)

  # 2. Train Model
  # Use datatrain, timevar, eventvar, ExpVars arguments
  trained_model_list <- RunCRModels(
    datatrain = baked_train,
    timevar = timevar,
    eventvar = eventvar,
    ExpVars = setdiff(names(baked_train), c(timevar, eventvar)), # Derive ExpVars within pipeline
    models_to_run = c(model_type),
    eval_times = eval_times,
    event_of_interest = event_interest, # Pass event of interest
    n_cores = 1
  )

  # 3. Predict
  # Use timevar, eventvar arguments
  predictions <- PredictCRModels(
    trained_models = trained_model_list,
    new_data = baked_test,
    timevar = timevar,
    eventvar = eventvar,
    predict_horizon = eval_times,
    event_of_interest = "all" # Predict all CIFs
  )

  # 4. Evaluate (optional, example C-index for event_interest)
  # Use timevar, eventvar arguments
  c_index <- timedepConcordanceCR(
    predictions = predictions,
    test_data = baked_test,
    timevar = timevar,
    eventvar = eventvar,
    predict_horizon = eval_times,
    event_of_interest = event_interest
  )
  message(paste("Pipeline CR C-index for", model_type, "(Event", event_interest, "):"))
  print(c_index)

  return(list(model = trained_model_list, predictions = predictions, c_index = c_index))
}

# Run pipeline for one CR model, passing timevar/eventvar
pipeline_fg_cr <- run_cr_pipeline(
  train_df = train_data_cr,
  test_df = test_data_cr,
  timevar = "time",
  eventvar = "cr_status",
  feature_vars = feature_vars_cr,
  model_type = "FineGray",
  eval_times = eval_times_cr,
  event_interest = 1 # Focus on event 1
)
message("Competing Risks ML Pipeline example completed for FineGray.")

message("\n--- Vignette Test Script Finished ---")
