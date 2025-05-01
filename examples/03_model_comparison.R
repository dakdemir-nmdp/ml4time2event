# Example 03: Model Comparison using ml4time2event

# --- 1. Load Libraries ---
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
# Install dependencies for models you want to compare (e.g., glmnet, xgboost, etc.)
# if (!requireNamespace("glmnet", quietly = TRUE)) { install.packages("glmnet") }
# if (!requireNamespace("xgboost", quietly = TRUE)) { install.packages("xgboost") }
# if (!requireNamespace("randomForestSRC", quietly = TRUE)) { install.packages("randomForestSRC") }
# ... add others as needed

library(survival)
library(dplyr)
library(pec) # For plotting C-index

# --- Load the ml4time2event package ---
# IMPORTANT: Assumes package functions are available.
# Use devtools::load_all() or library(ml4time2event).

# --- 2. Load and Prepare Data ---
data(pbc, package = "survival")
pbc_baseline <- pbc[!duplicated(pbc$id), ] %>%
  filter(time > 0) %>%
  mutate(
    status_event = ifelse(status == 2, 1, 0),
    sex = factor(sex), stage = factor(stage)
  ) %>%
  select(-id, -status)

time_var <- "time"
event_var <- "status_event"
predictor_vars <- setdiff(names(pbc_baseline), c(time_var, event_var))

# Handle missing values (simple median/mode imputation)
for (col in predictor_vars) {
  if (is.numeric(pbc_baseline[[col]])) {
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- median(pbc_baseline[[col]], na.rm = TRUE)
  } else if (is.factor(pbc_baseline[[col]])) {
    mode_val <- names(which.max(table(pbc_baseline[[col]])))
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
  }
   if (is.character(pbc_baseline[[col]])) {
     mode_val <- names(which.max(table(pbc_baseline[[col]])))
     pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
     pbc_baseline[[col]] <- factor(pbc_baseline[[col]])
  }
}

# --- 3. Data Splitting ---
set.seed(123)
train_indices <- sample(1:nrow(pbc_baseline), size = floor(0.75 * nrow(pbc_baseline)))
train_data <- pbc_baseline[train_indices, ]
test_data  <- pbc_baseline[-train_indices, ]

# --- 4. Run Multiple Models using RunSurvModels ---
# Assuming RunSurvModels and its dependencies (SurvModel_*, Predict_*) are loaded
cat("Running multiple survival models...\n")

# Select models to run (ensure corresponding SurvModel_* functions exist)
# Example: Run RF (implicitly included), Cox, and glmnet
models_to_run <- c("coxph", "glmnet") # Add others like "xgboost", "gbm", "bart" if available and desired

# RunSurvModels automatically includes RF_Model and RF_Model2 (top vars)
# It uses top N variables (nvars) for the additional models specified in 'models'
all_models_output <- RunSurvModels(
  datatrain = train_data,
  ExpVars = predictor_vars,
  timevar = time_var,
  eventvar = event_var,
  models = models_to_run,
  ntreeRF = 100, # Example RF parameter
  nvars = 15     # Example: Use top 15 vars for Cox, glmnet, etc.
)

# --- 5. Predict with Individual Models ---
# Need prediction functions for each model type (Predict_SurvModel_*)
# We'll predict manually here for clarity, though PredictSurvModels does this internally for ensembling.

cat("\nPredicting with individual models on test data...\n")
test_predictions_list <- list()
test_times <- test_data[[time_var]]
test_events <- test_data[[event_var]]

# Helper function for prediction and evaluation
predict_and_evaluate <- function(model_name, model_obj, predict_func) {
  if (is.null(model_obj)) return(NULL)
  cat("  Predicting with", model_name, "...\n")
  pred_out <- tryCatch(predict_func(model_obj, newdata = test_data),
                       error = function(e) { warning("Pred failed for ", model_name); NULL })
  if (is.null(pred_out)) return(NULL)

  cat("  Evaluating", model_name, "...\n")
  c_index <- tryCatch(timedepConcordance(pred_out$Probs, pred_out$Times, test_times, test_events),
                      error = function(e) { warning("C-index failed for ", model_name); NULL })
  return(list(predictions = pred_out, c_index = c_index))
}

# Predict and evaluate each model that was successfully trained
results_list <- list()
results_list$RF <- predict_and_evaluate("RF", all_models_output$RF_Model, Predict_SurvModel_RF)
results_list$RF2 <- predict_and_evaluate("RF2 (Top Vars)", all_models_output$RF_Model2, Predict_SurvModel_RF)
results_list$Cox <- predict_and_evaluate("CoxPH (Top Vars)", all_models_output$CPH_Model, Predict_SurvModel_Cox)
results_list$glmnet <- predict_and_evaluate("glmnet (Top Vars)", all_models_output$glmnet_Model, Predict_SurvModel_glmnet)
# Add predictions for other models if they were run...

# Filter out models that failed prediction/evaluation
valid_results <- Filter(Negate(is.null), results_list)

# --- 6. Compare Model Performance ---
cat("\n--- Time-dependent C-index Comparison ---\n")

# Extract C-index results for plotting
cindex_objects <- lapply(valid_results, function(res) res$c_index)
names(cindex_objects) <- names(valid_results) # Keep model names

# Plot using pec::plot() if C-index objects are compatible
if (length(cindex_objects) > 0) {
  # Check if objects are of the expected class from pec::cindex
  if (all(sapply(cindex_objects, inherits, "pecCindex"))) {
     plot(cindex_objects, legend.title = "Model", xlab="Time (days)", ylab="Time-dependent C-index")
     title("Model Comparison: Time-dependent Concordance")
  } else {
     warning("Cannot plot C-index directly, objects might not be from pec::cindex.")
     # Manual plotting could be done by extracting AppErr$matrix from each object
     # Example: plot(results_list$RF$c_index$times, results_list$RF$c_index$AppErr$matrix, type='l', col='blue')
     #          lines(results_list$Cox$c_index$times, results_list$Cox$c_index$AppErr$matrix, type='l', col='red') ...
  }
} else {
  cat("No valid C-index results to plot.\n")
}

# --- 7. Calculate Integrated Metrics (Optional) ---
# Assuming integratedC function is available
cat("\n--- Integrated C-index (Example Time Range: 0 to max test time) ---\n")
max_time <- max(test_times)
integration_range <- c(0, max_time)

for (model_name in names(valid_results)) {
  c_index_obj <- valid_results[[model_name]]$c_index
  if (!is.null(c_index_obj)) {
    # Assuming integratedC takes times and scores (AppErr$matrix)
    integrated_c <- tryCatch(
      integratedC(times = c_index_obj$times, scores = c_index_obj$AppErr$matrix, minmax = integration_range),
      error = function(e) { warning("Integrated C failed for ", model_name); NA }
    )
    cat(sprintf("%-15s Integrated C: %.4f\n", model_name, integrated_c))
  }
}

# --- End of Example 03 ---
