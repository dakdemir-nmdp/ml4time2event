# Example 04: Ensemble Prediction using ml4time2event

# --- 1. Load Libraries ---
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
# Install dependencies for models used in the ensemble (e.g., glmnet, RF, etc.)
# if (!requireNamespace("glmnet", quietly = TRUE)) { install.packages("glmnet") }
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

# --- 4. Run Multiple Models (Prerequisite for Ensemble) ---
# We need the output from RunSurvModels, similar to Example 03.
# Assuming RunSurvModels and its dependencies are loaded.
cat("Running multiple survival models to generate inputs for ensemble...\n")
models_to_run <- c("coxph", "glmnet") # Example models
all_models_output <- RunSurvModels(
  datatrain = train_data,
  ExpVars = predictor_vars,
  timevar = time_var,
  eventvar = event_var,
  models = models_to_run,
  ntreeRF = 100,
  nvars = 15
)

# --- 5. Generate Ensemble Predictions using PredictSurvModels ---
# Assuming PredictSurvModels, survprobMatInterpolator, and survprobMatListAveraging are loaded.
cat("\nGenerating ensemble predictions on test data...\n")

# Define prediction times (e.g., quantiles of event times or specific points of interest)
test_times_observed <- test_data[[time_var]][test_data[[event_var]] == 1]
prediction_times <- quantile(test_times_observed, probs = seq(0, 1, 0.1), na.rm = TRUE)
# Ensure prediction times are unique and sorted
prediction_times <- sort(unique(c(0, prediction_times)))

# Call PredictSurvModels
ensemble_predictions_output <- PredictSurvModels(
  models = all_models_output,
  newdata = test_data,
  newtimes = prediction_times
)

# --- 6. Evaluate Ensemble Performance ---
cat("\nEvaluating ensemble model performance...\n")
test_times <- test_data[[time_var]]
test_events <- test_data[[event_var]]

if (!is.null(ensemble_predictions_output$NewProbs)) {
  # Calculate time-dependent C-index for the ensemble
  ensemble_c_index <- timedepConcordance(
    predsurv = ensemble_predictions_output$NewProbs,
    predsurvtimes = prediction_times, # Use the times for which ensemble probs were calculated
    obstimes = test_times,
    obsevents = test_events,
    ctimes = prediction_times # Evaluate C-index at these specific times
  )

  cat("\n--- Ensemble Time-dependent C-index Summary ---\n")
  print(summary(ensemble_c_index$AppErr$matrix))

  # Plot ensemble C-index
  plot(ensemble_c_index)
  title("Time-dependent Concordance Index (Ensemble Model)")

  # Calculate Integrated C-index for the ensemble
  # Assuming integratedC is available
  integration_range <- c(min(prediction_times[prediction_times>0]), max(prediction_times)) # Example range
  if(diff(integration_range) > 0) {
      ensemble_integrated_c <- integratedC(
        times = ensemble_c_index$times,
        scores = ensemble_c_index$AppErr$matrix,
        minmax = integration_range
      )
      cat(sprintf("\nEnsemble Integrated C-index (%.1f - %.1f days): %.4f\n",
                  integration_range[1], integration_range[2], ensemble_integrated_c))
  } else {
      cat("\nCannot calculate Integrated C-index: Invalid integration range.\n")
  }


  # Calculate Brier Score for the ensemble (Optional)
  # Assuming BrierScore and integratedBrier are available
  cat("\nCalculating Brier Score for Ensemble...\n")
  ensemble_brier <- BrierScore(
      predsurv = ensemble_predictions_output$NewProbs,
      predsurvtimes = prediction_times,
      obstimes = test_times,
      obsevents = test_events
  )
  # Plot Brier score
  plot(ensemble_brier)
  title("Time-dependent Brier Score (Ensemble Model)")

  # Calculate Integrated Brier Score (IBS)
  ibs <- integratedBrier(
      times = ensemble_brier$times,
      scores = ensemble_brier$score$model, # Accessing the score correctly
      minmax = integration_range
  )
   cat(sprintf("\nEnsemble Integrated Brier Score (%.1f - %.1f days): %.4f\n",
                  integration_range[1], integration_range[2], ibs))


} else {
  cat("\nEnsemble prediction failed. Cannot evaluate.\n")
}

# --- 7. Compare Ensemble with Individual Models (Optional) ---
# You could retrieve individual model C-indices (calculated in Example 03 or recalculated here)
# and plot them alongside the ensemble C-index for comparison.

# --- End of Example 04 ---
