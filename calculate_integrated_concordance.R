# Calculate Integrated Concordance for All CR Models
# Using saved models from the comprehensive competing risks analysis

# Load the package
devtools::load_all()

library(survival)
library(dplyr)

# Load the BMT competing risks dataset
data_path <- system.file("extdata", "bmtcrr_competing_risks.csv", package = "ml4time2event")
bmt_data <- read.csv(data_path)

# Define variables
timevar <- "ftime"
eventvar <- "Status"
expvars <- c("Sex", "D", "Phase", "Age", "Source")

# Split data (same as in vignette)
set.seed(123)
train_indices <- sample(1:nrow(bmt_data), 0.7 * nrow(bmt_data))
train_data <- bmt_data[train_indices, ]
test_data <- bmt_data[-train_indices, ]

# Load all saved models
cat("Loading saved competing risks models...\n")
models <- list()
models[["Cox"]] <- readRDS("cox_crmodel_bmt.rds")
models[["FineGray"]] <- readRDS("finegray_crmodel_bmt.rds")
models[["RF"]] <- readRDS("rf_crmodel_bmt.rds")
models[["XGBoost"]] <- readRDS("xgboost_crmodel_bmt.rds")
models[["GAM"]] <- readRDS("gam_crmodel_bmt.rds")
models[["BART"]] <- readRDS("bart_crmodel_bmt.rds")
models[["DeepSurv"]] <- readRDS("deepsurv_crmodel_bmt.rds")
models[["RuleFit"]] <- readRDS("rulefit_crmodel_bmt.rds")
models[["SurvReg"]] <- readRDS("survreg_crmodel_bmt.rds")

# Load ensemble separately (it already contains predictions)
ensemble_model <- readRDS("ensemble_crmodel_bmt.rds")

cat("Loaded", length(models), "individual models + 1 ensemble model\n")

# Define prediction functions (excluding ensemble)
predict_functions <- list(
  "Cox" = Predict_CRModel_Cox,
  "FineGray" = Predict_CRModel_FineGray,
  "RF" = Predict_CRModel_RF,
  "XGBoost" = Predict_CRModel_xgboost,
  "GAM" = Predict_CRModel_GAM,
  "BART" = Predict_CRModel_BART,
  "DeepSurv" = Predict_CRModel_DeepSurv,
  "RuleFit" = Predict_CRModel_rulefit,
  "SurvReg" = Predict_CRModel_SurvReg
)

# Generate predictions for individual models
cat("\nGenerating predictions for individual models...\n")
predictions <- list()

for (model_name in names(models)) {
  if (model_name %in% names(predict_functions)) {
    cat("Generating", model_name, "predictions...\n")
    tryCatch({
      pred <- predict_functions[[model_name]](models[[model_name]], test_data)
      predictions[[model_name]] <- pred
      cat(model_name, "predictions - Times length:", length(pred$Times),
          "CIF dimensions:", dim(pred$CIF), "\n")
    }, error = function(e) {
      cat(model_name, "prediction failed:", e$message, "\n")
    })
  }
}

# Add ensemble predictions (already stored in the model)
predictions[["Ensemble"]] <- list(
  CIF = ensemble_model$cif,
  Times = ensemble_model$times
)
cat("Ensemble predictions - Times length:", length(ensemble_model$times),
    "CIF dimensions:", dim(ensemble_model$cif), "\n")

# Prepare survival object for evaluation
actual_times <- test_data[[timevar]]
actual_events <- test_data[[eventvar]]
surv_obj <- Surv(actual_times, actual_events, type = "mstate")

# Calculate integrated concordance for all models (Event 1)
cat("\n=== INTEGRATED CONCORDANCE RESULTS (Event 1) ===\n")
integrated_concordances <- list()

for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]

    # integratedConcordanceCR expects predictions as matrix (rows=observations, cols=times)
    cif_matrix <- t(pred$CIF)  # Transpose: observations x times

    # Calculate integrated concordance for competing risks (event 1)
    integrated_conc <- integratedConcordanceCR(
      SurvObj = surv_obj,
      Predictions = cif_matrix,
      eval.times = pred$Times,
      cause = 1,
      TestMat = test_data[, expvars]
    )

    integrated_concordances[[model_name]] <- integrated_conc
    cat(model_name, "integrated concordance (Event 1):", round(integrated_conc, 3), "\n")
  }, error = function(e) {
    cat(model_name, "integrated concordance calculation failed:", e$message, "\n")
  })
}

# Display summary
cat("\n=== SUMMARY: INTEGRATED CONCORDANCE INDEX (Event 1) ===\n")
cat("Test set size:", nrow(test_data), "patients\n")
cat("Time range:", round(range(actual_times), 1), "months\n")
cat("Event distribution:\n")
event_table <- table(actual_events)
cat("  Censored (0):", event_table["0"], "(", round(100 * event_table["0"] / nrow(test_data), 1), "%)\n")
cat("  Relapse (1):", event_table["1"], "(", round(100 * event_table["1"] / nrow(test_data), 1), "%)\n")
cat("  TRM (2):", event_table["2"], "(", round(100 * event_table["2"] / nrow(test_data), 1), "%)\n\n")

cat("Integrated Concordance Index Results (higher is better):\n")
for (model_name in names(integrated_concordances)) {
  if (!is.na(integrated_concordances[[model_name]])) {
    cat(sprintf("  %-10s: %.3f\n", model_name, integrated_concordances[[model_name]]))
  }
}

# Find best performing model
valid_concs <- integrated_concordances[!is.na(unlist(integrated_concordances))]
if (length(valid_concs) > 0) {
  best_model <- names(valid_concs)[which.max(unlist(valid_concs))]
  best_score <- max(unlist(valid_concs))
  cat("\n✓ Best performing model:", best_model, "with integrated concordance:", round(best_score, 3), "\n")
}