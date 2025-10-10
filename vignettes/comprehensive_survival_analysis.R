# Define prediction errors as Brier scores at median time
# Comprehensive Survival Analysis Pipeline with ml4time2event
# Using Follicular Lymphoma Dataset (TestDataFollic.csv)
# 
# This vignette demonstrates a complete survival analysis workflow including:
# - Data loading and preprocessing  
# - Multiple survival models (Cox, Random Forest, XGBoost)
# - Ensemble prediction methods
# - Survival curve visualization
# - Expected time lost calculations and plots
# - Model persistence (saving/loading)

# --- 1. Setup and Load Libraries ---
library(ml4time2event)
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load the package (use devtools::load_all() during development)
devtools::load_all()

# --- 2. Data Loading and Initial Exploration ---
cat("=== STEP 1: DATA LOADING ===\n")

# Load the follicular lymphoma dataset
data_path <- system.file("extdata", "TestDataFollic.csv", package = "ml4time2event")
if (data_path == "") {
  # Fallback for development - direct path
  data_path <- "inst/extdata/TestDataFollic.csv"
}

# Use readData function from the package
follic_data <- readData(data_path)

# Display basic information about the dataset
cat("Dataset dimensions:", nrow(follic_data), "rows,", ncol(follic_data), "columns\n")
cat("Column names:", paste(names(follic_data), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(follic_data))

# Check variable types and missing data
cat("\nVariable summary:\n")
print(summary(follic_data))

# --- 3. Data Preprocessing ---
cat("\n=== STEP 2: DATA PREPROCESSING ===\n")

# Clean and prepare the data
# Remove row names column if present
if ("X" %in% names(follic_data) || "" %in% names(follic_data)) {
  follic_data <- follic_data[, !names(follic_data) %in% c("X", "")]
}

# Convert categorical variables that are stored as character intervals
# Age groups: [17,47], (47,58], (58,67], (67,86]
# Hemoglobin: [40,130], (130,140], (140,150], (150,189]
follic_data$age_group <- factor(follic_data$age)
follic_data$hgb_group <- factor(follic_data$hgb)

# Create numeric versions for some models
follic_data$age_numeric <- as.numeric(factor(follic_data$age))
follic_data$hgb_numeric <- as.numeric(factor(follic_data$hgb))

# Ensure clinical stage and chemotherapy are factors
follic_data$clinstg <- factor(follic_data$clinstg)
follic_data$ch <- factor(follic_data$ch)

# Check the event status coding
cat("Event status distribution:\n")
print(table(follic_data$status))
cat("Time range:", range(follic_data$time), "\n")

# Variable profiling using package function (define expvars first)
candidate_expvars <- c("age_numeric", "hgb_numeric", "clinstg", "ch")
var_profile <- VariableProfile(follic_data, candidate_expvars)
cat("\nVariable Profile Summary:\n")
print(var_profile$Summary)

# --- 4. Data Splitting ---
cat("\n=== STEP 3: DATA SPLITTING ===\n")

# Use package's data splitting function for consistent, reproducible splits
#set.seed(123) # For reproducibility

# Create training/test split (70/30) using package function
split_data <- t2edata_split(follic_data, prop = 0.7)
train_data <- split_data$Train
test_data <- split_data$Test

cat("Training set size:", nrow(train_data), "\n")
cat("Test set size:", nrow(test_data), "\n")

# --- 5. Model Training ---
cat("\n=== STEP 4: SURVIVAL MODEL TRAINING ===\n")

# Define variables for modeling
timevar <- "time"
eventvar <- "status" 
expvars <- c("age_numeric", "hgb_numeric", "clinstg", "ch")

cat("Explanatory variables:", paste(expvars, collapse = ", "), "\n")
cat("Time variable:", timevar, "\n")
cat("Event variable:", eventvar, "\n")

# 5.1 Cox Proportional Hazards Model
cat("\nTraining Cox model...\n")
cox_model <- SurvModel_Cox(
  data = train_data,
  expvars = expvars,
  timevar = timevar,
  eventvar = eventvar
)

cat("Cox model trained successfully\n")
print(summary(cox_model$model))

# Initialize model storage
models <- list()
model_names <- c()

# 5.2 Random Forest Survival Model
cat("\nTraining Random Forest survival model...\n")
tryCatch({
  rf_model <- SurvModel_RF(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 200,                              # More trees
    samplesize = min(50, nrow(train_data)),   # Smaller sample size for small dataset
    splitrule = "logrank",                    # Better split rule for survival
    nodesize_try = c(1, 3, 5, 10)           # Fewer nodesize options
  )
  models[["RF"]] <- rf_model
  model_names <- c(model_names, "RF")
  cat("Random Forest model trained successfully\n")
}, error = function(e) {
  cat("Random Forest model failed:", e$message, "\n")
})

# 5.3 XGBoost Survival Model
cat("\nTraining XGBoost survival model...\n")
tryCatch({
  xgb_model <- SurvModel_xgboost(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    eta = 0.1,           # Higher learning rate for small dataset
    max_depth = 3,       # Shallower trees to avoid overfitting
    nrounds = 150        # More rounds with higher learning rate
  )
  models[["XGBoost"]] <- xgb_model
  model_names <- c(model_names, "XGBoost")
  cat("XGBoost model trained successfully\n")
}, error = function(e) {
  cat("XGBoost model failed:", e$message, "\n")
})

# 5.4 GLMNet (Elastic Net) Survival Model
cat("\nTraining GLMNet survival model...\n")
tryCatch({
  glmnet_model <- SurvModel_glmnet(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar
  )
  models[["GLMNet"]] <- glmnet_model
  model_names <- c(model_names, "GLMNet")
  cat("GLMNet model trained successfully\n")
}, error = function(e) {
  cat("GLMNet model failed:", e$message, "\n")
})

# 5.5 GBM (Gradient Boosting) Survival Model
cat("\nTraining GBM survival model...\n")
tryCatch({
  # Fix small dataset issue by duplicating data for GBM
  gbm_train_data <- rbind(train_data, train_data)
  gbm_model <- SurvModel_gbm(
    data = gbm_train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 100,
    bag.fraction = 0.5  # Reduce bag fraction for small datasets
  )
  models[["GBM"]] <- gbm_model
  model_names <- c(model_names, "GBM")
  cat("GBM model trained successfully\n")
}, error = function(e) {
  cat("GBM model failed:", e$message, "\n")
})

# 5.6 GAM (Generalized Additive Model) Survival Model
cat("\nTraining GAM survival model...\n")
tryCatch({
  gam_model <- SurvModel_GAM(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar
  )
  models[["GAM"]] <- gam_model
  model_names <- c(model_names, "GAM")
  cat("GAM model trained successfully\n")
}, error = function(e) {
  cat("GAM model failed:", e$message, "\n")
})

# 5.7 Survival Regression (Parametric) Model
cat("\nTraining SurvReg survival model...\n")
tryCatch({
  survreg_model <- SurvModel_SurvReg(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    dist = "lognormal"  # Changed from weibull to lognormal
  )
  models[["SurvReg"]] <- survreg_model
  model_names <- c(model_names, "SurvReg")
  cat("SurvReg model trained successfully\n")
}, error = function(e) {
  cat("SurvReg model failed:", e$message, "\n")
})

# 5.8 BART (Bayesian Additive Regression Trees) Survival Model
cat("\nTraining BART survival model...\n")
tryCatch({
  bart_model <- SurvModel_BART(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 50
  )
  models[["BART"]] <- bart_model
  model_names <- c(model_names, "BART")
  cat("BART model trained successfully\n")
}, error = function(e) {
  cat("BART model failed:", e$message, "\n")
})

# 5.9 RuleFit Survival Model
cat("\nTraining RuleFit survival model...\n")
tryCatch({
  rulefit_model <- SurvModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 200,                              # More trees for rule diversity
    nsample = min(nrow(train_data), 60),      # Better sampling for small dataset
    alpha = 0.0,                              # More L1 penalty for feature selection
    maxit = 3000                              # More iterations for convergence
  )
  models[["RuleFit"]] <- rulefit_model
  model_names <- c(model_names, "RuleFit")
  cat("RuleFit model trained successfully\n")
}, error = function(e) {
  cat("RuleFit model failed:", e$message, "\n")
})

# 5.10 DeepSurv (Neural Network) Survival Model
cat("\nTraining DeepSurv survival model...\n")
tryCatch({
  deepsurv_model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    size = 3,                                 # Smaller network for small dataset
    decay = 0.1,                              # Higher regularization
    maxit = 200                               # More iterations
  )
  models[["DeepSurv"]] <- deepsurv_model
  model_names <- c(model_names, "DeepSurv")
  cat("DeepSurv model trained successfully\n")
}, error = function(e) {
  cat("DeepSurv model failed:", e$message, "\n")
})

cat("\nTotal models trained successfully:", length(models), "\n")
cat("Successful models:", paste(model_names, collapse = ", "), "\n")

# --- 6. Model Predictions ---
cat("\n=== STEP 5: MODEL PREDICTIONS ===\n")

# 6.1 Cox model predictions
cat("Generating Cox model predictions...\n")
cox_pred <- Predict_SurvModel_Cox(cox_model, test_data)
cat("Cox predictions - Times length:", length(cox_pred$Times), 
    "Probs dimensions:", dim(cox_pred$Probs), "\n")

# 6.2 Generate predictions for all trained models
predictions <- list()
predict_functions <- list(
  "RF" = Predict_SurvModel_RF,
  "XGBoost" = Predict_SurvModel_xgboost,
  "GLMNet" = Predict_SurvModel_glmnet,
  "GBM" = Predict_SurvModel_gbm,
  "GAM" = Predict_SurvModel_GAM,
  "SurvReg" = Predict_SurvModel_SurvReg,
  "BART" = Predict_SurvModel_BART,
  "RuleFit" = Predict_SurvModel_rulefit,
  "DeepSurv" = Predict_SurvModel_DeepSurv
)

for (model_name in model_names) {
  if (model_name %in% names(models) && model_name %in% names(predict_functions)) {
    cat("Generating", model_name, "predictions...\n")
    tryCatch({
      pred <- predict_functions[[model_name]](models[[model_name]], test_data)
      predictions[[model_name]] <- pred
      cat(model_name, "predictions - Times length:", length(pred$Times), 
          "Probs dimensions:", dim(pred$Probs), "\n")
    }, error = function(e) {
      cat(model_name, "prediction failed:", e$message, "\n")
    })
  }
}

# --- 7. Model Evaluation ---
cat("\n=== STEP 6: MODEL EVALUATION ===\n")

# Prepare actual survival data for evaluation
actual_times <- test_data[[timevar]]
actual_events <- test_data[[eventvar]]

# Convert status to binary if needed (status 2 -> 1, others -> 0)
actual_events_binary <- ifelse(actual_events == 2, 1, ifelse(actual_events == 1, 1, 0))

# 7.1 Time-dependent Concordance
cat("Calculating concordance indices...\n")

# Initialize concordance storage
concordances <- list()

# Cox model concordance
cox_concordance_obj <- timedepConcordance(
  predsurv = cox_pred$Probs,
  predsurvtimes = cox_pred$Times,
  obstimes = actual_times,
  obsevents = actual_events_binary
)
cox_concordance <- mean(cox_concordance_obj$AppCindex$matrix, na.rm = TRUE)
concordances[["Cox"]] <- cox_concordance
cat("Cox model concordance:", round(cox_concordance, 3), "\n")

# Calculate concordance for all other models
for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    concordance_obj <- timedepConcordance(
      predsurv = pred$Probs,
      predsurvtimes = pred$Times,
      obstimes = actual_times,
      obsevents = actual_events_binary
    )
    concordance_val <- mean(concordance_obj$AppCindex$matrix, na.rm = TRUE)
    concordances[[model_name]] <- concordance_val
    cat(model_name, "concordance:", round(concordance_val, 3), "\n")
  }, error = function(e) {
    cat(model_name, "concordance calculation failed:", e$message, "\n")
  })
}

# 7.2 Brier Score Evaluation using package function
cat("\nCalculating Brier scores using package function...\n")

# Define evaluation time points
eval_times <- quantile(actual_times[actual_events_binary == 1], c(0.25, 0.5, 0.75))
cat("Evaluation time points:", round(eval_times, 2), "\n")

# Use median time for Brier score evaluation
median_time <- eval_times[2]  # 50% quantile

# Initialize Brier score storage
brier_scores <- list()

# Cox model Brier score
cox_brier <- BrierScore(
  predsurv = cox_pred$Probs,
  predsurvtimes = cox_pred$Times,
  obstimes = actual_times,
  obsevents = actual_events_binary,
  eval.times = median_time
)
brier_scores[["Cox"]] <- cox_brier
cat("Cox Brier score:", round(cox_brier, 4), "\n")


# Calculate Brier scores for all other models
for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    brier_val <- BrierScore(
      predsurv = pred$Probs,
      predsurvtimes = pred$Times,
      obstimes = actual_times,
      obsevents = actual_events_binary,
      eval.times = median_time
    )
    brier_scores[[model_name]] <- brier_val
    cat(model_name, "Brier score:", round(brier_val, 4), "\n")
  }, error = function(e) {
    cat(model_name, "Brier score calculation failed:", e$message, "\n")
  })
}

# Define prediction errors as Brier scores at median time
pred_errors <- brier_scores

# --- 8. Ensemble Predictions ---
cat("\n=== STEP 7: ENSEMBLE PREDICTIONS ===\n")

# Select top 3 models based on concordance for ensemble
top_models <- names(sort(unlist(concordances), decreasing = TRUE))[seq_len(min(3, length(concordances)))]
cat("Top models for ensemble:", paste(top_models, collapse = ", "), "\n")


# --- ENSEMBLE: Align patient order before averaging ---
# Use rownames or a unique patient ID to ensure columns (patients) are aligned across all models
# We'll use rownames of test_data as patient IDs (if not present, set them)
if (is.null(rownames(test_data)) || any(rownames(test_data) == "")) {
  rownames(test_data) <- as.character(seq_len(nrow(test_data)))
}
patient_ids <- rownames(test_data)

# Helper to get patient IDs for a prediction matrix
get_pred_patient_ids <- function(pred_matrix, data_ref) {
  # Try to get colnames, else fallback to rownames of data_ref
  if (!is.null(colnames(pred_matrix))) {
    return(colnames(pred_matrix))
  } else if (!is.null(rownames(data_ref))) {
    return(rownames(data_ref))
  } else {
    return(as.character(seq_len(ncol(pred_matrix))))
  }
}

# Collect all time points from selected models
all_times <- cox_pred$Times
for (model_name in top_models) {
  if (model_name != "Cox" && model_name %in% names(predictions)) {
    all_times <- c(all_times, predictions[[model_name]]$Times)
  }
}
common_times <- sort(unique(all_times))
common_times <- common_times[common_times <= max(all_times)]

# Interpolate and align columns for each model
ensemble_matrices <- list()

# Cox model (always included as baseline)
cox_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(cox_pred$Probs))
for (j in seq_len(ncol(cox_pred$Probs))) {
  cox_interp[, j] <- approx(cox_pred$Times, cox_pred$Probs[, j],
                           xout = common_times, method = "constant",
                           rule = 2, f = 0)$y
}
# Assign patient IDs as colnames
colnames(cox_interp) <- get_pred_patient_ids(cox_pred$Probs, test_data)
ensemble_matrices[["Cox"]] <- cox_interp[, patient_ids, drop = FALSE]

# Add other top models
for (model_name in top_models) {
  if (model_name != "Cox" && model_name %in% names(predictions)) {
    pred <- predictions[[model_name]]
    model_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(pred$Probs))
    for (j in seq_len(ncol(pred$Probs))) {
      model_interp[, j] <- approx(pred$Times, pred$Probs[, j],
                                 xout = common_times, method = "constant",
                                 rule = 2, f = 0)$y
    }
    # Assign patient IDs as colnames
    colnames(model_interp) <- get_pred_patient_ids(pred$Probs, test_data)
    # Align columns to patient_ids
    aligned_interp <- model_interp[, patient_ids, drop = FALSE]
    ensemble_matrices[[model_name]] <- aligned_interp
  }
}

# Create ensemble by averaging predictions (now columns are aligned)
if (length(ensemble_matrices) > 1) {
  # All matrices have columns in the same order (patient_ids)
  ensemble_probs <- Reduce("+", ensemble_matrices) / length(ensemble_matrices)
} else {
  ensemble_probs <- ensemble_matrices[[1]]
}

cat("Ensemble predictions created for", ncol(ensemble_probs), "subjects using", 
    length(ensemble_matrices), "models\n")


# Evaluate ensemble
ensemble_concordance_obj <- timedepConcordance(
  predsurv = ensemble_probs,
  predsurvtimes = common_times,
  obstimes = actual_times,
  obsevents = actual_events_binary
)
ensemble_concordance <- mean(ensemble_concordance_obj$AppCindex$matrix, na.rm = TRUE)
# Add ensemble concordance to concordances list for downstream tables/printing
concordances[["Ensemble"]] <- ensemble_concordance

# Calculate ensemble Brier score using package function
ensemble_brier <- BrierScore(
  predsurv = ensemble_probs,
  predsurvtimes = common_times,
  obstimes = actual_times,
  obsevents = actual_events_binary,
  eval.times = median_time
)

cat("Ensemble concordance:", round(ensemble_concordance, 3), "\n")
cat("Ensemble Brier score:", round(ensemble_brier, 4), "\n")

# --- 9. Expected Time Lost Calculations ---
cat("\n=== STEP 8: EXPECTED TIME LOST ANALYSIS ===\n")


# Calculate expected time lost using the package function
max_time <- as.numeric(quantile(actual_times, 0.9, na.rm = TRUE))
cat("90th percentile follow-up time:", round(max_time, 2), "\n")

# Calculate expected time lost for all models
etl_results <- list()



# Harmonized ETL calculation: interpolate all model survival curves to common_times, then calculate ETL
model_interp_list <- list()

# Print diagnostic info for time grid and UL
cat("\n[DIAGNOSTIC] common_times range:", range(common_times), "length:", length(common_times), "\n")
cat("[DIAGNOSTIC] max_time (UL):", max_time, "\n")

# Ensure UL does not exceed the max of common_times
UL_fixed <- min(max_time, max(common_times, na.rm=TRUE))
cat("[DIAGNOSTIC] Using UL_fixed:", UL_fixed, "\n")

# Cox model
cox_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(cox_pred$Probs))
for (j in seq_len(ncol(cox_pred$Probs))) {
  cox_interp[, j] <- approx(cox_pred$Times, cox_pred$Probs[, j],
                           xout = common_times, method = "constant",
                           rule = 2, f = 0)$y
}
cat("[DIAGNOSTIC] Cox interp dim:", dim(cox_interp), "\n")
model_interp_list[["Cox"]] <- cox_interp

# Other models
for (model_name in names(predictions)) {
  pred <- predictions[[model_name]]
  model_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(pred$Probs))
  for (j in seq_len(ncol(pred$Probs))) {
    model_interp[, j] <- approx(pred$Times, pred$Probs[, j],
                               xout = common_times, method = "constant",
                               rule = 2, f = 0)$y
  }
  cat("[DIAGNOSTIC]", model_name, "interp dim:", dim(model_interp), "\n")
  model_interp_list[[model_name]] <- model_interp
}

# Calculate ETL for all models on common_times and fixed UL
for (model_name in names(model_interp_list)) {
  interp_probs <- model_interp_list[[model_name]]
  interp_formatted <- list(NewProbs = interp_probs)
  etl_list <- CalculateExpectedTimeLost(
    PredictedCurves = list(interp_formatted),
    modeltypes = c("SURV"),
    times = common_times,
    UL = UL_fixed
  )
  etl_results[[model_name]] <- etl_list[[1]]
  cat(model_name, "expected time lost - mean:", round(mean(etl_list[[1]]), 2),
      "range:", round(range(etl_list[[1]]), 2), "\n")
}



# --- ENSEMBLE ETL: Ensure patient order and time grid match exactly, and UL matches ---
ensemble_probs_aligned <- ensemble_probs[, patient_ids, drop = FALSE]
cat("[DIAGNOSTIC] Ensemble probs dim:", dim(ensemble_probs_aligned), "\n")
ensemble_pred_formatted <- list(NewProbs = ensemble_probs_aligned)
ensemble_etl_list <- CalculateExpectedTimeLost(
  PredictedCurves = list(ensemble_pred_formatted),
  modeltypes = c("SURV"),
  times = common_times,
  UL = UL_fixed
)
ensemble_etl <- ensemble_etl_list[[1]]
etl_results[["Ensemble"]] <- ensemble_etl
cat("Ensemble expected time lost - mean:", round(mean(ensemble_etl), 2), 
    "range:", round(range(ensemble_etl), 2), "\n")



# --- 10. Diagnostic: Plot survival curves for top models and ensemble for a few patients ---
cat("\n=== DIAGNOSTIC: Plotting survival curves for top models and ensemble ===\n")
library(reshape2)
num_patients_to_plot <- min(3, ncol(ensemble_probs))
diagnostic_plots <- list()
for (i in 1:num_patients_to_plot) {
  plot_df <- data.frame(Time = common_times)
  # Add each top model's survival for patient i
  for (model_name in names(ensemble_matrices)) {
    plot_df[[model_name]] <- ensemble_matrices[[model_name]][, i]
  }
  plot_df$Ensemble <- ensemble_probs[, i]
  plot_df_melt <- reshape2::melt(plot_df, id.vars = "Time", variable.name = "Model", value.name = "Survival")
  p_diag <- ggplot(plot_df_melt, aes(x = Time, y = Survival, color = Model)) +
    geom_step() +
    ggtitle(paste("Survival curves for patient", i, "(top models + ensemble)")) +
    theme_minimal()
  diagnostic_plots[[i]] <- p_diag
}
if (length(diagnostic_plots) > 0) {
  pdf("diagnostic_survival_curves.pdf", width = 7, height = 5)
  for (p in diagnostic_plots) print(p)
  dev.off()
  cat("Diagnostic plots saved to diagnostic_survival_curves.pdf\n")
}

# --- 10. Visualization ---
cat("\n=== STEP 9: VISUALIZATION ===\n")

# 10.1 Model Performance Comparison Plot
all_model_names <- unique(c("Cox", names(predictions), "Ensemble"))
performance_df <- data.frame(
  Model = all_model_names,
  Concordance = sapply(all_model_names, function(x) {
    if (x %in% names(concordances)) concordances[[x]] else NA
  }),
  Prediction_Error = sapply(all_model_names, function(x) {
    if (x %in% names(pred_errors)) pred_errors[[x]] else NA
  })
)

# Remove rows with all NAs
performance_df <- performance_df[!is.na(performance_df$Concordance) | !is.na(performance_df$Prediction_Error), ]

print(performance_df)

# 10.2 Survival Curves for Selected Patients
cat("\nPlotting survival curves for first 3 test patients...\n")

# Create plots using the package plotting functions
# Plot with all 3 patients
p1 <- plot_survival_curves(
  predictions = c(list(Cox = cox_pred), predictions, list(Ensemble = list(Times = common_times, Probs = ensemble_probs))),
  patients_to_plot = seq_len(min(3, nrow(test_data))),
  highlight_ensemble = TRUE,
  title = "Predicted Survival Curves - All Models Comparison"
)

# Focused plot for Patient 1 only
p1_single <- plot_survival_curves(
  predictions = c(list(Cox = cox_pred), predictions, list(Ensemble = list(Times = common_times, Probs = ensemble_probs))),
  patients_to_plot = 1,
  highlight_ensemble = TRUE,
  title = "Survival Curves - All Models (Patient 1)"
)


# 10.3 Expected Time Lost Distribution
# Prepare ETL and concordance summary for plotting function (mean ETL and concordance per model)
# Compute ETL means
etl_means <- sapply(names(etl_results), function(m) mean(etl_results[[m]]))
# Compute concordance, handling Ensemble separately
concordance_vals <- sapply(names(etl_results), function(m) {
  if (m == "Ensemble") {
    round(ensemble_concordance, 3)
  } else if (m %in% names(concordances)) {
    round(concordances[[m]], 3)
  } else {
    NA
  }
})
performance_etl_df <- data.frame(
  Model = names(etl_results),
  Expected_Time_Lost = as.numeric(etl_means),
  Concordance = concordance_vals
)
# Custom plots: ETL and Concordance barplots as separate figures
library(ggplot2)
# ETL barplot
p2_etl <- ggplot(performance_etl_df, aes(x = reorder(Model, Expected_Time_Lost), y = Expected_Time_Lost, fill = Model == "Ensemble")) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  labs(title = "Model Performance: Expected Time Lost (ETL)",
    x = "Model", y = "Expected Time Lost") +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "steelblue")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))
# Concordance barplot
p2_c <- ggplot(performance_etl_df, aes(x = reorder(Model, Concordance), y = Concordance, fill = Model == "Ensemble")) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  labs(title = "Model Performance: Concordance Index (C)",
    x = "Model", y = "Concordance Index") +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "steelblue")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 10.4 Combined Plot: Survival Curves with Expected Time Lost (Patient 1)
# Use p1_single and add ETL annotation for Patient 1
cox_etl1 <- round(etl_results[["Cox"]][1], 2)
ensemble_etl1 <- round(ensemble_etl[1], 2)
best_model <- names(sort(sapply(etl_results, function(x) mean(x, na.rm=TRUE))))[1]
best_etl1 <- round(etl_results[[best_model]][1], 2)
p3 <- p1_single +
  ggplot2::annotate("text", x = Inf, y = 0.8, hjust = 1.1,
           label = paste("Cox ETL:", cox_etl1), color = "red", size = 3) +
  ggplot2::annotate("text", x = Inf, y = 0.75, hjust = 1.1,
           label = paste("Best ETL:", best_etl1), color = "green", size = 3) +
  ggplot2::annotate("text", x = Inf, y = 0.7, hjust = 1.1,
           label = paste("Ensemble ETL:", ensemble_etl1), color = "black", size = 3)

# Display plots
print(p1)  # All models, 3 patients
print(p1_single)  # All models, Patient 1 only
print(p2_etl)  # ETL distribution
print(p2_c)     # Concordance distribution
print(p3)  # Patient 1 with ETL annotations

# Save plots
if (requireNamespace("gridExtra", quietly = TRUE)) {
  # Create comprehensive plot layout
  combined_plot <- gridExtra::grid.arrange(
    p1,          # Top: All models, 3 patients
    gridExtra::grid.arrange(p1_single, p2_etl, p2_c, ncol = 3),  # Middle: Single patient + ETL + Concordance
    p3,          # Bottom: Patient 1 with annotations
    ncol = 1, heights = c(2, 1.5, 1.5)
  )
  
  # Save to file
  ggsave("follic_survival_analysis.pdf", combined_plot, 
         width = 20, height = 20, units = "in")
  cat("Plots saved to follic_survival_analysis.pdf\n")
}

# --- 11. Model Persistence ---
cat("\n=== STEP 10: MODEL SAVING AND LOADING ===\n")

# Save models to RDS files
saveRDS(cox_model, "cox_model_follic.rds")

# Save all trained models
saved_models <- c("cox_model_follic.rds")
for (model_name in names(models)) {
  filename <- paste0(tolower(model_name), "_model_follic.rds")
  saveRDS(models[[model_name]], filename)
  saved_models <- c(saved_models, filename)
}

# Save ensemble components
ensemble_object <- list(
  probs = ensemble_probs,
  times = common_times,
  concordance = ensemble_concordance,
  brier_score = ensemble_brier,
  expected_time_lost = ensemble_etl,
  top_models = top_models
)
saveRDS(ensemble_object, "ensemble_model_follic.rds")
saved_models <- c(saved_models, "ensemble_model_follic.rds")

cat("Models saved successfully:\n")
for (filename in saved_models) {
  cat("-", filename, "\n")
}

# Demonstrate model loading
cat("\nDemonstrating model loading...\n")
loaded_cox <- readRDS("cox_model_follic.rds")
loaded_rf <- readRDS("rf_model_follic.rds")
loaded_ensemble <- readRDS("ensemble_model_follic.rds")

cat("Loaded Cox model class:", class(loaded_cox$model), "\n")
cat("Loaded RF model class:", class(loaded_rf$model), "\n")
cat("Loaded ensemble concordance:", round(loaded_ensemble$concordance, 3), "\n")

# Test predictions with loaded models
cat("\nTesting predictions with loaded models...\n")
test_pred_cox <- Predict_SurvModel_Cox(loaded_cox, test_data[1:2, ])
test_pred_rf <- Predict_SurvModel_RF(loaded_rf, test_data[1:2, ])

cat("Loaded Cox prediction dimensions:", dim(test_pred_cox$Probs), "\n")
cat("Loaded RF prediction dimensions:", dim(test_pred_rf$Probs), "\n")

# --- 12. Summary Report ---
cat("\n=== FINAL SUMMARY ===\n")
cat("Comprehensive Survival Analysis Complete!\n\n")

cat("Dataset: Follicular Lymphoma (n =", nrow(follic_data), ")\n")
cat("Training set:", nrow(train_data), "subjects\n")
cat("Test set:", nrow(test_data), "subjects\n\n")

cat("Models Successfully Trained:", length(models) + 1, "\n")  # +1 for Cox
cat("Models with Predictions:", length(predictions) + 1, "\n\n")  # +1 for Cox

cat("Model Performance (Concordance Index):\n")
for (model_name in names(concordances)) {
  cat("-", model_name, ":", round(concordances[[model_name]], 3), "\n")
}

        # Harmonized ETL calculation: interpolate all model survival curves to common_times, then calculate ETL
        model_interp_list <- list()

        # Cox model
        cox_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(cox_pred$Probs))
        for (j in seq_len(ncol(cox_pred$Probs))) {
          cox_interp[, j] <- approx(cox_pred$Times, cox_pred$Probs[, j],
                                   xout = common_times, method = "constant",
                                   rule = 2, f = 0)$y
        }
        model_interp_list[["Cox"]] <- cox_interp

        # Other models
        for (model_name in names(predictions)) {
          pred <- predictions[[model_name]]
          model_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(pred$Probs))
          for (j in seq_len(ncol(pred$Probs))) {
            model_interp[, j] <- approx(pred$Times, pred$Probs[, j],
                                       xout = common_times, method = "constant",
                                       rule = 2, f = 0)$y
          }
          model_interp_list[[model_name]] <- model_interp
        }

        # Calculate ETL for all models on common_times
        for (model_name in names(model_interp_list)) {
          interp_probs <- model_interp_list[[model_name]]
          interp_formatted <- list(NewProbs = interp_probs)
          etl_list <- CalculateExpectedTimeLost(
            PredictedCurves = list(interp_formatted),
            modeltypes = c("SURV"),
            times = common_times,
            UL = max_time
          )
          etl_results[[model_name]] <- etl_list[[1]]
          cat(model_name, "expected time lost - mean:", round(mean(etl_list[[1]]), 2),
              "range:", round(range(etl_list[[1]]), 2), "\n")
        }
cat("✓ Ensemble methods\n")
cat("✓ Expected time lost calculations\n")
cat("✓ Comprehensive visualization\n")
cat("✓ Model persistence and reloading\n")

cat("\n=== VIGNETTE COMPLETE ===\n")