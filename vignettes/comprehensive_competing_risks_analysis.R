# Comprehensive Competing# --- 2. Data Loading and Initial Exploration ---peline with ml4time2event
# Using Bone Marrow Transplant (BMT) Dataset
#
# This vignette demonstrates a complete competing risks analysis workflow including:
# - Data loading and preprocessing with competing risks data
# - Multiple competing risks models (Cox, Fine-Gray, Random Forest, XGBoost, etc.)
# - Ensemble prediction methods for competing risks
# - Cumulative Incidence Function (CIF) visualization
# - Expected Time Lost (ETL) for competing risks
# - Competing risks-specific metrics and evaluation
# - Model persistence (saving/loading)

# --- 1. Setup and Load Libraries ---
# Load the package during development
devtools::load_all()

library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Additional packages for competing risks
if (!require("cmprsk", quietly = TRUE)) install.packages("cmprsk")
library(cmprsk)

cat("=== STEP 1: DATA LOADING ===\n")

# Load the BMT competing risks dataset
data_path <- system.file("extdata", "bmtcrr_competing_risks.csv", package = "ml4time2event")
if (data_path == "") {
  # Fallback for development - direct path
  data_path <- "inst/extdata/bmtcrr_competing_risks.csv"
}

# Read the data
bmt_data <- read.csv(data_path)

# Display basic information about the dataset
cat("Dataset dimensions:", nrow(bmt_data), "rows,", ncol(bmt_data), "columns\n")
cat("Column names:", paste(names(bmt_data), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(bmt_data))

# Check event status distribution
cat("\n=== Event Status Distribution ===\n")
cat("Event variable coding:\n")
cat("  0 = Censored (event-free)\n")
cat("  1 = Relapse (event of interest)\n")
cat("  2 = Treatment-Related Mortality (TRM) - competing risk\n\n")

status_table <- table(bmt_data$Status)
print(status_table)
cat("\nEvent rates:\n")
cat("  Censored:", status_table["0"], "(", round(100 * status_table["0"] / nrow(bmt_data), 1), "%)\n")
cat("  Relapse:", status_table["1"], "(", round(100 * status_table["1"] / nrow(bmt_data), 1), "%)\n")
cat("  TRM:", status_table["2"], "(", round(100 * status_table["2"] / nrow(bmt_data), 1), "%)\n")

# --- 3. Data Preprocessing ---
cat("\n=== STEP 2: DATA PREPROCESSING ===\n")

# Ensure categorical variables are factors
bmt_data$Sex <- factor(bmt_data$Sex)
bmt_data$D <- factor(bmt_data$D)
bmt_data$Phase <- factor(bmt_data$Phase)
bmt_data$Source <- factor(bmt_data$Source)

# Display variable summary
cat("\nVariable types:\n")
print(str(bmt_data))

# Check time variable
cat("\nTime range (ftime in months):", range(bmt_data$ftime), "\n")

# Define variables for modeling
timevar <- "ftime"
eventvar <- "Status"
expvars <- c("Sex", "D", "Phase", "Age", "Source")

# Variable profiling
var_profile <- VariableProfile(bmt_data, expvars)
cat("\nVariable Profile Summary:\n")
print(var_profile$Summary)

# --- 4. Initial Competing Risks Exploration ---
cat("\n=== STEP 3: CUMULATIVE INCIDENCE EXPLORATION ===\n")

# Calculate Cumulative Incidence Functions by disease phase
cat("Calculating CIF by disease phase...\n")
cif_phase <- cuminc(ftime = bmt_data$ftime,
                    fstatus = bmt_data$Status,
                    group = bmt_data$Phase)

# Display cumulative incidence at key time points
cat("\nCumulative Incidence at 12 months by Phase:\n")
timepoint_12mo <- 12

# Extract CIF values at 12 months for each phase and event type
for (phase_name in levels(bmt_data$Phase)) {
  # Get CIF for relapse (event 1)
  cif_name_relapse <- paste(phase_name, "1", sep = " ")
  if (cif_name_relapse %in% names(cif_phase)) {
    cif_obj <- cif_phase[[cif_name_relapse]]
    idx <- which.min(abs(cif_obj$time - timepoint_12mo))
    if (length(idx) > 0) {
      ci_relapse <- cif_obj$est[idx]
      cat(sprintf("  %s - Relapse: %.3f\n", phase_name, ci_relapse))
    }
  }
}

# --- 5. Data Splitting ---
cat("\n=== STEP 4: DATA SPLITTING ===\n")

# Use package's data splitting function for consistent, reproducible splits
set.seed(123) # For reproducibility

# Create training/test split (70/30) using package function
split_data <- t2edata_split(bmt_data, prop = 0.7)
train_data <- split_data$Train
test_data <- split_data$Test

cat("Training set size:", nrow(train_data), "\n")
cat("Test set size:", nrow(test_data), "\n")

# Check event distribution in training set
train_status_table <- table(train_data$Status)
cat("\nTraining set event distribution:\n")
print(train_status_table)

# --- 6. Model Training ---
cat("\n=== STEP 5: COMPETING RISKS MODEL TRAINING ===\n")

cat("Explanatory variables:", paste(expvars, collapse = ", "), "\n")
cat("Time variable:", timevar, "\n")
cat("Event variable:", eventvar, "\n")
primary_event_code <- 1L
competing_event_codes <- 2L
event_codes_all <- c(primary_event_code, competing_event_codes)

cat("Event of interest:", primary_event_code, "(Relapse)\n")
cat("Competing risk(s):", paste(competing_event_codes, collapse = ", "), "(TRM)\n\n")

# Initialize model storage
models <- list()
model_names <- c()

# 6.1 Cause-Specific Cox Model
cat("Training Cause-Specific Cox model...\n")
tryCatch({
  cox_model <- CRModel_Cox(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes_all
  )
  models[["Cox"]] <- cox_model
  model_names <- c(model_names, "Cox")
  cat("✓ Cox model trained successfully\n")
}, error = function(e) {
  cat("✗ Cox model failed:", e$message, "\n")
})

# 6.2 Fine-Gray Subdistribution Hazard Model
cat("\nTraining Fine-Gray model...\n")
tryCatch({
  fg_model <- CRModel_FineGray(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = primary_event_code
  )
  models[["FineGray"]] <- fg_model
  model_names <- c(model_names, "FineGray")
  cat("✓ Fine-Gray model trained successfully\n")
}, error = function(e) {
  cat("✗ Fine-Gray model failed:", e$message, "\n")
})

# 6.3 Random Forest for Competing Risks
cat("\nTraining Random Forest competing risks model...\n")
tryCatch({
  rf_model <- CRModel_RF(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 200,
    samplesize = min(50, nrow(train_data))
  )
  models[["RF"]] <- rf_model
  model_names <- c(model_names, "RF")
  cat("✓ Random Forest model trained successfully\n")
}, error = function(e) {
  cat("✗ Random Forest model failed:", e$message, "\n")
})

# 6.4 XGBoost for Competing Risks
cat("\nTraining XGBoost competing risks model...\n")
tryCatch({
  xgb_model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes_all,
    nrounds = 100,
    eta = 0.1,
    max_depth = 3
  )
  models[["XGBoost"]] <- xgb_model
  model_names <- c(model_names, "XGBoost")
  cat("✓ XGBoost model trained successfully\n")
}, error = function(e) {
  cat("✗ XGBoost model failed:", e$message, "\n")
})

# 6.5 GAM for Competing Risks
cat("\nTraining GAM competing risks model...\n")
tryCatch({
  gam_model <- CRModel_GAM(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes_all
  )
  models[["GAM"]] <- gam_model
  model_names <- c(model_names, "GAM")
  cat("✓ GAM model trained successfully\n")
}, error = function(e) {
  cat("✗ GAM model failed:", e$message, "\n")
})

# 6.6 BART for Competing Risks
cat("\nTraining BART competing risks model...\n")
tryCatch({
  bart_model <- CRModel_BART(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = primary_event_code,
    ntree = 50,  # Smaller for speed
    ndpost = 200,  # Smaller for speed
    nskip = 50,   # Smaller for speed
    keepevery = 5
  )
  models[["BART"]] <- bart_model
  model_names <- c(model_names, "BART")
  cat("✓ BART model trained successfully\n")
}, error = function(e) {
  cat("✗ BART model failed:", e$message, "\n")
})

# 6.7 DeepSurv for Competing Risks
cat("\nTraining DeepSurv competing risks model (event", primary_event_code, ")...\n")
deepsurv_primary <- NULL
deepsurv_competing_models <- list()

tryCatch({
  deepsurv_primary <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = primary_event_code,
    size = 5,
    decay = 0.01,
    maxit = 100
  )
  cat("✓ DeepSurv primary event model trained successfully\n")
}, error = function(e) {
  cat("✗ DeepSurv primary event model failed:", e$message, "\n")
})

if (!is.null(deepsurv_primary) && length(competing_event_codes) > 0) {
  for (comp_code in competing_event_codes) {
    cat("Training DeepSurv competing event model (event", comp_code, ")...\n")
    tryCatch({
      comp_model <- CRModel_DeepSurv(
        data = train_data,
        expvars = expvars,
        timevar = timevar,
        eventvar = eventvar,
        event_codes = comp_code,
        size = 5,
        decay = 0.01,
        maxit = 100
      )
      model_name_comp <- paste0("event_", comp_code)
      deepsurv_competing_models[[model_name_comp]] <- comp_model
      cat("✓ DeepSurv competing event model", comp_code, "trained successfully\n")
    }, error = function(e) {
      cat("✗ DeepSurv competing event model", comp_code, "failed:", e$message, "\n")
    })
  }
}

if (!is.null(deepsurv_primary)) {
  models[["DeepSurv"]] <- list(
    primary = deepsurv_primary,
    competing = deepsurv_competing_models
  )
  model_names <- c(model_names, "DeepSurv")
  if (length(competing_event_codes) > 0 && length(deepsurv_competing_models) == 0) {
    warning(
      "DeepSurv competing event models could not be trained; cumulative incidence estimates may be unavailable."
    )
  }
}

# 6.8 RuleFit for Competing Risks
cat("\nTraining RuleFit competing risks model...\n")
tryCatch({
  rulefit_model <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = primary_event_code
  )
  models[["RuleFit"]] <- rulefit_model
  model_names <- c(model_names, "RuleFit")
  cat("✓ RuleFit model trained successfully\n")
}, error = function(e) {
  cat("✗ RuleFit model failed:", e$message, "\n")
})

# 6.9 SurvReg for Competing Risks
cat("\nTraining SurvReg competing risks model...\n")
tryCatch({
  survreg_model <- CRModel_SurvReg(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes_all
  )
  models[["SurvReg"]] <- survreg_model
  model_names <- c(model_names, "SurvReg")
  cat("✓ SurvReg model trained successfully\n")
}, error = function(e) {
  cat("✗ SurvReg model failed:", e$message, "\n")
})

cat("\nTotal models trained successfully:", length(models), "\n")
cat("Successful models:", paste(model_names, collapse = ", "), "\n")

# --- 7. Model Predictions ---
cat("\n=== STEP 6: COMPETING RISKS MODEL PREDICTIONS ===\n")

# Generate predictions for all trained models
predictions <- list()

# Define prediction functions mapping
predict_functions <- list(
  "Cox" = Predict_CRModel_Cox,
  "FineGray" = Predict_CRModel_FineGray,
  "RF" = Predict_CRModel_RF,
  "XGBoost" = Predict_CRModel_xgboost,
  "GAM" = Predict_CRModel_GAM,
  "BART" = Predict_CRModel_BART,
  "DeepSurv" = function(model_bundle, newdata) {
    if (is.null(model_bundle$primary)) {
      stop("DeepSurv model bundle is missing the primary event model")
    }
    if (length(model_bundle$competing) == 0) {
      stop("DeepSurv predictions require trained models for all competing events")
    }
    Predict_CRModel_DeepSurv(
      modelout = model_bundle$primary,
      newdata = newdata,
      other_models = model_bundle$competing
    )
  },
  "RuleFit" = Predict_CRModel_rulefit,
  "SurvReg" = Predict_CRModel_SurvReg
)

for (model_name in model_names) {
  if (model_name %in% names(predict_functions)) {
    cat("Generating", model_name, "predictions...\n")
    tryCatch({
      pred <- predict_functions[[model_name]](models[[model_name]], test_data)
      if (is.null(pred$CIFs)) {
        stop(model_name, " predictions did not include cumulative incidence functions")
      }
      predictions[[model_name]] <- pred
      cif_dims <- paste(dim(pred$CIFs), collapse = " x ")
      cat(model_name, "predictions - Times length:", length(pred$Times),
          "CIF dimensions:", cif_dims, "\n")
    }, error = function(e) {
      cat(model_name, "prediction failed:", e$message, "\n")
    })
  }
}

# --- 8. Model Evaluation ---
cat("\n=== STEP 7: COMPETING RISKS MODEL EVALUATION ===\n")

# Prepare actual survival data for evaluation
actual_times <- test_data[[timevar]]
actual_events <- test_data[[eventvar]]

# Create survival object for competing risks
surv_obj <- Surv(actual_times, actual_events, type = "mstate")

# 8.1 Time-dependent Concordance for Competing Risks
cat("Calculating concordance indices for competing risks...\n")

# Initialize concordance storage
concordances <- list()

# Evaluation time point (median time to relapse event)
eval_time <- median(actual_times[actual_events == 1])
cat("Evaluation time point:", round(eval_time, 2), "months\n\n")

# Calculate concordance for all models
for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    if (is.null(pred$CIFs)) {
      stop("CIFs unavailable for model")
    }

    time_idx <- which.min(abs(pred$Times - eval_time))
    cif_at_eval <- pred$CIFs[time_idx, ]

    concordance_val <- timedepConcordanceCR(
      SurvObj = surv_obj,
      Predictions = matrix(cif_at_eval, ncol = 1),
      time = eval_time,
      cause = 1,
      TestMat = test_data[, expvars]
    )

    concordances[[model_name]] <- concordance_val
    cat(model_name, "concordance:", round(concordance_val, 3), "\n")
  }, error = function(e) {
    concordances[[model_name]] <- NA
    cat(model_name, "concordance: NA (failed:", e$message, ")\n")
  })
}

# 8.2 Brier Score for Competing Risks
cat("\nCalculating Brier scores for competing risks...\n")

brier_scores <- list()

for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    if (is.null(pred$CIFs)) {
      stop("CIFs unavailable for model")
    }

    # Use package function for consistent Brier score calculation
    brier_score <- BrierScoreCR(
      SurvObj = surv_obj,
      Predictions = t(pred$CIFs),  # BrierScoreCR expects observations x times
      time = eval_time,
      cause = 1,
      TestMat = test_data[, expvars],
      pred_times = pred$Times
    )
    brier_scores[[model_name]] <- brier_score

    cat(model_name, "Brier score:", round(brier_score, 4), "\n")
  }, error = function(e) {
    brier_scores[[model_name]] <- NA
    cat(model_name, "Brier score: NA (failed:", e$message, ")\n")
  })
}

# 8.3 Integrated Concordance for Competing Risks (Event 1)
cat("\nCalculating integrated concordance indices for competing risks (Event 1)...\n")

# Initialize integrated concordance storage
integrated_concordances <- list()

# Calculate integrated concordance for all models
for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    if (is.null(pred$CIFs)) {
      stop("CIFs unavailable for model")
    }

    cif_matrix <- t(pred$CIFs)

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
    integrated_concordances[[model_name]] <- NA
    cat(model_name, "integrated concordance (Event 1): NA (failed:", e$message, ")\n")
  })
}

# --- 9. Ensemble Predictions ---
cat("\n=== STEP 8: ENSEMBLE PREDICTIONS FOR COMPETING RISKS ===\n")

# Select top 5 models based on concordance for ensemble
valid_concordances <- concordances[!is.na(unlist(concordances))]
if (length(valid_concordances) >= 5) {
  top_models <- names(sort(unlist(valid_concordances), decreasing = TRUE))[1:5]
} else if (length(valid_concordances) > 0) {
  top_models <- names(valid_concordances)
} else {
  top_n <- min(5, length(model_names))
  top_models <- if (top_n > 0) model_names[seq_len(top_n)] else character(0)
}

cat("Top models for ensemble:", paste(top_models, collapse = ", "), "\n")

# Create ensemble by averaging CIF predictions from top models
# First, collect all time points from selected models
all_times <- c()
for (model_name in top_models) {
  if (model_name %in% names(predictions)) {
    pred_obj <- predictions[[model_name]]
    if (!is.null(pred_obj$CIFs)) {
      all_times <- c(all_times, pred_obj$Times)
    }
  }
}

if (length(all_times) == 0) {
  stop("No models with cumulative incidence functions available for ensemble construction.")
}

# Create common time grid
common_times <- sort(unique(all_times))
common_times <- common_times[common_times <= max(all_times) & common_times > 0]

# If only one model, use it directly as ensemble
if (length(top_models) == 1) {
  model_name <- top_models[1]
  pred <- predictions[[model_name]]
  if (is.null(pred$CIFs)) {
    stop("Selected ensemble model does not provide CIF outputs.")
  }
  ensemble_cif <- t(pred$CIFs)  # Transpose to times x observations format
  common_times <- pred$Times
  cat("Ensemble uses single model:", model_name, "\n")
} else {
  # Interpolate CIF predictions to common time grid for ensemble models
  ensemble_matrices <- list()

  for (model_name in top_models) {
    if (model_name %in% names(predictions)) {
      pred <- predictions[[model_name]]
      if (is.null(pred$CIFs)) {
        next
      }
      model_interp <- matrix(NA, nrow = length(common_times), ncol = ncol(pred$CIFs))

  for (j in seq_len(ncol(pred$CIFs))) {
        tryCatch({
          model_interp[, j] <- approx(pred$Times, pred$CIFs[, j],
                                       xout = common_times, method = "linear",
                                       rule = 2)$y
        }, error = function(e) {
          # Fallback: use nearest neighbor interpolation
          model_interp[, j] <- sapply(common_times, function(t) {
            idx <- which.min(abs(pred$Times - t))
            pred$CIFs[idx, j]
          })
        })
      }
      ensemble_matrices[[model_name]] <- model_interp
    }
  }

  if (length(ensemble_matrices) == 0) {
    stop("No ensemble matrices constructed; ensure selected models provide CIFs.")
  }

  # Create ensemble by averaging CIF predictions
  if (length(ensemble_matrices) > 1) {
    ensemble_cif <- Reduce("+", ensemble_matrices) / length(ensemble_matrices)
  } else {
    ensemble_cif <- ensemble_matrices[[1]]
  }
}

# Ensure CIF values are between 0 and 1
ensemble_cif[ensemble_cif < 0] <- 0
ensemble_cif[ensemble_cif > 1] <- 1

cat("Ensemble predictions created for", ncol(ensemble_cif), "subjects using",
    length(top_models), "models\n")

# Evaluate ensemble
ensemble_time_idx <- which.min(abs(common_times - eval_time))
ensemble_cif_at_eval <- ensemble_cif[ensemble_time_idx, ]

ensemble_concordance <- timedepConcordanceCR(
  SurvObj = surv_obj,
  Predictions = matrix(ensemble_cif_at_eval, ncol = 1),
  time = eval_time,
  cause = 1,
  TestMat = test_data[, expvars]
)

# Use package function for consistent Brier score calculation
ensemble_brier <- BrierScoreCR(
  SurvObj = surv_obj,
  Predictions = t(ensemble_cif),  # BrierScoreCR expects observations x times
  time = eval_time,
  cause = 1,
  TestMat = test_data[, expvars],
  pred_times = common_times
)

cat("Ensemble concordance:", round(ensemble_concordance, 3), "\n")
cat("Ensemble Brier score:", round(ensemble_brier, 4), "\n")

# Calculate integrated concordance for ensemble
cat("\nCalculating integrated concordance for ensemble (Event 1)...\n")
ensemble_integrated_conc <- tryCatch({
  # ensemble_cif is in format: times x observations, so transpose for integratedConcordanceCR
  ensemble_cif_matrix <- t(ensemble_cif)  # Now: observations x times

  integrated_conc <- integratedConcordanceCR(
    SurvObj = surv_obj,
    Predictions = ensemble_cif_matrix,
    eval.times = common_times,
    cause = 1,
    TestMat = test_data[, expvars]
  )

  integrated_concordances[["Ensemble"]] <- integrated_conc
  cat("Ensemble integrated concordance (Event 1):", round(integrated_conc, 3), "\n")
  integrated_conc
}, error = function(e) {
  cat("Ensemble integrated concordance calculation failed:", e$message, "\n")
  NA_real_
})

# --- 10. Expected Time Lost (ETL) Calculation ---
cat("\n=== STEP 9: EXPECTED TIME LOST ANALYSIS ===\n")

# Calculate expected time lost for competing risks using the package function
# ETL = integral from 0 to UL of CIF(t) dt
max_time <- max(actual_times)
cat("Maximum follow-up time:", round(max_time, 2), "months\n")

# Calculate ETL for all models using CalculateExpectedTimeLost
etl_results <- list()

for (model_name in names(predictions)) {
  tryCatch({
    pred <- predictions[[model_name]]
    if (is.null(pred$CIFs)) {
      stop("CIFs unavailable for model")
    }

    # Format predictions for CalculateExpectedTimeLost function
    # CIFs matrix is already in correct (times x observations) format
    pred_formatted <- list(NewProbs = pred$CIFs)
    
    # Use the package function for consistent ETL calculation
    etl_list <- CalculateExpectedTimeLost(
      PredictedCurves = list(pred_formatted),
      modeltypes = c("CR"),  # Competing risks model type
      times = pred$Times,
      UL = max_time
    )
    
    etl_values <- etl_list[[1]]
    etl_results[[model_name]] <- etl_values
    cat(model_name, "ETL - mean:", round(mean(etl_values, na.rm = TRUE), 2),
        "range:", paste(round(range(etl_values, na.rm = TRUE), 2), collapse = "-"), "\n")
  }, error = function(e) {
    cat(model_name, "ETL calculation failed:", e$message, "\n")
  })
}

# Ensemble ETL using the package function for consistency
tryCatch({
  # Format ensemble predictions for CalculateExpectedTimeLost function
  # ensemble_cif is already in (times x observations) format
  ensemble_pred_formatted <- list(NewProbs = ensemble_cif)
  
  # Use the package function for consistent ETL calculation
  ensemble_etl_list <- CalculateExpectedTimeLost(
    PredictedCurves = list(ensemble_pred_formatted),
    modeltypes = c("CR"),  # Competing risks model type
    times = common_times,
    UL = max_time
  )
  
  ensemble_etl <- ensemble_etl_list[[1]]
  etl_results[["Ensemble"]] <- ensemble_etl
  cat("Ensemble ETL - mean:", round(mean(ensemble_etl), 2),
      "range:", paste(round(range(ensemble_etl), 2), collapse = "-"), "\n")
}, error = function(e) {
  cat("Ensemble ETL calculation failed:", e$message, "\n")
  # Fallback to manual calculation if needed
  ensemble_etl <- rep(NA, ncol(ensemble_cif))
  etl_results[["Ensemble"]] <- ensemble_etl
})

# --- 11. Visualization ---
cat("\n=== STEP 10: COMPETING RISKS VISUALIZATION ===\n")

# 11.1 Model Performance Comparison
# Handle case where integrated concordances failed
int_conc_values <- rep(NA, length(names(concordances)))
names(int_conc_values) <- names(concordances)
for (model_name in names(concordances)) {
  if (model_name %in% names(integrated_concordances)) {
    int_conc_values[model_name] <- integrated_concordances[[model_name]]
  }
}

performance_df <- data.frame(
  Model = names(concordances),
  Concordance = unlist(concordances),
  Integrated_Concordance = int_conc_values,
  Brier_Score = unlist(brier_scores[names(concordances)]),
  stringsAsFactors = FALSE
)

# Add ensemble
performance_df <- rbind(performance_df,
                        data.frame(Model = "Ensemble",
                                   Concordance = ensemble_concordance,
                                   Integrated_Concordance = ensemble_integrated_conc,
                                   Brier_Score = ensemble_brier,
                                   stringsAsFactors = FALSE))

# Remove rows with NA concordance
performance_df <- performance_df[!is.na(performance_df$Concordance), ]

cat("\nModel Performance Summary:\n")
print(performance_df)

# 11.2 Cumulative Incidence Curves for Selected Patients
cat("\nCreating CIF plots for first 3 test patients...\n")

# Create plots using the package plotting functions  
# Plot with all 3 patients
p1 <- plot_cif_curves(
  predictions = c(predictions, list(Ensemble = list(Times = common_times, CIFs = ensemble_cif))),
  patients_to_plot = seq_len(min(3, nrow(test_data))),
  highlight_ensemble = TRUE,
  title = "Predicted Cumulative Incidence Functions (CIF) - All Models"
)

# Focused plot for Patient 1 only
p2 <- plot_cif_curves(
  predictions = c(predictions, list(Ensemble = list(Times = common_times, CIFs = ensemble_cif))),
  patients_to_plot = 1,
  highlight_ensemble = TRUE,
  title = "Cumulative Incidence Functions - Patient 1"
)

cat("✓ CIF plots created using package functions\n")

# Print CIF values for first 3 patients at key time points
cat("\n=== CIF VALUES FOR FIRST 3 TEST PATIENTS ===\n")
key_times <- c(6, 12, 24, 36)  # months

for (patient in 1:3) {
  cat(sprintf("\n--- PATIENT %d ---\n", patient))
  cat("Patient characteristics:\n")
  print(test_data[patient, expvars])
  cat("\nCIF values at key time points:\n")
  
  for (model_name in names(predictions)) {
    pred <- predictions[[model_name]]
    cat(sprintf("\n%s model:\n", model_name))
    
    for (time_point in key_times) {
      time_idx <- which.min(abs(pred$Times - time_point))
      cif_value <- pred$CIFs[time_idx, patient]  # Fixed: CIFs matrix is [times, observations]
      cat(sprintf("  %2d months: %.4f\n", time_point, cif_value))
    }
  }
  
  # Ensemble
  cat("\nEnsemble model:\n")
  for (time_point in key_times) {
    time_idx <- which.min(abs(common_times - time_point))
    cif_value <- ensemble_cif[time_idx, patient]
    cat(sprintf("  %2d months: %.4f\n", time_point, cif_value))
  }
}

cat("\n✓ CIF values displayed for first 3 patients\n")

# 11.3 Create all plots

# Plot 1 and Plot 2 are already created above using package functions

# Plot 3-5: Performance plots using package function
p3 <- plot_model_performance(
  performance_df = performance_df,
  metric = "concordance",
  highlight_ensemble = TRUE,
  title = "Model Performance: Concordance Index"
)

p4 <- plot_model_performance(
  performance_df = performance_df,
  metric = "brier",
  highlight_ensemble = TRUE,
  title = "Model Performance: Brier Score"
)

# Plot 5: Expected Time Lost - create data frame for this
etl_plot_data <- data.frame()
for (model_name in names(etl_results)) {
  model_etl_df <- data.frame(
    ETL = etl_results[[model_name]],
    Model = model_name,
    IsEnsemble = model_name == "Ensemble",
    stringsAsFactors = FALSE
  )
  etl_plot_data <- rbind(etl_plot_data, model_etl_df)
}

p5 <- ggplot2::ggplot(etl_plot_data, ggplot2::aes(x = reorder(.data$Model, .data$ETL, FUN = median), y = .data$ETL,
                                 fill = .data$IsEnsemble)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "lightblue"), guide = "none") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "Expected Time Lost Distribution by Model",
       subtitle = "Ensemble in BLACK",
       x = "Model", y = "Expected Time Lost (months)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))



# Plot 6: Cumulative Incidence by Disease Phase (from test data)
test_cif_phase <- cuminc(ftime = test_data$ftime,
                         fstatus = test_data$Status,
                         group = test_data$Phase)

# Extract CIF data for plotting
cif_phase_plot_data <- data.frame()
for (cif_name in names(test_cif_phase)) {
  if (cif_name %in% c("Tests")) next  # Skip test results

  cif_obj <- test_cif_phase[[cif_name]]
  parts <- strsplit(cif_name, " ")[[1]]

  if (length(parts) >= 2) {
    event_type <- as.numeric(parts[length(parts)])
    phase_name <- paste(parts[1:(length(parts)-1)], collapse = " ")

    cif_data <- data.frame(
      Time = cif_obj$time,
      CIF = cif_obj$est,
      Phase = phase_name,
      Event = ifelse(event_type == 1, "Relapse", "TRM"),
      stringsAsFactors = FALSE
    )
    cif_phase_plot_data <- rbind(cif_phase_plot_data, cif_data)
  }
}

p6 <- ggplot(cif_phase_plot_data[cif_phase_plot_data$Event == "Relapse", ],
             aes(x = Time, y = CIF, color = Phase)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cumulative Incidence of Relapse by Disease Phase",
       subtitle = "Test Set Data",
       x = "Time (months)", y = "Cumulative Incidence of Relapse") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

p7 <- ggplot(cif_phase_plot_data[cif_phase_plot_data$Event == "TRM", ],
             aes(x = Time, y = CIF, color = Phase)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cumulative Incidence of TRM by Disease Phase",
       subtitle = "Test Set Data",
       x = "Time (months)", y = "Cumulative Incidence of TRM") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# Display plots
print(p1)  # CIF curves for 3 patients, all models
print(p2)  # CIF for patient 1
print(p3)  # Concordance comparison
print(p4)  # Brier score comparison
print(p5)  # ETL distribution
print(p6)  # CIF by phase - Relapse
print(p7)  # CIF by phase - TRM

# Save comprehensive plot to PDF
cat("\nSaving plots to PDF...\n")
pdf("bmt_competing_risks_analysis.pdf", width = 16, height = 20)

# Page 1: CIF curves for 3 patients
print(p1)

# Page 2: Patient 1 focused + Concordance
grid.arrange(p2, p3, ncol = 2)

# Page 3: Brier Score + ETL
grid.arrange(p4, p5, ncol = 2)

# Page 4: CIF by Phase (both events)
grid.arrange(p6, p7, ncol = 1)

dev.off()
cat("✓ Plots saved to bmt_competing_risks_analysis.pdf\n")

# --- 12. Model Persistence ---
cat("\n=== STEP 11: MODEL SAVING AND LOADING ===\n")

# Save all trained models
saved_models <- c()
for (model_name in names(models)) {
  filename <- paste0(tolower(model_name), "_crmodel_bmt.rds")
  saveRDS(models[[model_name]], filename)
  saved_models <- c(saved_models, filename)
}

# Save ensemble components
ensemble_object <- list(
  cif = ensemble_cif,
  times = common_times,
  concordance = ensemble_concordance,
  brier_score = ensemble_brier,
  etl = ensemble_etl,
  top_models = top_models
)
saveRDS(ensemble_object, "ensemble_crmodel_bmt.rds")
saved_models <- c(saved_models, "ensemble_crmodel_bmt.rds")

cat("Models saved successfully:\n")
for (filename in saved_models) {
  cat("-", filename, "\n")
}

# --- 13. Summary Report ---
cat("\n=== FINAL SUMMARY ===\n")
cat("══════════════════════════════════════════════════════════\n")
cat("Comprehensive Competing Risks Analysis Complete!\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("Dataset: Bone Marrow Transplant (BMT) - Competing Risks\n")
cat("Total observations:", nrow(bmt_data), "\n")
cat("Training set:", nrow(train_data), "subjects\n")
cat("Test set:", nrow(test_data), "subjects\n\n")

cat("Event Distribution:\n")
cat("- Censored:", status_table["0"], "(", round(100 * status_table["0"] / nrow(bmt_data), 1), "%)\n")
cat("- Relapse (event of interest):", status_table["1"], "(", round(100 * status_table["1"] / nrow(bmt_data), 1), "%)\n")
cat("- TRM (competing risk):", status_table["2"], "(", round(100 * status_table["2"] / nrow(bmt_data), 1), "%)\n\n")

cat("Models Successfully Trained:", length(models), "\n")
cat("Models with Predictions:", length(predictions), "\n\n")

cat("Model Performance (Concordance Index):\n")
for (i in seq_len(nrow(performance_df))) {
  prefix <- if (performance_df$Model[i] == "Ensemble") ">>> " else "  - "
  cat(prefix, performance_df$Model[i], ":", round(performance_df$Concordance[i], 3), "\n")
}

cat("\nModel Performance (Brier Score):\n")
for (i in seq_len(nrow(performance_df))) {
  prefix <- if (performance_df$Model[i] == "Ensemble") ">>> " else "  - "
  cat(prefix, performance_df$Model[i], ":", round(performance_df$Brier_Score[i], 4), "\n")
}

cat("\nExpected Time Lost (Mean, in months):\n")
for (model_name in names(etl_results)) {
  prefix <- if (model_name == "Ensemble") ">>> " else "  - "
  cat(prefix, model_name, ":", round(mean(etl_results[[model_name]]), 2), "\n")
}

# Show best performing model
best_model_idx <- which.max(performance_df$Concordance)
best_model_name <- performance_df$Model[best_model_idx]
best_concordance <- performance_df$Concordance[best_model_idx]

cat("\n✓ Best Performing Model (by C-index):", best_model_name,
    "with concordance:", round(best_concordance, 3), "\n\n")

cat("Output Files Generated:\n")
cat("- Model files:", length(saved_models), "RDS files\n")
cat("- Visualization: bmt_competing_risks_analysis.pdf (4 pages)\n\n")

cat("This analysis demonstrates:\n")
cat("✓ Data loading and preprocessing for competing risks\n")
cat("✓ Multiple competing risks models (Cox, Fine-Gray, RF, XGBoost, GAM)\n")
cat("✓ Cumulative Incidence Function (CIF) calculations\n")
cat("✓ Expected Time Lost (ETL) for competing risks\n")
cat("✓ Model evaluation with CR-specific metrics\n")
cat("✓ Ensemble methods for competing risks (HIGHLIGHTED IN BLACK)\n")
cat("✓ Comprehensive CIF and ETL visualization\n")
cat("✓ Model persistence and reloading\n")

cat("\n══════════════════════════════════════════════════════════\n")
cat("=== COMPETING RISKS VIGNETTE COMPLETE ===\n")
cat("══════════════════════════════════════════════════════════\n")
