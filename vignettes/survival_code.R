## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----load-libraries-----------------------------------------------------------
library(ml4time2event)
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

devtools::load_all()


## ----load-data----------------------------------------------------------------
lung_data <- ml4time2event::get_lung_survival_data()

glimpse(lung_data)

summary(lung_data[c("time", "status", "age", "ph.karno", "wt.loss")])


## ----explore-data-------------------------------------------------------------
status_counts <- dplyr::count(lung_data, status)
status_counts

event_rate <- mean(lung_data$status == 1)
time_range <- range(lung_data$time)

data.frame(
  metric = c("Event rate", "Minimum follow-up (days)", "Maximum follow-up (days)"),
  value = c(round(event_rate, 3), round(time_range[1], 1), round(time_range[2], 1))
)


## ----variable-profile---------------------------------------------------------
candidate_expvars <- c(
  "age", "sex", "ph.ecog", "ph.karno", "pat.karno",
  "meal.cal", "wt.loss", "age_group", "performance_good"
)
var_profile <- VariableProfile(lung_data, candidate_expvars)
var_profile$Summary


## ----split-data---------------------------------------------------------------
set.seed(123)
split_data <- t2edata_split(lung_data, prop = 0.7)
train_data <- split_data[["Train"]]
test_data <- split_data[["Test"]]

data.frame(
  dataset = c("Training", "Test"),
  observations = c(nrow(train_data), nrow(test_data))
)


## ----train-models-------------------------------------------------------------
timevar <- "time"
eventvar <- "status"
expvars <- c("age", "sex", "ph.ecog", "ph.karno", "meal.cal", "wt.loss")

data.frame(
  setting = c("Time variable", "Event variable", "Predictors"),
  value = c(timevar, eventvar, paste(expvars, collapse = ", "))
)

models <- RunSurvModels(
  datatrain = train_data,
  ExpVars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  models = c("glmnet", "coxph", "rulefit", "xgboost", "gam", "gbm", "ExpSurvReg", "bart", "deepsurv"),
  ntreeRF = 300,
  nvars = 20,
  ntimes = 50  # New parameter: number of time points for prediction grid
)

# Individual model components now include:
# - For Cox model: ntimes, alpha (elastic net), nfolds (cross-validation)
# - For XGBoost: ntimes (prediction grid), verbose (progress reporting)
# - Consistent time grid specification across all models

models


## ----generate-predictions-----------------------------------------------------
complete_test_data <- test_data[complete.cases(test_data[, expvars]), ]

# Note: With the updated API, time grids are now automatically managed by individual
# models. The ntimes parameter controls prediction grid resolution (default: 50).
# Individual model predictions include a 'Times' attribute for easy access.

ensemble_predictions <- PredictSurvModels(
  models = models,
  newdata = complete_test_data,
  new_times = seq(1, max(complete_test_data$time), length.out = 50),
  ensemble_method = "average"
)

individual_predictions <- ensemble_predictions$ModelPredictions
ensemble_probs <- ensemble_predictions$NewProbs
common_times <- ensemble_predictions$NewTimes
models_used <- ensemble_predictions$models_used

# Debug: Check ExpSurvReg predictions
cat("Names in individual_predictions:", paste(names(individual_predictions), collapse = ", "), "\n")
if ("survregexp_Model" %in% names(individual_predictions)) {
  cat("ExpSurvReg (survregexp_Model) predictions summary:\n")
  print(summary(individual_predictions$survregexp_Model))
  cat("Any NA in ExpSurvReg predictions:", any(is.na(individual_predictions$survregexp_Model)), "\n")
  cat("First few rows of ExpSurvReg predictions:\n")
  print(head(individual_predictions$survregexp_Model))
} else {
  cat("survregexp_Model not found in individual_predictions\n")
}

data.frame(
  metric = c(
    "Test subjects (complete cases)",
    "Prediction time points",
    "Models in ensemble"
  ),
  value = c(
    nrow(complete_test_data),
    length(common_times),
    paste(models_used, collapse = ", ")
  )
)


## ----prepare-predictions------------------------------------------------------
pretty_model_name <- function(x) {
  mapping <- c(
    RF_Model = "Random Forest",
    RF_Model2 = "Random Forest (Top Vars)",
    glmnet_Model = "GLMNet",
    CPH_Model = "Cox PH",
    bart_Model = "BART",
    deepsurv_Model = "DeepSurv",
    gam_Model = "GAM",
    gbm_Model = "GBM",
    survregexp_Model = "ExpSurvReg",
    survregweib_Model = "WeibSurvReg",
    xgboost_Model = "XGBoost",
    RuleFit_Model = "RuleFit"
  )
  if (x %in% names(mapping)) mapping[[x]] else tools::toTitleCase(gsub("_", " ", gsub("_Model$", "", x)))
}

model_prediction_objects <- lapply(individual_predictions, function(mat) {
  list(Times = common_times, Probs = mat)
})

names(model_prediction_objects) <- vapply(names(individual_predictions), pretty_model_name, character(1))

all_predictions <- c(
  model_prediction_objects,
  list(Ensemble = list(Times = common_times, Probs = ensemble_probs))
)

names(all_predictions)


## ----evaluate-models----------------------------------------------------------
actual_times <- complete_test_data[[timevar]]
actual_events <- complete_test_data[[eventvar]]
actual_events_binary <- as.integer(actual_events == 1)

median_time <- median(actual_times[actual_events_binary == 1], na.rm = TRUE)
max_time <- as.numeric(quantile(actual_times, 0.9, na.rm = TRUE))

etl_results <- list()

performance_summary <- do.call(rbind, lapply(names(all_predictions), function(model_name) {
  pred_obj <- all_predictions[[model_name]]

  concordance <- tryCatch({
    # Use predictions at median time for concordance
    time_idx <- which.min(abs(pred_obj$Times - median_time))
    pred_at_median <- pred_obj$Probs[time_idx, ]
    temp_data <- data.frame(time = actual_times, status = actual_events, pred = -pred_at_median)
    conc_obj <- survival::concordance(Surv(time, status) ~ pred, data = temp_data)
    conc_obj$concordance
  }, error = function(e) NA_real_)

  # BrierScore now uses 'eval_times' parameter (standardized across survival and CR)
  brier <- tryCatch({
    # Probs matrices are oriented as time x observation throughout the workflow
    brier_obj <- BrierScore(
      predsurv = pred_obj$Probs,
      pred_times = pred_obj$Times,
      obstimes = actual_times,
      obsevents = actual_events_binary,
      eval_times = median_time  # New standardized parameter name
    )
    if (!is.null(brier_obj$AppErr)) brier_obj$AppErr$model[1] else NA_real_
  }, error = function(e) NA_real_)

  etl_values <- tryCatch({
    etl_list <- CalculateExpectedTimeLost(
      PredictedCurves = list(list(NewProbs = pred_obj$Probs)),
      modeltypes = c("SURV"),
      times = pred_obj$Times,
      UL = max_time
    )
    etl_list[[1]]
  }, error = function(e) rep(NA_real_, ncol(pred_obj$Probs)))

  etl_results[[model_name]] <<- etl_values

  data.frame(
    Model = model_name,
    Concordance = concordance,
    Brier_Score = brier,
    Mean_ETL = mean(etl_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

performance_summary

