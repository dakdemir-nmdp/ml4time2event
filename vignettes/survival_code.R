knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(ml4time2event)
library(survival)
library(dplyr)
library(ggplot2)

lung_data <- ml4time2event::get_lung_survival_data()

glimpse(lung_data)

summary(lung_data[c("time", "status", "age", "ph.karno", "wt.loss")])

# Event distribution
table(lung_data$status)

# Summary statistics
summary(lung_data[c("time", "status")])

candidate_expvars <- c(
  "age", "sex", "ph.ecog", "ph.karno", "pat.karno",
  "meal.cal", "wt.loss", "age_group", "performance_good"
)
var_profile <- VariableProfile(lung_data, candidate_expvars)
var_profile$Summary

split_data <- t2edata_split(lung_data, prop = 0.7)
train_data <- split_data$Train
test_data <- split_data$Test

timevar <- "time"
eventvar <- "status"
expvars <- c("age", "sex", "ph.ecog", "ph.karno", "meal.cal", "wt.loss")

models <- RunSurvModels(
  datatrain = train_data,
  ExpVars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  models = c("glmnet", "coxph", "rulefit", "xgboost", "gam", "gbm", "ExpSurvReg", "bart", "shallownn", "ttah"),
  ntimes = 50
)

models

complete_test_data <- test_data[complete.cases(test_data[, expvars]), ]
prediction_times <- sort(unique(complete_test_data[[timevar]]))

# Generate predictions
predictions <- PredictSurvModels(
  models = models,
  newdata = complete_test_data,
  new_times = prediction_times,
  ensemble_method = "average"
)

leaderboard <- EvaluateSurvModels(
  models = models,
  data = complete_test_data,
  timevar = timevar,
  eventvar = eventvar,
  eval_times = prediction_times,
  prediction_times = prediction_times,
  ensemble_method = "average"
)

leaderboard$display_name <- format_model_name(leaderboard$model, model_type = "survival")
print(leaderboard[, c("display_name", "integrated_c", "integrated_brier", "rank")])

plot_survival_curves(
  predictions = predictions,
  patients_to_plot = 1:3,
  highlight_ensemble = TRUE
)

# Save ensemble
ensemble_file <- tempfile("lung_cancer_ensemble_", fileext = ".rds")
SaveEnsemble(models, file = ensemble_file)

# Load ensemble
loaded_ensemble <- LoadEnsemble(file = ensemble_file)
print(loaded_ensemble)
