knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(ml4time2event)
library(survival)
library(dplyr)
library(ggplot2)

if (!requireNamespace("cmprsk", quietly = TRUE)) {
  stop("Install the 'cmprsk' package to run this vignette.")
}
library(cmprsk)

bmt_data <- ml4time2event::get_bmt_competing_risks_data()

glimpse(bmt_data)

summary(bmt_data[c("ftime", "status", "age")])

# Event counts: 0=censored, 1=relapse, 2=treatment-related mortality
table(bmt_data$status)

candidate_expvars <- c("sex", "d", "phase", "age", "source")
var_profile <- VariableProfile(bmt_data, candidate_expvars)
var_profile$Summary

set.seed(123)
split_data <- t2edata_split(bmt_data, prop = 0.7)
train_data <- split_data$Train
test_data <- split_data$Test

timevar <- "ftime"
eventvar <- "status"
expvars <- c("sex", "d", "phase", "age", "source")
primary_event_code <- 1L

models <- RunCRModels(
  datatrain = train_data,
  ExpVars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  models = c("cox", "FG", "xgboost", "gam", "bart", "rulefit", "survreg", "ttah")
)

models

predictions <- PredictCRModels(
  models = models,
  newdata = test_data,
  new_times = seq(0, max(test_data$ftime), length.out = 50),
  ensemble_method = "average"
)

leaderboard <- EvaluateCRModels(
  models = models,
  data = test_data,
  timevar = timevar,
  eventvar = eventvar,
  eval_times = seq(0, max(test_data[[timevar]]), length.out = 50),
  ensemble_method = "average",
  cause = primary_event_code
)

leaderboard$display_name <- format_model_name(leaderboard$model, model_type = "competing_risks")
print(leaderboard[, c("display_name", "integrated_c", "integrated_brier", "rank")])

plot_cif_curves(
  predictions = predictions,
  patients_to_plot = 1:3,
  highlight_ensemble = TRUE
)

# Save ensemble
ensemble_file <- tempfile("cr_ensemble_model_", fileext = ".rds")
SaveEnsemble(models, file = ensemble_file)

# Load ensemble
loaded_ensemble <- LoadEnsemble(file = ensemble_file)
print(loaded_ensemble)
