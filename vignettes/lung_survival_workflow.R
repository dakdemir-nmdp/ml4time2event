## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----load-libraries-----------------------------------------------------------
library(ml4time2event)
library(dplyr)
library(survival)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
}


## ----load-data----------------------------------------------------------------
lung_data <- ml4time2event::get_lung_survival_data()

dplyr::glimpse(lung_data)

candidate_expvars <- c(
  "age", "sex", "ph.ecog", "ph.karno", "pat.karno",
  "meal.cal", "wt.loss", "age_group", "performance_good"
)

var_profile <- VariableProfile(lung_data, candidate_expvars)
var_profile[c("age", "sex", "ph.ecog")]


## ----split-data---------------------------------------------------------------
set.seed(2025)
split_data <- t2edata_split(lung_data, prop = 0.75, strata = "status")
train_data <- split_data$Train
test_data <- split_data$Test

dplyr::tibble(
  dataset = c("Training", "Test"),
  observations = c(nrow(train_data), nrow(test_data)),
  events = c(sum(train_data$status == 1), sum(test_data$status == 1))
)


## ----build-recipe-------------------------------------------------------------
timevar <- "time"
eventvar <- "status"
idvars <- character(0)

base_recipe <- t2emodel_data_recipe_init(
  timevar = timevar,
  eventvar = eventvar,
  expvar = candidate_expvars,
  idvars = idvars,
  traindata = train_data
)

preprocess_recipe <- minimal_data_recipe(
  model_recipe = base_recipe,
  pmiss = 0.25,
  pother = 0.05,
  dummy = FALSE
)

prepped_recipe <- prep_data_recipe(preprocess_recipe, training = train_data)
summary(prepped_recipe)

train_prepared <- bake_data_recipe(prepped_recipe, data = train_data)
test_prepared <- bake_data_recipe(prepped_recipe, data = test_data)

processed_expvars <- setdiff(colnames(train_prepared), c(timevar, eventvar))
processed_expvars


## ----train-model--------------------------------------------------------------
cox_model <- SurvModel_Cox(
  data = train_prepared,
  expvars = processed_expvars,
  timevar = timevar,
  eventvar = eventvar,
  varsel = "penalized",
  alpha = 0.5,
  nfolds = 5,
  ntimes = 75
)

cox_model$expvars


## ----predict------------------------------------------------------------------
prediction_horizon <- seq(
  from = 0,
  to = max(test_prepared[[timevar]], na.rm = TRUE),
  length.out = 50
)

cox_predictions <- Predict_SurvModel_Cox(
  modelout = cox_model,
  newdata = test_prepared,
  new_times = prediction_horizon
)

head(cox_predictions$Times)
cox_predictions$Probs[seq_len(5), seq_len(min(3, ncol(cox_predictions$Probs)))]


## ----evaluate-----------------------------------------------------------------
test_times <- test_prepared[[timevar]]
test_events <- as.integer(test_prepared[[eventvar]] == 1)

eval_times <- c(180, 365, 730)

brier_results <- BrierScore(
  predsurv = cox_predictions$Probs,
  pred_times = cox_predictions$Times,
  obstimes = test_times,
  obsevents = test_events,
  eval_times = eval_times
)

dplyr::tibble(
  time = brier_results$AppErr$time,
  brier_score = brier_results$AppErr$model
)

integrated_brier <- integratedBrier(
  predsurv = cox_predictions$Probs,
  pred_times = cox_predictions$Times,
  obstimes = test_times,
  obsevents = test_events,
  eval_times = eval_times
)
integrated_brier

cindex_curve <- timedepConcordance(
  predsurv = cox_predictions$Probs,
  pred_times = cox_predictions$Times,
  obstimes = test_times,
  obsevents = test_events
)

dplyr::tibble(
  time = cindex_curve$time,
  c_index = cindex_curve$AppCindex$matrix
)

integrated_c <- integratedC(
  predsurv = cox_predictions$Probs,
  pred_times = cox_predictions$Times,
  obstimes = test_times,
  obsevents = test_events
)
integrated_c


## ----inspect-individual-------------------------------------------------------
dplyr::tibble(
  time = cox_predictions$Times,
  survival_probability = cox_predictions$Probs[, 1]
)
