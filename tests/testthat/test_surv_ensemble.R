library(testthat)
library(here)
# library(data.table) # Removed
library(survival)

# Assuming the functions are available in the environment
source(here("R/surv_ensemble.R"))
# Need the individual model wrappers if they are called internally, or mock them.
# source(here("R/models/surv_cox.R")) # Example dependency
# source(here("R/models/surv_random_forest.R")) # Example dependency

context("Testing surv_ensemble functions")

# --- Test Data Setup ---
# Reusing the setup from previous Surv tests
set.seed(789)
n_obs_surv <- 50
surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)), # 0=censored, 1=event
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)),
  x3 = rnorm(n_obs_surv, mean = 2),
  stringsAsFactors = FALSE
)
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))

# --- Mock Model and Predict Functions ---
# Mock functions to simulate survival model fitting and prediction

# Mock Survival Model 1 (e.g., Cox)
Mock_SurvModel1 <- function(formula, data, ...) {
  # Returns a simple list structure mimicking a model object
  list(model_type = "MockSurv1", formula = formula, coef = rnorm(3))
}
Predict_Mock_SurvModel1 <- function(model, data, times, ...) {
  # Returns mock survival predictions (matrix subjects x times)
  n_test <- nrow(data)
  n_times <- length(times)
  # Simulate decreasing survival probabilities
  preds <- matrix(runif(n_test * n_times, 0.5, 0.9), nrow = n_test, ncol = n_times)
  preds <- t(apply(preds, 1, cummin)) # Ensure non-increasing
  return(preds)
}

# Mock Survival Model 2 (e.g., RF)
Mock_SurvModel2 <- function(formula, data, ...) {
  # Returns a different structure
  list(model_info = list(type = "MockSurv2"), formula = formula, importance = rnorm(3))
}
Predict_Mock_SurvModel2 <- function(model, data, times, ...) {
   # Returns mock survival predictions
  n_test <- nrow(data)
  n_times <- length(times)
  preds <- matrix(runif(n_test * n_times, 0.4, 0.8), nrow = n_test, ncol = n_times)
  preds <- t(apply(preds, 1, cummin)) # Ensure non-increasing
  return(preds)
}

# Define the list of models for the ensemble wrappers
# The names MUST match the function names (without "SurvModel_")
surv_model_list <- list(
  Mock_SurvModel1 = Mock_SurvModel1,
  Mock_SurvModel2 = Mock_SurvModel2
)
# Prediction functions need similar naming convention
surv_predict_list <- list(
  Predict_Mock_SurvModel1 = Predict_Mock_SurvModel1,
  Predict_Mock_SurvModel2 = Predict_Mock_SurvModel2
)


# --- Tests for RunSurvModels ---

test_that("RunSurvModels runs specified models and returns a list", {
  # Assuming RunSurvModels takes the list as an argument for testability
  fitted_models <- RunSurvModels(formula = surv_formula, data = train_data_surv,
                               models = c("Mock_SurvModel1", "Mock_SurvModel2"), # Names match keys in list
                               model_list = surv_model_list) # Pass the mock list

  expect_type(fitted_models, "list")
  # Expecting one model object per model type
  expect_length(fitted_models, 2)
  expect_named(fitted_models, c("Mock_SurvModel1", "Mock_SurvModel2"), ignore.order = TRUE)

  # Check individual model structures (based on mocks)
  expect_equal(fitted_models$Mock_SurvModel1$model_type, "MockSurv1")
  expect_equal(fitted_models$Mock_SurvModel2$model_info$type, "MockSurv2")
})

test_that("RunSurvModels handles single model", {
   fitted_models <- RunSurvModels(formula = surv_formula, data = train_data_surv,
                               models = "Mock_SurvModel1",
                               model_list = surv_model_list)
   expect_length(fitted_models, 1)
   expect_named(fitted_models, "Mock_SurvModel1")
})

test_that("RunSurvModels requires formula, data, models, model_list", {
  expect_error(RunSurvModels(data = train_data_surv, models="M1", model_list=surv_model_list), "formula")
  expect_error(RunSurvModels(formula = surv_formula, models="M1", model_list=surv_model_list), "data")
  expect_error(RunSurvModels(formula = surv_formula, data=train_data_surv, model_list=surv_model_list), "models")
  # expect_error(RunSurvModels(formula = surv_formula, data=train_data_surv, models="M1"), "model_list") # If required
})


# --- Tests for PredictSurvModels ---

test_that("PredictSurvModels generates predictions for specified models", {
  # Fit models first
  fitted_models <- RunSurvModels(formula = surv_formula, data = train_data_surv,
                               models = c("Mock_SurvModel1", "Mock_SurvModel2"),
                               model_list = surv_model_list)

  # Get predictions using the mock prediction functions
  all_predictions <- PredictSurvModels(models = fitted_models, data = test_data_surv, times = time_points_surv,
                                     predict_list = surv_predict_list) # Pass the mock list

  expect_type(all_predictions, "list")
  # Expecting predictions for each model type
  expect_length(all_predictions, 2) # Mock1, Mock2
  expect_named(all_predictions, c("Mock_SurvModel1", "Mock_SurvModel2"), ignore.order = TRUE)

  # Check structure for Mock_SurvModel1 predictions
  preds_m1 <- all_predictions$Mock_SurvModel1
  expect_true(is.matrix(preds_m1))
  expect_equal(nrow(preds_m1), nrow(test_data_surv))
  expect_equal(ncol(preds_m1), length(time_points_surv))
  expect_true(all(preds_m1 >= 0 & preds_m1 <= 1))
  if (length(time_points_surv) > 1) {
      expect_true(all(apply(preds_m1, 1, function(r) all(diff(r) <= 1e-9))))
  }


  # Check structure for Mock_SurvModel2 predictions
  preds_m2 <- all_predictions$Mock_SurvModel2
  expect_true(is.matrix(preds_m2))
  expect_equal(nrow(preds_m2), nrow(test_data_surv))
  expect_equal(ncol(preds_m2), length(time_points_surv))
  expect_true(all(preds_m2 >= 0 & preds_m2 <= 1))
   if (length(time_points_surv) > 1) {
      expect_true(all(apply(preds_m2, 1, function(r) all(diff(r) <= 1e-9))))
  }

})

test_that("PredictSurvModels requires models, data, times", {
   fitted_models <- RunSurvModels(formula = surv_formula, data = train_data_surv,
                               models = "Mock_SurvModel1", model_list = surv_model_list)
   expect_error(PredictSurvModels(data = test_data_surv, times = time_points_surv, predict_list = surv_predict_list), "models")
   expect_error(PredictSurvModels(models = fitted_models, times = time_points_surv, predict_list = surv_predict_list), "data")
   expect_error(PredictSurvModels(models = fitted_models, data = test_data_surv, predict_list = surv_predict_list), "times")
   # expect_error(PredictSurvModels(models = fitted_models, data = test_data_surv, times = time_points_surv), "predict_list") # If required
})
