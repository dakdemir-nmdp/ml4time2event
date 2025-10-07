library(testthat)
library(here)
# library(data.table) # Removed
library(survival)

# Assuming the functions are available in the environment
source(here("R/utils/cr_ensemble.R"))
# Need the individual model wrappers if they are called internally, or mock them.
# source(here("R/models/cr_cox.R")) # Example dependency
# source(here("R/models/cr_random_forest.R")) # Example dependency

context("Testing cr_ensemble functions")

# --- Test Data Setup ---
# Reusing the setup from previous CR tests
set.seed(123)
n_obs <- 50
cr_data <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)
cr_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices <- 1:40
test_indices <- 41:50
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]
time_points <- c(quantile(train_data$time[train_data$status != 0], 0.25),
                 median(train_data$time[train_data$status != 0]),
                 quantile(train_data$time[train_data$status != 0], 0.75))

# --- Mock Model and Predict Functions ---
# Mock functions to simulate model fitting and prediction

# Mock Model 1 (e.g., Cox)
Mock_CRModel1 <- function(formula, data, failcode = NULL, ...) {
  # Returns a simple list structure mimicking a model object
  list(model_type = "Mock1", failcode = failcode, formula = formula, coef = rnorm(3))
}
Predict_Mock_CRModel1 <- function(model, data, times, ...) {
  # Returns mock predictions (list of matrices, one per cause)
  n_test <- nrow(data)
  n_times <- length(times)
  # Predict only for the failcode the model was 'built' for
  pred_list <- list()
  if (model$failcode == 1) {
      pred_list[[1]] <- matrix(runif(n_test * n_times, 0, 0.4), nrow = n_test, ncol = n_times)
      pred_list[[2]] <- matrix(NA, nrow = n_test, ncol = n_times) # No prediction for other cause
  } else {
      pred_list[[1]] <- matrix(NA, nrow = n_test, ncol = n_times)
      pred_list[[2]] <- matrix(runif(n_test * n_times, 0, 0.4), nrow = n_test, ncol = n_times)
  }
  return(pred_list)
}

# Mock Model 2 (e.g., RF)
Mock_CRModel2 <- function(formula, data, failcode = NULL, ...) {
  # Returns a different structure
  list(model_info = list(type = "Mock2", failcode = failcode), formula = formula, importance = rnorm(3))
}
Predict_Mock_CRModel2 <- function(model, data, times, ...) {
   # Returns mock predictions
  n_test <- nrow(data)
  n_times <- length(times)
  pred_list <- list()
   if (model$model_info$failcode == 1) {
      pred_list[[1]] <- matrix(runif(n_test * n_times, 0.1, 0.5), nrow = n_test, ncol = n_times)
      pred_list[[2]] <- matrix(NA, nrow = n_test, ncol = n_times)
  } else {
      pred_list[[1]] <- matrix(NA, nrow = n_test, ncol = n_times)
      pred_list[[2]] <- matrix(runif(n_test * n_times, 0.1, 0.5), nrow = n_test, ncol = n_times)
  }
  return(pred_list)
}

# Define the list of models for the ensemble wrappers
# The names MUST match the function names (without "CRModel_")
cr_model_list <- list(
  Mock_CRModel1 = Mock_CRModel1,
  Mock_CRModel2 = Mock_CRModel2
)
# Prediction functions need similar naming convention
cr_predict_list <- list(
  Predict_Mock_CRModel1 = Predict_Mock_CRModel1,
  Predict_Mock_CRModel2 = Predict_Mock_CRModel2
)


# --- Tests for RunCRModels ---

test_that("RunCRModels runs specified models and returns a list", {
  # Need to pass the mock list to the function if it doesn't discover them globally
  # Assuming RunCRModels takes the list as an argument for testability
  fitted_models <- RunCRModels(formula = cr_formula, data = train_data,
                               models = c("Mock_CRModel1", "Mock_CRModel2"), # Names match keys in list
                               model_list = cr_model_list, # Pass the mock list
                               failcodes = c(1, 2)) # Specify failcodes

  expect_type(fitted_models, "list")
  # Expecting models for each type and each failcode
  expect_length(fitted_models, 4) # Mock1*2 failcodes + Mock2*2 failcodes
  expect_named(fitted_models, c("Mock_CRModel1_1", "Mock_CRModel1_2", "Mock_CRModel2_1", "Mock_CRModel2_2"), ignore.order = TRUE)

  # Check individual model structures (based on mocks)
  expect_equal(fitted_models$Mock_CRModel1_1$model_type, "Mock1")
  expect_equal(fitted_models$Mock_CRModel1_1$failcode, 1)
  expect_equal(fitted_models$Mock_CRModel2_2$model_info$type, "Mock2")
  expect_equal(fitted_models$Mock_CRModel2_2$model_info$failcode, 2)
})

test_that("RunCRModels handles single model and single failcode", {
   fitted_models <- RunCRModels(formula = cr_formula, data = train_data,
                               models = "Mock_CRModel1",
                               model_list = cr_model_list,
                               failcodes = 1)
   expect_length(fitted_models, 1)
   expect_named(fitted_models, "Mock_CRModel1_1")
   expect_equal(fitted_models$Mock_CRModel1_1$failcode, 1)
})

test_that("RunCRModels requires formula, data, models, model_list, failcodes", {
  expect_error(RunCRModels(data = train_data, models="M1", model_list=cr_model_list, failcodes=1), "formula")
  expect_error(RunCRModels(formula = cr_formula, models="M1", model_list=cr_model_list, failcodes=1), "data")
  expect_error(RunCRModels(formula = cr_formula, data=train_data, model_list=cr_model_list, failcodes=1), "models")
  # expect_error(RunCRModels(formula = cr_formula, data=train_data, models="M1", failcodes=1), "model_list") # If required
  expect_error(RunCRModels(formula = cr_formula, data=train_data, models="M1", model_list=cr_model_list), "failcodes")
})


# --- Tests for PredictCRModels ---

test_that("PredictCRModels generates predictions for specified models", {
  # Fit models first
  fitted_models <- RunCRModels(formula = cr_formula, data = train_data,
                               models = c("Mock_CRModel1", "Mock_CRModel2"),
                               model_list = cr_model_list,
                               failcodes = c(1, 2))

  # Get predictions using the mock prediction functions
  all_predictions <- PredictCRModels(models = fitted_models, data = test_data, times = time_points,
                                     predict_list = cr_predict_list) # Pass the mock list

  expect_type(all_predictions, "list")
  # Expecting predictions for each model type (combining failcodes)
  expect_length(all_predictions, 2) # Mock1, Mock2
  expect_named(all_predictions, c("Mock_CRModel1", "Mock_CRModel2"), ignore.order = TRUE)

  # Check structure for Mock_CRModel1 predictions
  preds_m1 <- all_predictions$Mock_CRModel1
  expect_type(preds_m1, "list")
  expect_length(preds_m1, 2) # Cause 1, Cause 2
  expect_true(is.matrix(preds_m1[[1]]))
  expect_true(is.matrix(preds_m1[[2]]))
  expect_equal(nrow(preds_m1[[1]]), nrow(test_data))
  expect_equal(ncol(preds_m1[[1]]), length(time_points))
  # Check that NAs from mock predict functions were handled/combined correctly
  expect_false(any(is.na(preds_m1[[1]]))) # Cause 1 preds should exist
  expect_false(any(is.na(preds_m1[[2]]))) # Cause 2 preds should exist

  # Check structure for Mock_CRModel2 predictions
  preds_m2 <- all_predictions$Mock_CRModel2
  expect_type(preds_m2, "list")
  expect_length(preds_m2, 2)
  expect_true(is.matrix(preds_m2[[1]]))
  expect_true(is.matrix(preds_m2[[2]]))
  expect_equal(nrow(preds_m2[[1]]), nrow(test_data))
  expect_equal(ncol(preds_m2[[1]]), length(time_points))
  expect_false(any(is.na(preds_m2[[1]])))
  expect_false(any(is.na(preds_m2[[2]])))

  # Check probabilities sum <= 1
  expect_true(all(preds_m1[[1]] + preds_m1[[2]] <= 1.0001, na.rm = TRUE))
  expect_true(all(preds_m2[[1]] + preds_m2[[2]] <= 1.0001, na.rm = TRUE))

})

test_that("PredictCRModels requires models, data, times", {
   fitted_models <- RunCRModels(formula = cr_formula, data = train_data,
                               models = "Mock_CRModel1", model_list = cr_model_list, failcodes = 1)
   expect_error(PredictCRModels(data = test_data, times = time_points, predict_list = cr_predict_list), "models")
   expect_error(PredictCRModels(models = fitted_models, times = time_points, predict_list = cr_predict_list), "data")
   expect_error(PredictCRModels(models = fitted_models, data = test_data, predict_list = cr_predict_list), "times")
   # expect_error(PredictCRModels(models = fitted_models, data = test_data, times = time_points), "predict_list") # If required
})
