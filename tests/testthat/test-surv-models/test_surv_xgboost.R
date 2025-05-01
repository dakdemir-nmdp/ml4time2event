library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(xgboost)  # Required for xgboost

# Assuming the functions are available in the environment
source(here("R/models/surv_xgboost.R"))

context("Testing surv_xgboost functions")

# --- Test Data Setup ---
# Reusing the setup from test_surv_random_forest.R
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
# xgboost requires numeric matrix, Surv object needs special handling (e.g., aft or cox objective)
# The wrapper function SurvModel_xgboost should handle formula conversion and data prep (e.g., xgb.DMatrix)
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_xgboost ---

test_that("SurvModel_xgboost runs and returns an xgb.Booster object", {
  skip_if_not_installed("xgboost")

  # xgboost requires nrounds
  model_xgb_surv <- SurvModel_xgboost(formula = surv_formula, data = train_data_surv, nrounds = 10) # Provide nrounds

  # Check output type
  expect_s3_class(model_xgb_surv, "xgb.Booster")

  # Check basic model properties if possible
  # expect_true(!is.null(model_xgb_surv$handle)) # Internal handle
})

test_that("SurvModel_xgboost handles additional parameters (e.g., eta)", {
  skip_if_not_installed("xgboost")
  # Test with different eta
  model_xgb_params_surv <- SurvModel_xgboost(formula = surv_formula, data = train_data_surv, nrounds = 10, eta = 0.1)
  expect_s3_class(model_xgb_params_surv, "xgb.Booster")
  # Checking passed parameters in xgboost object can be tricky, might need to check call if stored
})

test_that("SurvModel_xgboost requires formula and data", {
  skip_if_not_installed("xgboost")
  expect_error(SurvModel_xgboost(formula = surv_formula, nrounds = 10), "argument \"data\" is missing")
  expect_error(SurvModel_xgboost(data = train_data_surv, nrounds = 10), "argument \"formula\" is missing")
  # Check if nrounds is required by the wrapper
  # expect_error(SurvModel_xgboost(formula = surv_formula, data = train_data_surv), "nrounds")
})


# --- Tests for Predict_SurvModel_xgboost ---

test_that("Predict_SurvModel_xgboost returns predictions in correct format", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("survival") # Needed for baseline hazard/survival estimation

  model_xgb_surv <- SurvModel_xgboost(formula = surv_formula, data = train_data_surv, nrounds = 10)
  # Prediction likely involves predict(..., type="response") and baseline hazard/survival
  # Assuming Predict_SurvModel_xgboost handles this.
  predictions <- Predict_SurvModel_xgboost(model = model_xgb_surv, data = test_data_surv, times = time_points_surv)

  # Check output structure (should be a matrix of survival probabilities)
  expect_true(is.matrix(predictions))

  # Check dimensions
  # Rows = number of test observations, Cols = number of time points
  expect_equal(nrow(predictions), nrow(test_data_surv))
  expect_equal(ncol(predictions), length(time_points_surv))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions >= 0 & predictions <= 1, na.rm = TRUE))

  # Check that survival probabilities are non-increasing over time for each subject
  if (length(time_points_surv) > 1) {
    all_non_increasing <- all(apply(predictions, 1, function(row) all(diff(row) <= 1e-9))) # Allow for small tolerance
    expect_true(all_non_increasing)
  }
})

test_that("Predict_SurvModel_xgboost handles single time point", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("survival")

  model_xgb_surv <- SurvModel_xgboost(formula = surv_formula, data = train_data_surv, nrounds = 10)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_xgboost(model = model_xgb_surv, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_xgboost requires model, data, and times", {
  skip_if_not_installed("xgboost")
  model_xgb_surv <- SurvModel_xgboost(formula = surv_formula, data = train_data_surv, nrounds = 10)
  expect_error(Predict_SurvModel_xgboost(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_xgboost(model = model_xgb_surv, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_xgboost(model = model_xgb_surv, data = test_data_surv), "argument \"times\" is missing")
})
