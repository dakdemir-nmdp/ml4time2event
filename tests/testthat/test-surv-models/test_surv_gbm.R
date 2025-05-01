library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(gbm)      # Required for gbm

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_gbm.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing surv_gbm functions")

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
# gbm with distribution='coxph' requires status to be 1 for event, 0 for censored
# The Surv object handles this representation.
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_gbm ---

test_that("SurvModel_gbm runs and returns a gbm object", {
  skip_if_not_installed("gbm")

  # gbm requires n.trees parameter
  model_gbm_surv <- SurvModel_gbm(formula = surv_formula, data = train_data_surv, n.trees = 10) # Provide n.trees

  # Check output type
  expect_s3_class(model_gbm_surv, "gbm")

  # Check basic model properties
  expect_equal(model_gbm_surv$distribution$name, "coxph") # Assuming it uses coxph distribution
  expect_equal(model_gbm_surv$n.trees, 10)
  expect_equal(model_gbm_surv$nTrain, nrow(train_data_surv))
})

test_that("SurvModel_gbm handles additional parameters", {
  skip_if_not_installed("gbm")
  # Test with different interaction.depth
  model_gbm_params_surv <- SurvModel_gbm(formula = surv_formula, data = train_data_surv, n.trees = 10, interaction.depth = 2)
  expect_s3_class(model_gbm_params_surv, "gbm")
  expect_equal(model_gbm_params_surv$interaction.depth, 2)
})

test_that("SurvModel_gbm requires formula and data", {
  skip_if_not_installed("gbm")
  expect_error(SurvModel_gbm(formula = surv_formula, n.trees = 10), "argument \"data\" is missing")
  expect_error(SurvModel_gbm(data = train_data_surv, n.trees = 10), "argument \"formula\" is missing")
  # Also check if n.trees is required by the wrapper
  # expect_error(SurvModel_gbm(formula = surv_formula, data = train_data_surv), "n.trees")
})


# --- Tests for Predict_SurvModel_gbm ---

test_that("Predict_SurvModel_gbm returns predictions in correct format", {
  skip_if_not_installed("gbm")

  model_gbm_surv <- SurvModel_gbm(formula = surv_formula, data = train_data_surv, n.trees = 10)
  # Prediction for gbm survival models often involves predict(..., type="response") giving relative risk,
  # then combining with baseline hazard (e.g., from survfit on the gbm object).
  # Assuming Predict_SurvModel_gbm handles this internally.
  predictions <- Predict_SurvModel_gbm(model = model_gbm_surv, data = test_data_surv, times = time_points_surv)

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

test_that("Predict_SurvModel_gbm handles single time point", {
  skip_if_not_installed("gbm")

  model_gbm_surv <- SurvModel_gbm(formula = surv_formula, data = train_data_surv, n.trees = 10)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_gbm(model = model_gbm_surv, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_gbm requires model, data, and times", {
  skip_if_not_installed("gbm")
  model_gbm_surv <- SurvModel_gbm(formula = surv_formula, data = train_data_surv, n.trees = 10)
  expect_error(Predict_SurvModel_gbm(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_gbm(model = model_gbm_surv, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_gbm(model = model_gbm_surv, data = test_data_surv), "argument \"times\" is missing")
})
