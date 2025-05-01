library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(glmnet)   # Required for glmnet/cv.glmnet

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_glmnet.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing surv_glmnet functions")

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
# glmnet requires a numeric matrix for x, Surv object for y
# The wrapper function SurvModel_glmnet should handle formula conversion
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_glmnet ---

test_that("SurvModel_glmnet runs and returns a glmnet/cv.glmnet object", {
  skip_if_not_installed("glmnet")

  # Assuming the wrapper uses cv.glmnet by default
  model_glmnet_surv <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, nlambda = 5)

  # Check output type (cv.glmnet or glmnet)
  expect_true(inherits(model_glmnet_surv, "cv.glmnet") || inherits(model_glmnet_surv, "glmnet"))

  # Check basic model properties if possible (depends on whether cv.glmnet or glmnet is returned)
  if (inherits(model_glmnet_surv, "cv.glmnet")) {
     expect_true(!is.null(model_glmnet_surv$glmnet.fit))
     expect_true(!is.null(model_glmnet_surv$lambda.min))
  } else {
     expect_true(!is.null(model_glmnet_surv$beta))
  }
})

test_that("SurvModel_glmnet handles additional parameters (e.g., alpha)", {
  skip_if_not_installed("glmnet")
  # Test with alpha = 0 (Ridge)
  model_glmnet_ridge <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, alpha = 0)
  expect_true(inherits(model_glmnet_ridge, "cv.glmnet") || inherits(model_glmnet_ridge, "glmnet"))
  # Check if alpha was passed (might not be stored directly in the object easily accessible)

   # Test with alpha = 1 (Lasso)
  model_glmnet_lasso <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, alpha = 1, nlambda = 5)
  expect_true(inherits(model_glmnet_lasso, "cv.glmnet") || inherits(model_glmnet_lasso, "glmnet"))
})

test_that("SurvModel_glmnet requires formula and data", {
  skip_if_not_installed("glmnet")
  expect_error(SurvModel_glmnet(formula = surv_formula, nlambda = 5), "argument \"data\" is missing")
  expect_error(SurvModel_glmnet(data = train_data_surv, nlambda = 5), "argument \"formula\" is missing")
})


# --- Tests for Predict_SurvModel_glmnet ---

test_that("Predict_SurvModel_glmnet returns predictions in correct format", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survival") # Needed for baseline hazard estimation

  model_glmnet_surv <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, nlambda = 5)
  # Prediction likely involves predict(..., type="response") and baseline hazard
  # Assuming Predict_SurvModel_glmnet handles this.
  predictions <- Predict_SurvModel_glmnet(model = model_glmnet_surv, data = test_data_surv, times = time_points_surv)

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

test_that("Predict_SurvModel_glmnet handles single time point", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survival")

  model_glmnet_surv <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, nlambda = 5)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_glmnet(model = model_glmnet_surv, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_glmnet requires model, data, and times", {
  skip_if_not_installed("glmnet")
  model_glmnet_surv <- SurvModel_glmnet(formula = surv_formula, data = train_data_surv, nlambda = 5)
  expect_error(Predict_SurvModel_glmnet(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_glmnet(model = model_glmnet_surv, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_glmnet(model = model_glmnet_surv, data = test_data_surv), "argument \"times\" is missing")
})
