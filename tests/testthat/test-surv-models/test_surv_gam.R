library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(mgcv)     # Required for gam

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_gam.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing surv_gam functions")

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
# GAM formula might need smooth terms, e.g., s(x1)
# Assuming the wrapper handles a standard formula for now
surv_formula <- Surv(time, status) ~ s(x1, k=3) + x2 + s(x3, k=3) # Example GAM formula
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_GAM ---

test_that("SurvModel_GAM runs and returns a gam object", {
  skip_if_not_installed("mgcv")

  # Need to know how SurvModel_GAM handles the Surv object for GAM fitting
  # Common approaches: Cox using gam(..., family=cox.ph), or specific survival families like weibull, cox.ph etc.
  # Assuming it uses a compatible family like cox.ph or similar logic internally.
  model_gam_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv)

  # Check output type
  expect_s3_class(model_gam_surv, "gam")

  # Check basic model properties
  expect_true(!is.null(model_gam_surv$coefficients))
  # Check if family is appropriate (e.g., cox.ph if that's used)
  # expect_equal(model_gam_surv$family$family, "cox.ph") # Example check
})

test_that("SurvModel_GAM handles additional parameters (if applicable)", {
  skip_if_not_installed("mgcv")
  # Example: if the wrapper allowed passing 'method' argument
  # model_gam_params_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv, method = "REML")
  # expect_s3_class(model_gam_params_surv, "gam")
  # expect_equal(model_gam_params_surv$method, "REML")

  # For now, just check it runs
   model_gam_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv)
   expect_s3_class(model_gam_surv, "gam")
})

test_that("SurvModel_GAM requires formula and data", {
  skip_if_not_installed("mgcv")
  expect_error(SurvModel_GAM(formula = surv_formula), "argument \"data\" is missing")
  expect_error(SurvModel_GAM(data = train_data_surv), "argument \"formula\" is missing")
})


# --- Tests for Predict_SurvModel_GAM ---

test_that("Predict_SurvModel_GAM returns predictions in correct format", {
  skip_if_not_installed("mgcv")

  model_gam_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv)
  # Prediction for GAM survival models is complex.
  # It might involve predict.gam(type="lp") and then manual calculation of survival function,
  # or using specific helper functions/packages.
  # Assuming Predict_SurvModel_GAM handles this internally.
  predictions <- Predict_SurvModel_GAM(model = model_gam_surv, data = test_data_surv, times = time_points_surv)

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

test_that("Predict_SurvModel_GAM handles single time point", {
  skip_if_not_installed("mgcv")

  model_gam_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_GAM(model = model_gam_surv, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_GAM requires model, data, and times", {
  skip_if_not_installed("mgcv")
  model_gam_surv <- SurvModel_GAM(formula = surv_formula, data = train_data_surv)
  expect_error(Predict_SurvModel_GAM(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_GAM(model = model_gam_surv, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_GAM(model = model_gam_surv, data = test_data_surv), "argument \"times\" is missing")
})
