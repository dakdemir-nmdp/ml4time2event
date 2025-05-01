library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # Required for Surv and survreg

# Assuming the functions are available in the environment
source(here("R/models/surv_reg.R"))

context("Testing surv_reg functions")

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
# survreg uses Surv object directly
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_SurvReg ---

test_that("SurvModel_SurvReg runs and returns a survreg object", {
  skip_if_not_installed("survival")

  # Assuming default distribution is weibull if not specified
  model_survreg <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv)

  # Check output type
  expect_s3_class(model_survreg, "survreg")

  # Check basic model properties
  expect_true(!is.null(model_survreg$coefficients))
  expect_equal(model_survreg$dist, "weibull") # Check default distribution
})

test_that("SurvModel_SurvReg handles different distributions", {
  skip_if_not_installed("survival")
  # Test with exponential distribution
  model_survreg_exp <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv, dist = "exponential")
  expect_s3_class(model_survreg_exp, "survreg")
  expect_equal(model_survreg_exp$dist, "exponential")

  # Test with lognormal
  model_survreg_lnorm <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv, dist = "lognormal")
  expect_s3_class(model_survreg_lnorm, "survreg")
  expect_equal(model_survreg_lnorm$dist, "lognormal")
})

test_that("SurvModel_SurvReg requires formula and data", {
  skip_if_not_installed("survival")
  expect_error(SurvModel_SurvReg(formula = surv_formula), "argument \"data\" is missing")
  expect_error(SurvModel_SurvReg(data = train_data_surv), "argument \"formula\" is missing")
})


# --- Tests for Predict_SurvModel_SurvReg ---

test_that("Predict_SurvModel_SurvReg returns predictions in correct format", {
  skip_if_not_installed("survival")

  model_survreg <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv)
  # Prediction involves predict(..., type="quantile", p=...) and transforming to survival probability
  # Assuming Predict_SurvModel_SurvReg handles this.
  predictions <- Predict_SurvModel_SurvReg(model = model_survreg, data = test_data_surv, times = time_points_surv)

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

test_that("Predict_SurvModel_SurvReg handles single time point", {
  skip_if_not_installed("survival")

  model_survreg <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_SurvReg(model = model_survreg, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_SurvReg requires model, data, and times", {
  skip_if_not_installed("survival")
  model_survreg <- SurvModel_SurvReg(formula = surv_formula, data = train_data_surv)
  expect_error(Predict_SurvModel_SurvReg(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_SurvReg(model = model_survreg, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_SurvReg(model = model_survreg, data = test_data_surv), "argument \"times\" is missing")
})
