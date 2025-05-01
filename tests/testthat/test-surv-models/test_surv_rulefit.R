library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
# library(pre) # Required for the functions being tested - load if needed, or skip

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_rulefit.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing surv_rulefit functions")

# --- Test Data Setup ---
# Reusing the setup from test_surv_random_forest.R
set.seed(789)
n_obs_surv <- 50
surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)), # 0=censored, 1=event
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)), # pre handles factors
  x3 = rnorm(n_obs_surv, mean = 2),
  stringsAsFactors = FALSE
)
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_rulefit ---

test_that("SurvModel_rulefit runs and returns a pre object", {
  skip_if_not_installed("pre")

  # Assuming the wrapper fits a standard 'pre' model for survival
  model_rulefit_surv <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv)

  # Check output type
  expect_s3_class(model_rulefit_surv, "pre")

  # Check basic model properties if possible
  expect_true(!is.null(model_rulefit_surv$rules))
  expect_true(!is.null(model_rulefit_surv$glmnet.fit))
})

test_that("SurvModel_rulefit handles additional parameters (if applicable)", {
  skip_if_not_installed("pre")
  # Example: if the wrapper allowed passing 'maxdepth'
  # model_rulefit_params <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv, maxdepth = 2)
  # expect_s3_class(model_rulefit_params, "pre")
  # Check if parameter was used (might require inspecting internal structure)

  # For now, just check it runs
  model_rulefit_surv <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv, maxdepth = 2)
  expect_s3_class(model_rulefit_surv, "pre")
})

test_that("SurvModel_rulefit requires formula and data", {
  skip_if_not_installed("pre")
  expect_error(SurvModel_rulefit(formula = surv_formula), "argument \"data\" is missing")
  expect_error(SurvModel_rulefit(data = train_data_surv), "argument \"formula\" is missing")
})


# --- Tests for Predict_SurvModel_rulefit ---

test_that("Predict_SurvModel_rulefit returns predictions in correct format", {
  skip_if_not_installed("pre")
  skip_if_not_installed("survival") # Needed for baseline hazard/survival estimation

  model_rulefit_surv <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv)
  # Prediction likely involves predict(..., type="response") and baseline hazard/survival
  # Assuming Predict_SurvModel_rulefit handles this.
  predictions <- Predict_SurvModel_rulefit(model = model_rulefit_surv, data = test_data_surv, times = time_points_surv)

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

test_that("Predict_SurvModel_rulefit handles single time point", {
  skip_if_not_installed("pre")
  skip_if_not_installed("survival")

  model_rulefit_surv <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv)
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_rulefit(model = model_rulefit_surv, data = test_data_surv, times = single_time)

  expect_true(is.matrix(predictions))
  expect_equal(ncol(predictions), 1) # Should have 1 column
  expect_equal(nrow(predictions), nrow(test_data_surv))
})

test_that("Predict_SurvModel_rulefit requires model, data, and times", {
  skip_if_not_installed("pre")
  model_rulefit_surv <- SurvModel_rulefit(formula = surv_formula, data = train_data_surv, maxdepth = 2)
  expect_error(Predict_SurvModel_rulefit(data = test_data_surv, times = time_points_surv), "argument \"model\" is missing")
  expect_error(Predict_SurvModel_rulefit(model = model_rulefit_surv, times = time_points_surv), "argument \"data\" is missing")
  expect_error(Predict_SurvModel_rulefit(model = model_rulefit_surv, data = test_data_surv), "argument \"times\" is missing")
})
