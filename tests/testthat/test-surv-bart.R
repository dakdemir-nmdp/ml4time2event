library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
# library(dbarts) # Required for the functions being tested - load if needed, or skip

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/general_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_bart.R")

context("Testing surv_bart functions")

# --- Test Data Setup ---
# Reusing the setup from test_surv_random_forest.R
set.seed(789)
n_obs_surv <- 50
surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)), # 0=censored, 1=event
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)), # BART might handle factors
  x3 = rnorm(n_obs_surv, mean = 2),
  stringsAsFactors = FALSE
)
surv_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_BART ---

test_that("SurvModel_BART runs and returns a model object", {
  skip_if_not_installed("BART")

  # Assuming the wrapper fits a standard BART survival model
  model_bart_surv <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status")

  # Check output type
  expect_type(model_bart_surv, "list")
  # Check for expected elements within the BART model object
  expect_true(!is.null(model_bart_surv$model)) # BART returns model object
  expect_true(!is.null(model_bart_surv$expvars))
})

test_that("SurvModel_BART handles additional parameters (if applicable)", {
  skip_if_not_installed("BART")
  # Example: if the wrapper allowed passing 'ntree'
  model_bart_params <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status", ntree = 10)
  expect_type(model_bart_params, "list")
  # Check if parameter was used (might require inspecting internal structure or call)
  # expect_equal(model_bart_params$call$ntree, 50) # Example check

  # For now, just check it runs
  model_bart_surv <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status")
  expect_type(model_bart_surv, "list")
})

test_that("SurvModel_BART requires data, expvars, timevar, and eventvar", {
  skip_if_not_installed("BART")
  expect_error(SurvModel_BART(), "argument \"data\" is missing")
  expect_error(SurvModel_BART(data = train_data_surv), "argument \"expvars\" is missing")
  expect_error(SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3")), "argument \"timevar\" is missing")
  expect_error(SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time"), "argument \"eventvar\" is missing")
})


# --- Tests for Predict_SurvModel_BART ---

test_that("Predict_SurvModel_BART returns predictions in correct format", {
  skip_if_not_installed("BART")

  model_bart_surv <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status")
  # Prediction for BART survival models might involve predict() on the fitted object
  # and potentially custom logic to get survival probabilities at specific times.
  # Assuming Predict_SurvModel_BART handles this.
  predictions <- Predict_SurvModel_BART(modelout = model_bart_surv, newdata = test_data_surv, new_times = time_points_surv)

  # Check output structure
  expect_type(predictions, "list")

  # Check dimensions
  # Rows = number of time points, Cols = number of test observations
  expect_equal(nrow(predictions$Probs), length(time_points_surv))
  expect_equal(ncol(predictions$Probs), nrow(test_data_surv))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions$Probs >= 0 & predictions$Probs <= 1, na.rm = TRUE))

  # Check that survival probabilities are non-increasing over time for each subject
  if (length(time_points_surv) > 1) {
    all_non_increasing <- all(apply(predictions$Probs, 2, function(col) all(diff(col) <= 1e-9))) # Allow for small tolerance
    expect_true(all_non_increasing)
  }
})

test_that("Predict_SurvModel_BART handles single time point", {
  skip_if_not_installed("BART")

  model_bart_surv <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status")
  single_time <- median(train_data_surv$time[train_data_surv$status == 1])
  predictions <- Predict_SurvModel_BART(modelout = model_bart_surv, newdata = test_data_surv, new_times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions$Probs))
  expect_equal(ncol(predictions$Probs), nrow(test_data_surv))
  expect_equal(nrow(predictions$Probs), 1) # Single time point
})

test_that("Predict_SurvModel_BART requires model and data", {
  skip_if_not_installed("BART")
  model_bart_surv <- SurvModel_BART(data = train_data_surv, expvars = c("x1", "x2", "x3"), timevar = "time", eventvar = "status")
  expect_error(Predict_SurvModel_BART(newdata = test_data_surv), "argument \"modelout\" is missing")
  expect_error(Predict_SurvModel_BART(modelout = model_bart_surv), "argument \"newdata\" is missing")
  # new_times is optional, so no error expected when omitted
})
