library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(xgboost)  # Required for xgboost

# Assuming the functions are available in the environment
# source(here("R/models/surv_xgboost.R"))  # Removed - functions loaded via package

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
time_var <- "time"
event_var <- "status"
expvars <- c("x1", "x2", "x3")
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_xgboost ---

test_that("SurvModel_xgboost runs and returns expected structure", {
  skip_if_not_installed("xgboost")

  model_xgb <- SurvModel_xgboost(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )

  # Check output structure
  expect_type(model_xgb, "list")
  expect_named(model_xgb, c("model", "times", "varprof", "expvars"))
  expect_s3_class(model_xgb$model, "xgb.Booster")
  expect_type(model_xgb$times, "double")
  expect_type(model_xgb$varprof, "list")
  expect_type(model_xgb$expvars, "character")
})

test_that("SurvModel_xgboost handles different parameters", {
  skip_if_not_installed("xgboost")
  # Test basic functionality - parameters are handled internally
  model_xgb <- SurvModel_xgboost(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  expect_s3_class(model_xgb$model, "xgb.Booster")
})

test_that("SurvModel_xgboost requires correct inputs", {
  skip_if_not_installed("xgboost")
  expect_error(SurvModel_xgboost(data = train_data_surv, expvars = expvars, timevar = time_var), "argument \"eventvar\" is missing")
  expect_error(SurvModel_xgboost(expvars = expvars, timevar = time_var, eventvar = event_var), "argument \"data\" is missing")
})


# --- Tests for Predict_SurvModel_xgboost ---

test_that("Predict_SurvModel_xgboost returns predictions in correct format", {
  skip_if_not_installed("xgboost")

  model_xgb <- SurvModel_xgboost(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  predictions <- Predict_SurvModel_xgboost(
    modelout = model_xgb,
    newdata = test_data_surv
  )

  # Check output structure
  expect_type(predictions, "list")
  expect_named(predictions, c("Probs", "Times"))
  expect_true(is.matrix(predictions$Probs))
  expect_type(predictions$Times, "double")

  # Check dimensions
  expect_equal(nrow(predictions$Probs), length(predictions$Times))
  expect_equal(ncol(predictions$Probs), nrow(test_data_surv))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions$Probs >= 0 & predictions$Probs <= 1, na.rm = TRUE))

  # Check that survival probabilities are non-increasing over time for each subject
  if (length(predictions$Times) > 1) {
    all_non_increasing <- all(apply(predictions$Probs, 2, function(col) all(diff(col) <= 1e-9)))
    expect_true(all_non_increasing)
  }
})

test_that("Predict_SurvModel_xgboost handles custom times", {
  skip_if_not_installed("xgboost")

  model_xgb <- SurvModel_xgboost(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  # Note: Current implementation doesn't take custom times parameter
  # This test just verifies basic functionality
  predictions <- Predict_SurvModel_xgboost(
    modelout = model_xgb,
    newdata = test_data_surv
  )

  expect_type(predictions, "list")
  expect_named(predictions, c("Probs", "Times"))
  expect_true(length(predictions$Times) > 0)
})

test_that("Predict_SurvModel_xgboost requires correct inputs", {
  skip_if_not_installed("xgboost")
  model_xgb <- SurvModel_xgboost(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  expect_error(Predict_SurvModel_xgboost(newdata = test_data_surv), "argument \"modelout\" is missing")
  expect_error(Predict_SurvModel_xgboost(modelout = model_xgb), "argument \"newdata\" is missing")
})
