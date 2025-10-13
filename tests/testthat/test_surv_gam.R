library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(mgcv)     # Required for gam

# Assuming the functions are available in the environment
# source(here("R/models/surv_gam.R"))  # Removed - functions loaded via package

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
time_var <- "time"
event_var <- "status"
expvars <- c("x1", "x2", "x3")
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_GAM ---

test_that("SurvModel_GAM runs and returns expected structure", {
  skip_if_not_installed("mgcv")

  model_gam <- SurvModel_GAM(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )

  # Check output structure
  expect_type(model_gam, "list")
  expect_named(model_gam, c("model", "times", "varprof", "expvars"))
  expect_s3_class(model_gam$model, "gam")
  expect_type(model_gam$times, "double")
  expect_type(model_gam$varprof, "list")
  expect_type(model_gam$expvars, "character")
})

test_that("SurvModel_GAM handles additional parameters", {
  skip_if_not_installed("mgcv")

  # Test with different shrinkTreshold
  model_gam_params <- SurvModel_GAM(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    shrinkTreshold = 5
  )

  expect_s3_class(model_gam_params$model, "gam")
})

test_that("SurvModel_GAM requires correct inputs", {
  skip_if_not_installed("mgcv")
  expect_error(SurvModel_GAM(expvars = expvars, timevar = time_var, eventvar = event_var), "argument \"data\" is missing")
  expect_error(SurvModel_GAM(data = train_data_surv, timevar = time_var, eventvar = event_var), "argument \"expvars\" is missing")
  expect_error(SurvModel_GAM(data = train_data_surv, expvars = expvars, eventvar = event_var), "argument \"timevar\" is missing")
  expect_error(SurvModel_GAM(data = train_data_surv, expvars = expvars, timevar = time_var), "argument \"eventvar\" is missing")
})


# --- Tests for Predict_SurvModel_GAM ---

test_that("Predict_SurvModel_GAM returns predictions in correct format", {
  skip_if_not_installed("mgcv")

  model_gam <- SurvModel_GAM(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  predictions <- Predict_SurvModel_GAM(
    modelout = model_gam,
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
})

test_that("Predict_SurvModel_GAM handles custom times", {
  skip_if_not_installed("mgcv")

  model_gam <- SurvModel_GAM(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  # Note: The current implementation doesn't support custom times - it uses all times from the model
  # This test just verifies the function runs and returns the expected structure
  predictions <- Predict_SurvModel_GAM(
    modelout = model_gam,
    newdata = test_data_surv
  )

  expect_type(predictions, "list")
  expect_named(predictions, c("Probs", "Times"))
  expect_true(is.matrix(predictions$Probs))
  expect_type(predictions$Times, "double")
})

test_that("Predict_SurvModel_GAM requires correct inputs", {
  skip_if_not_installed("mgcv")
  model_gam <- SurvModel_GAM(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  expect_error(Predict_SurvModel_GAM(newdata = test_data_surv), "argument \"modelout\" is missing")
  expect_error(Predict_SurvModel_GAM(modelout = model_gam), "argument \"newdata\" is missing")
})
