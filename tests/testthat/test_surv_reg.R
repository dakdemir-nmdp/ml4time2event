library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # Required for Surv and survreg

# Assuming the functions are available in the environment
# source(here("R/models/surv_reg.R"))  # Removed - functions loaded via package

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
time_var <- "time"
event_var <- "status"
expvars <- c("x1", "x2", "x3")
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_SurvReg ---

test_that("SurvModel_SurvReg runs and returns expected structure", {
  skip_if_not_installed("survival")

  # Test with default distribution (exponential)
  model_survreg <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )

  # Check output structure
  expect_type(model_survreg, "list")
  expect_named(model_survreg, c("survregOut", "times", "varprof"))
  expect_s3_class(model_survreg$survregOut, "survreg")
  expect_type(model_survreg$times, "double")
  expect_type(model_survreg$varprof, "list")
})

test_that("SurvModel_SurvReg handles different distributions", {
  skip_if_not_installed("survival")
  # Test with weibull distribution
  model_survreg_weibull <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    dist = "weibull"
  )
  expect_s3_class(model_survreg_weibull$survregOut, "survreg")
  expect_equal(model_survreg_weibull$survregOut$dist, "weibull")

  # Test with lognormal
  model_survreg_lnorm <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    dist = "lognormal"
  )
  expect_s3_class(model_survreg_lnorm$survregOut, "survreg")
  expect_equal(model_survreg_lnorm$survregOut$dist, "lognormal")
})

test_that("SurvModel_SurvReg requires correct inputs", {
  skip_if_not_installed("survival")
  expect_error(SurvModel_SurvReg(data = train_data_surv, expvars = expvars, timevar = time_var), "argument \"eventvar\" is missing")
  expect_error(SurvModel_SurvReg(expvars = expvars, timevar = time_var, eventvar = event_var), "argument \"data\" is missing")
})


# --- Tests for Predict_SurvModel_SurvReg ---

test_that("Predict_SurvModel_SurvReg returns predictions in correct format", {
  skip_if_not_installed("survival")

  model_survreg <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  predictions <- Predict_SurvModel_SurvReg(
    modelout = model_survreg,
    newdata = test_data_surv,
    newtimes = time_points_surv
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

test_that("Predict_SurvModel_SurvReg handles default times", {
  skip_if_not_installed("survival")

  model_survreg <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  predictions <- Predict_SurvModel_SurvReg(
    modelout = model_survreg,
    newdata = test_data_surv
  )

  expect_type(predictions, "list")
  expect_named(predictions, c("Probs", "Times"))
  expect_true(length(predictions$Times) > 0)
})

test_that("Predict_SurvModel_SurvReg requires correct inputs", {
  skip_if_not_installed("survival")
  model_survreg <- SurvModel_SurvReg(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var
  )
  expect_error(Predict_SurvModel_SurvReg(newdata = test_data_surv), "argument \"modelout\" is missing")
  expect_error(Predict_SurvModel_SurvReg(modelout = model_survreg), "argument \"newdata\" is missing")
})
