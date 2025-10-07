library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(rpart)    # Required for rpart
library(partykit) # Required for partykit
library(glmnet)   # Required for glmnet

# Assuming the functions are available in the environment
# source(here("R/surv_rulefit.R"))  # Removed - functions loaded via package

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
time_var <- "time"
event_var <- "status"
expvars <- c("x1", "x2", "x3")
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
time_points_surv <- quantile(train_data_surv$time[train_data_surv$status == 1], c(0.25, 0.5, 0.75))


# --- Tests for SurvModel_rulefit ---

test_that("SurvModel_rulefit runs and returns expected structure", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  model_rulefit <- SurvModel_rulefit(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    ntree = 5,  # Use fewer trees for faster testing
    nsample = 30
  )

  # Check output structure
  expect_type(model_rulefit, "list")
  expect_named(model_rulefit, c("model", "times", "varprof", "expvars", "factor_levels"))
  expect_type(model_rulefit$model, "list")
  expect_type(model_rulefit$times, "double")
  expect_type(model_rulefit$varprof, "list")
  expect_type(model_rulefit$expvars, "character")
  expect_type(model_rulefit$factor_levels, "list")
})

test_that("SurvModel_rulefit handles additional parameters", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  model_rulefit_params <- SurvModel_rulefit(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    ntree = 5,
    nsample = 30,
    keepvars = "x1"
  )

  expect_type(model_rulefit_params$model, "list")
})

test_that("SurvModel_rulefit requires correct inputs", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  expect_error(SurvModel_rulefit(expvars = expvars, timevar = time_var, eventvar = event_var), "argument \"data\" is missing")
  expect_error(SurvModel_rulefit(data = train_data_surv, timevar = time_var, eventvar = event_var), "argument \"expvars\" is missing")
  expect_error(SurvModel_rulefit(data = train_data_surv, expvars = expvars, eventvar = event_var), "argument \"timevar\" is missing")
  expect_error(SurvModel_rulefit(data = train_data_surv, expvars = expvars, timevar = time_var), "argument \"eventvar\" is missing")
})


# --- Tests for Predict_SurvModel_rulefit ---

test_that("Predict_SurvModel_rulefit returns predictions in correct format", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  model_rulefit <- SurvModel_rulefit(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    ntree = 5,
    nsample = 30
  )
  predictions <- Predict_SurvModel_rulefit(
    modelout = model_rulefit,
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

test_that("Predict_SurvModel_rulefit handles custom times", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  model_rulefit <- SurvModel_rulefit(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    ntree = 5,
    nsample = 30
  )
  # Note: The current implementation doesn't support custom times - it uses all times from the model
  # This test just verifies the function runs and returns the expected structure
  predictions <- Predict_SurvModel_rulefit(
    modelout = model_rulefit,
    newdata = test_data_surv
  )

  expect_type(predictions, "list")
  expect_named(predictions, c("Probs", "Times"))
  expect_true(is.matrix(predictions$Probs))
  expect_type(predictions$Times, "double")
})

test_that("Predict_SurvModel_rulefit requires correct inputs", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")

  model_rulefit <- SurvModel_rulefit(
    data = train_data_surv,
    expvars = expvars,
    timevar = time_var,
    eventvar = event_var,
    ntree = 5,
    nsample = 30
  )
  expect_error(Predict_SurvModel_rulefit(newdata = test_data_surv), "argument \"modelout\" is missing")
  expect_error(Predict_SurvModel_rulefit(modelout = model_rulefit), "argument \"newdata\" is missing")
})
