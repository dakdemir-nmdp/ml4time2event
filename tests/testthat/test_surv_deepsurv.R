# ==============================================================================
# Test Suite for DeepSurv Neural Network Survival Model
# ==============================================================================

library(testthat)
library(survival)
library(nnet)

# Source required files
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/general_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_deepsurv.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_ensemble.R")

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated survival data for testing
set.seed(42)
n_train <- 200
n_test <- 50

# Training data
train_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = rbinom(n_train, 1, 0.7),
  x1 = rnorm(n_train),
  x2 = rnorm(n_train),
  x3 = rnorm(n_train, mean = 1),
  x4 = rnorm(n_train, mean = -1),
  cat1 = factor(sample(c("A", "B", "C"), n_train, replace = TRUE)),
  cat2 = factor(sample(c("Low", "High"), n_train, replace = TRUE))
)

# Test data
test_data <- data.frame(
  time = rexp(n_test, rate = 0.1),
  event = rbinom(n_test, 1, 0.7),
  x1 = rnorm(n_test),
  x2 = rnorm(n_test),
  x3 = rnorm(n_test, mean = 1),
  x4 = rnorm(n_test, mean = -1),
  cat1 = factor(sample(c("A", "B", "C"), n_test, replace = TRUE),
                levels = c("A", "B", "C")),
  cat2 = factor(sample(c("Low", "High"), n_test, replace = TRUE),
                levels = c("Low", "High"))
)

# Define variables
expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# Tests for SurvModel_DeepSurv - Basic Functionality
# ==============================================================================

test_that("SurvModel_DeepSurv fits basic model", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  # Check output structure
  expect_s3_class(model, "ml4t2e_surv_deepsurv")
  expect_true(is.list(model))
  expect_named(model, c("model", "times", "varprof", "expvars", "factor_levels"))

  # Check model components
  expect_s3_class(model$model, "nnet")
  expect_true(is.numeric(model$times))
  expect_true(is.list(model$varprof))
  expect_equal(model$expvars, expvars_numeric)
})

test_that("SurvModel_DeepSurv handles factor variables", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  expect_s3_class(model, "ml4t2e_surv_deepsurv")
  expect_s3_class(model$model, "nnet")

  # Check varprof captures factor levels
  expect_true(is.table(model$varprof$cat1))
  expect_true(is.table(model$varprof$cat2))
  expect_setequal(names(model$varprof$cat1), c("A", "B", "C"))
  expect_setequal(names(model$varprof$cat2), c("Low", "High"))
})

test_that("SurvModel_DeepSurv accepts custom parameters", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    size = 5,
    decay = 0.1,
    maxit = 100
  )

  expect_s3_class(model, "ml4t2e_surv_deepsurv")
  expect_s3_class(model$model, "nnet")
  expect_equal(model$model$n[2], 5)  # Check hidden layer size
})

# ==============================================================================
# Tests for Predict_SurvModel_DeepSurv - Basic Functionality
# ==============================================================================

test_that("Predict_SurvModel_DeepSurv returns correct output structure", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_DeepSurv(model, test_data)

  # Check output structure
  expect_true(is.list(preds))
  expect_named(preds, c("Probs", "Times"))

  # Check dimensions
  expect_true(is.matrix(preds$Probs))
  expect_true(is.numeric(preds$Times))
  expect_equal(nrow(preds$Probs), length(preds$Times))
  expect_equal(ncol(preds$Probs), nrow(test_data))
})

test_that("Predict_SurvModel_DeepSurv predictions are valid probabilities", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_DeepSurv(model, test_data)

  # Check probabilities are between 0 and 1
  expect_true(all(preds$Probs >= 0))
  expect_true(all(preds$Probs <= 1))

  # Check survival probabilities are monotonically decreasing over time for each observation
  if (length(preds$Times) > 1) {
    all_non_increasing <- all(apply(preds$Probs, 2, function(col) all(diff(col) <= 1e-9)))
    expect_true(all_non_increasing)
  }
})

test_that("Predict_SurvModel_DeepSurv handles custom time points", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  custom_times <- c(1, 5, 10, 15)
  preds <- Predict_SurvModel_DeepSurv(model, test_data, newtimes = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$Probs), length(custom_times))
})

test_that("Predict_SurvModel_DeepSurv includes time 0", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_DeepSurv(model, test_data)

  expect_true(preds$Times[1] == 0)
})

# ==============================================================================
# Tests for DeepSurv in Ensemble Framework
# ==============================================================================

test_that("DeepSurv works in RunSurvModels ensemble", {
  fitted_models <- RunSurvModels(
    datatrain = train_data,
    ExpVars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    models = c("glmnet", "deepsurv")
  )

  # Check DeepSurv model was fitted
  expect_true(!is.null(fitted_models$deepsurv_Model))
  expect_s3_class(fitted_models$deepsurv_Model, "ml4t2e_surv_deepsurv")
})

test_that("DeepSurv works in PredictSurvModels ensemble", {
  fitted_models <- RunSurvModels(
    datatrain = train_data,
    ExpVars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    models = c("glmnet", "deepsurv")
  )

  predictions <- PredictSurvModels(
    models = fitted_models,
    newdata = test_data,
    newtimes = c(5, 10, 15)
  )

  # Check DeepSurv predictions are available
  expect_true("newprobsdeepsurv_Model" %in% names(predictions$ModelPredictions))
  expect_true(!is.null(predictions$NewProbs))
})

# ==============================================================================
# Tests for Input Validation
# ==============================================================================

test_that("SurvModel_DeepSurv validates inputs", {
  # Missing data
  expect_error(
    SurvModel_DeepSurv(expvars = expvars_numeric, timevar = "time", eventvar = "event")
  )

  # Missing required parameters
  expect_error(
    SurvModel_DeepSurv(data = train_data, timevar = "time", eventvar = "event")
  )
  expect_error(
    SurvModel_DeepSurv(data = train_data, expvars = expvars_numeric, eventvar = "event")
  )
  expect_error(
    SurvModel_DeepSurv(data = train_data, expvars = expvars_numeric, timevar = "time")
  )
})

test_that("Predict_SurvModel_DeepSurv validates inputs", {
  model <- SurvModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  # Missing newdata
  expect_error(Predict_SurvModel_DeepSurv(model))

  # Invalid newtimes
  expect_error(
    Predict_SurvModel_DeepSurv(model, test_data, newtimes = c(-1, 5, 10))
  )
})
