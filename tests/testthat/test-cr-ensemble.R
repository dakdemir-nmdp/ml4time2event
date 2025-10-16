# ==============================================================================
# Test Suite for Competing Risks Ensemble Models
# ==============================================================================

library(testthat)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated competing risks data for testing
set.seed(42)
n_train <- 50  # Reduced for faster testing
n_test <- 20   # Reduced for faster testing

# Training data - competing risks format
train_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)), # 0=censored, 1=cause1, 2=cause2
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
  event = sample(0:2, n_test, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
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
expvars <- c("x1", "x2", "x3", "cat1", "cat2")
time_points <- c(1, 5, 10, 15)

# ==============================================================================
# Tests for RunCRModels - Basic Functionality
# ==============================================================================

test_that("RunCRModels fits basic models", {
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("fastcmprsk")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,  # Small number for testing
    varsel = FALSE
  )

  # Check output structure
  expect_type(models, "list")
  expect_s3_class(models, "CREnsemble")
  expect_true("model_status" %in% names(models))
  expect_named(models, c("input", "model_status", "RF_Model", "RF_Model2", "Cox_Model"), ignore.order = TRUE)

  # Check input structure
  expect_type(models$input, "list")
  expect_named(models$input, c("ExpVars", "ExpVars2", "timevar", "eventvar"))

  # Check that RF models are always fitted
  expect_true(!is.null(models$RF_Model))
  expect_true(!is.null(models$RF_Model2))

  # Check requested model
  expect_true(!is.null(models$Cox_Model))
})

test_that("RunCRModels handles multiple models", {
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("fastcmprsk")
  skip_if_not_installed("BART")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("FG", "cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  expect_named(models, c("input", "model_status", "RF_Model", "RF_Model2", "FG_Model", "Cox_Model"), ignore.order = TRUE)
  expect_true(!is.null(models$FG_Model))
  expect_true(!is.null(models$Cox_Model))
})

test_that("RunCRModels handles variable selection", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = TRUE
  )

  # With variable selection, ExpVars2 should be a subset
  expect_true(length(models$input$ExpVars2) <= length(models$input$ExpVars))
  expect_true(all(models$input$ExpVars2 %in% models$input$ExpVars))
})

test_that("RunCRModels handles model failures gracefully", {
  skip_if_not_installed("randomForestSRC")

  # Create data that might cause some models to fail
  bad_train_data <- train_data
  bad_train_data$time[1:5] <- NA  # Introduce NAs

  models <- RunCRModels(
    datatrain = bad_train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Should still return a valid structure even if some models fail
  expect_type(models, "list")
  expect_true("model_status" %in% names(models))
})

# ==============================================================================
# Tests for PredictCRModels - Basic Functionality
# ==============================================================================

test_that("PredictCRModels works with basic models", {
  skip_if_not_installed("randomForestSRC")

  # Fit models
  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Make predictions
  preds <- PredictCRModels(models, test_data, time_points)

  # Check output structure
  expect_type(preds, "list")
  expect_true(all(c("ModelPredictions", "NewProbs", "models_used", "ensemble_method") %in% names(preds)))

  # Check ensemble predictions
  expect_true(is.matrix(preds$NewProbs))
  expect_equal(nrow(preds$NewProbs), length(time_points))
  expect_equal(ncol(preds$NewProbs), nrow(test_data))
  expect_true(all(preds$NewProbs >= 0 & preds$NewProbs <= 1))
  expect_equal(preds$ensemble_method, "average")
  expect_true(length(preds$models_used) > 0)

  # Check individual model predictions
  expect_type(preds$ModelPredictions, "list")
  expect_true(length(preds$ModelPredictions) > 0)
})

test_that("PredictCRModels handles multiple models", {
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("fastcmprsk")

  # Fit models
  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("FG", "cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Make predictions
  preds <- PredictCRModels(models, test_data, time_points)

  expect_true(length(preds$ModelPredictions) >= 4)  # RF, RF2, FG, Cox
  expect_true(is.matrix(preds$NewProbs))
  expect_true(all(preds$NewProbs >= 0 & preds$NewProbs <= 1))
  expect_equal(preds$ensemble_method, "average")
})

test_that("PredictCRModels handles missing models gracefully", {
  skip_if_not_installed("randomForestSRC")

  # Create a models object with some NULL entries
  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Manually set one model to NULL to simulate failure
  models$Cox_Model <- NULL

  # Should still work
  preds <- PredictCRModels(models, test_data, time_points)
  expect_type(preds, "list")
  expect_true(is.matrix(preds$NewProbs))
  expect_equal(preds$ensemble_method, "average")
})

# ==============================================================================
# Tests for Ensemble Logic
# ==============================================================================

test_that("Ensemble averaging produces valid CIFs", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  preds <- PredictCRModels(models, test_data, time_points)

  # CIFs should be between 0 and 1
  expect_true(all(preds$NewProbs >= 0 & preds$NewProbs <= 1))

  # CIFs should be non-decreasing over time for each observation
  for (i in seq_len(ncol(preds$NewProbs))) {
    cif_curve <- preds$NewProbs[, i]
    expect_true(all(diff(cif_curve) >= -1e-10))  # Allow small numerical errors
  }
})

test_that("Ensemble handles edge cases", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = character(0),  # No additional models
    ntreeRF = 10,
    varsel = FALSE
  )

  preds <- PredictCRModels(models, test_data, time_points)

  # Should still work with just RF models
  expect_true(is.matrix(preds$NewProbs))
  expect_true(all(preds$NewProbs >= 0 & preds$NewProbs <= 1))
})

# ==============================================================================
# Tests for Input Validation
# ==============================================================================

test_that("RunCRModels validates inputs", {
  # Missing datatrain
  expect_error(RunCRModels(ExpVars = expvars, timevar = "time", eventvar = "event"))

  # Missing ExpVars
  expect_error(RunCRModels(datatrain = train_data, timevar = "time", eventvar = "event"))

  # Missing timevar
  expect_error(RunCRModels(datatrain = train_data, ExpVars = expvars, eventvar = "event"))

  # Missing eventvar
  expect_error(RunCRModels(datatrain = train_data, ExpVars = expvars, timevar = "time"))

  # Invalid datatrain
  expect_error(RunCRModels(datatrain = "not data", ExpVars = expvars, timevar = "time", eventvar = "event"))

  # Empty ExpVars
  expect_error(RunCRModels(datatrain = train_data, ExpVars = character(0), timevar = "time", eventvar = "event"))
})

test_that("PredictCRModels validates inputs", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = train_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Missing models
  expect_error(PredictCRModels(newdata = test_data, newtimes = time_points))

  # Missing newdata
  expect_error(PredictCRModels(models = models, newtimes = time_points))

  # Invalid models
  expect_error(PredictCRModels(models = "not models", newdata = test_data, newtimes = time_points))

  # Invalid newdata
  expect_error(PredictCRModels(models = models, newdata = "not data", newtimes = time_points))
})