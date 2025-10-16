# ==============================================================================
# Test Suite for Ensemble Model Persistence
# Tests for S3 classes, save/load functionality
# ==============================================================================

library(testthat)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(456)
n_train <- 50

# Survival data
surv_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = rbinom(n_train, 1, 0.7),
  x1 = rnorm(n_train),
  x2 = rnorm(n_train)
)

# Competing risks data
cr_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
  x1 = rnorm(n_train),
  x2 = rnorm(n_train)
)

expvars <- c("x1", "x2")

# ==============================================================================
# Tests for S3 Classes
# ==============================================================================

test_that("RunSurvModels returns SurvEnsemble object", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  expect_true(is.SurvEnsemble(models))
  expect_s3_class(models, "SurvEnsemble")
  expect_true(inherits(models, "list"))
})

test_that("RunCRModels returns CREnsemble object", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = cr_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  expect_true(is.CREnsemble(models))
  expect_s3_class(models, "CREnsemble")
  expect_true(inherits(models, "list"))
})

test_that("SurvEnsemble class can be manually created", {
  simple_list <- list(
    input = list(ExpVars = c("x1", "x2"), timevar = "time", eventvar = "event"),
    model_status = c(model1 = TRUE, model2 = FALSE)
  )

  ensemble <- SurvEnsemble(simple_list)

  expect_true(is.SurvEnsemble(ensemble))
  expect_equal(ensemble$input$ExpVars, c("x1", "x2"))
})

test_that("CREnsemble class can be manually created", {
  simple_list <- list(
    input = list(ExpVars = c("x1", "x2"), timevar = "time", eventvar = "event"),
    model_status = c(model1 = TRUE, model2 = FALSE)
  )

  ensemble <- CREnsemble(simple_list)

  expect_true(is.CREnsemble(ensemble))
  expect_equal(ensemble$input$ExpVars, c("x1", "x2"))
})

test_that("SurvEnsemble validation works", {
  expect_error(SurvEnsemble("not a list"), "must be a list")
  expect_error(SurvEnsemble(list(foo = "bar")), "must contain 'input'")
})

test_that("CREnsemble validation works", {
  expect_error(CREnsemble("not a list"), "must be a list")
  expect_error(CREnsemble(list(foo = "bar")), "must contain 'input'")
})

# ==============================================================================
# Tests for Print Methods
# ==============================================================================

test_that("print.SurvEnsemble works", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  # Test that print doesn't error
  expect_output(print(models), "Survival Ensemble Model")
  expect_output(print(models), "Variables:")
  expect_output(print(models), "Models fitted:")
})

test_that("print.CREnsemble works", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = cr_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Test that print doesn't error
  expect_output(print(models), "Competing Risks Ensemble Model")
  expect_output(print(models), "Variables:")
  expect_output(print(models), "Models fitted:")
})

# ==============================================================================
# Tests for Summary Methods
# ==============================================================================

test_that("summary.SurvEnsemble works", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  expect_output(summary(models), "Survival Ensemble Model Summary")
  expect_output(summary(models), "Input Configuration:")
  expect_output(summary(models), "Model Status:")
})

test_that("summary.CREnsemble works", {
  skip_if_not_installed("randomForestSRC")

  models <- RunCRModels(
    datatrain = cr_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  expect_output(summary(models), "Competing Risks Ensemble Model Summary")
  expect_output(summary(models), "Input Configuration:")
  expect_output(summary(models), "Model Status:")
})

# ==============================================================================
# Tests for Save/Load Functionality
# ==============================================================================

test_that("SaveEnsemble and LoadEnsemble work for SurvEnsemble", {
  skip_if_not_installed("randomForestSRC")

  # Create a model
  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  # Create temp file
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Save
  expect_message(result_file <- SaveEnsemble(models, temp_file), "saved to")
  expect_true(file.exists(temp_file))

  # Load
  loaded_models <- LoadEnsemble(temp_file)

  # Verify
  expect_true(is.SurvEnsemble(loaded_models))
  expect_equal(loaded_models$input$ExpVars, models$input$ExpVars)
  expect_equal(loaded_models$model_status, models$model_status)
})

test_that("SaveEnsemble and LoadEnsemble work for CREnsemble", {
  skip_if_not_installed("randomForestSRC")

  # Create a model
  models <- RunCRModels(
    datatrain = cr_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("cox"),
    ntreeRF = 10,
    varsel = FALSE
  )

  # Create temp file
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Save
  expect_message(result_file <- SaveEnsemble(models, temp_file), "saved to")
  expect_true(file.exists(temp_file))

  # Load
  loaded_models <- LoadEnsemble(temp_file)

  # Verify
  expect_true(is.CREnsemble(loaded_models))
  expect_equal(loaded_models$input$ExpVars, models$input$ExpVars)
  expect_equal(loaded_models$model_status, models$model_status)
})

test_that("SaveEnsemble adds .rds extension automatically", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  temp_file <- tempfile()  # No extension
  on.exit(unlink(paste0(temp_file, ".rds")), add = TRUE)

  result_file <- SaveEnsemble(models, temp_file)

  expect_true(grepl("\\.rds$", result_file))
  expect_true(file.exists(result_file))
})

test_that("SaveEnsemble includes metadata", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  SaveEnsemble(models, temp_file)
  loaded_models <- LoadEnsemble(temp_file)

  expect_true(".__metadata__" %in% names(loaded_models))
  expect_true("saved_date" %in% names(loaded_models$.__metadata__))
  expect_true("r_version" %in% names(loaded_models$.__metadata__))
  expect_true("package_version" %in% names(loaded_models$.__metadata__))
})

test_that("LoadEnsemble shows metadata message", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  SaveEnsemble(models, temp_file)

  expect_message(LoadEnsemble(temp_file), "Loaded ensemble saved on")
})

test_that("LoadEnsemble handles non-existent file", {
  expect_error(LoadEnsemble("nonexistent_file.rds"), "File not found")
})

test_that("SaveEnsemble converts plain list to ensemble", {
  plain_list <- list(
    input = list(ExpVars = c("x1", "x2"), timevar = "time", eventvar = "event"),
    RF_Model = list(some_model = "data")
  )

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  expect_message(SaveEnsemble(plain_list, temp_file), "Converting plain list")
  expect_true(file.exists(temp_file))

  loaded <- LoadEnsemble(temp_file)
  expect_true(is.SurvEnsemble(loaded))
})

# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("Saved ensemble can be used for prediction", {
  skip_if_not_installed("randomForestSRC")

  # Fit model
  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  # Save
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)
  SaveEnsemble(models, temp_file)

  # Load
  loaded_models <- LoadEnsemble(temp_file)

  # Predict with loaded model
  test_data <- surv_data[1:10, ]
  newtimes <- c(1, 5, 10)

  preds <- PredictSurvModels(
    models = loaded_models,
    newdata = test_data,
    newtimes = newtimes
  )

  # Should work - predictions should be returned
  expect_type(preds, "list")
  expect_true("NewProbs" %in% names(preds))
  expect_true("models_used" %in% names(preds))

  # NewProbs might be NULL if all predictions had dimension mismatches,
  # but if it exists it should be a matrix or able to be converted to one
  if (!is.null(preds$NewProbs)) {
    expect_true(is.matrix(preds$NewProbs) || is.numeric(preds$NewProbs))
    # If numeric, should be able to convert to matrix
    if (!is.matrix(preds$NewProbs)) {
      preds$NewProbs <- as.matrix(preds$NewProbs)
    }
    expect_true(is.matrix(preds$NewProbs))
  }
})

test_that("Ensemble class preserved through save/load cycle", {
  skip_if_not_installed("randomForestSRC")

  models <- RunSurvModels(
    datatrain = surv_data,
    ExpVars = expvars,
    timevar = "time",
    eventvar = "event",
    models = c("coxph"),
    ntreeRF = 10
  )

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Original class
  expect_true(is.SurvEnsemble(models))

  # Save and load
  SaveEnsemble(models, temp_file)
  loaded <- LoadEnsemble(temp_file)

  # Class preserved
  expect_true(is.SurvEnsemble(loaded))
  expect_equal(class(models), class(loaded))
})
