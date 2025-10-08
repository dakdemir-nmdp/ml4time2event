# ==============================================================================
# Test Suite for Competing Risks DeepSurv Model
# ==============================================================================

library(testthat)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated competing risks data for testing
set.seed(42)
n_train <- 200
n_test <- 50

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
expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# Tests for CRModel_DeepSurv - Basic Functionality
# ==============================================================================

test_that("CRModel_DeepSurv fits basic model", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,  # Small network for testing
    maxit = 50  # Limited iterations for testing
  )

  # Check output structure
  expect_type(model, "list")
  expect_named(model, c("model", "times", "varprof", "expvars",
                       "factor_levels", "failcode", "model_type"))

  # Check model type
  expect_equal(model$model_type, "cr_deepsurv")
  expect_s3_class(model$model, "nnet")

  # Check times
  expect_true(is.numeric(model$times))
  expect_true(length(model$times) > 0)
  expect_true(all(model$times >= 0))

  # Check varprof
  expect_true(is.list(model$varprof))
  expect_equal(length(model$varprof), length(expvars_numeric))

  # Check failcode
  expect_equal(model$failcode, 1)
})

test_that("CRModel_DeepSurv handles factor variables", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  expect_s3_class(model$model, "nnet")
  expect_true(is.list(model$varprof))
})

test_that("CRModel_DeepSurv handles different failcodes", {
  # Test failcode = 2
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 2,
    size = 3,
    maxit = 50
  )

  expect_equal(model$failcode, 2)
  expect_s3_class(model$model, "nnet")
})

# ==============================================================================
# Tests for CRModel_DeepSurv - Input Validation
# ==============================================================================

test_that("CRModel_DeepSurv validates inputs", {
  # Missing data
  expect_error(CRModel_DeepSurv(expvars = expvars_numeric, timevar = "time", eventvar = "event"),
               "argument \"data\" is missing")

  # Missing expvars
  expect_error(CRModel_DeepSurv(data = train_data, timevar = "time", eventvar = "event"),
               "argument \"expvars\" is missing")

  # Missing timevar
  expect_error(CRModel_DeepSurv(data = train_data, expvars = expvars_numeric, eventvar = "event"),
               "argument \"timevar\" is missing")

  # Missing eventvar
  expect_error(CRModel_DeepSurv(data = train_data, expvars = expvars_numeric, timevar = "time"),
               "argument \"eventvar\" is missing")

  # Invalid data type
  expect_error(CRModel_DeepSurv(data = "not data", expvars = expvars_numeric,
                               timevar = "time", eventvar = "event"),
               "'data' must be a data frame")

  # Empty expvars
  expect_error(CRModel_DeepSurv(data = train_data, expvars = character(0),
                               timevar = "time", eventvar = "event"),
               "'expvars' must be a non-empty character vector")

  # Invalid failcode
  expect_error(CRModel_DeepSurv(data = train_data, expvars = expvars_numeric,
                               timevar = "time", eventvar = "event", failcode = 0),
               "'failcode' must be a positive integer")

  # Missing timevar column
  expect_error(CRModel_DeepSurv(data = train_data, expvars = expvars_numeric,
                               timevar = "missing_time", eventvar = "event"),
               "'timevar' not found in data")

  # Missing eventvar column
  expect_error(CRModel_DeepSurv(data = train_data, expvars = expvars_numeric,
                               timevar = "time", eventvar = "missing_event"),
               "'eventvar' not found in data")

  # Missing expvar column
  expect_error(CRModel_DeepSurv(data = train_data, expvars = c(expvars_numeric, "missing_var"),
                               timevar = "time", eventvar = "event"),
               "not found in data")
})

test_that("CRModel_DeepSurv handles insufficient data", {
  small_data <- train_data[1:5, ]  # Very small dataset

  expect_error(CRModel_DeepSurv(data = small_data, expvars = expvars_numeric,
                               timevar = "time", eventvar = "event"),
               "Insufficient data after removing missing values")
})

test_that("CRModel_DeepSurv handles no events of interest", {
  # Create data with no events of failcode = 1
  no_event_data <- train_data
  no_event_data$event[no_event_data$event == 1] <- 2  # Change all events to competing

  expect_error(CRModel_DeepSurv(data = no_event_data, expvars = expvars_numeric,
                               timevar = "time", eventvar = "event", failcode = 1),
               "No events of type 1 in training data")
})

# ==============================================================================
# Tests for Predict_CRModel_DeepSurv - Basic Functionality
# ==============================================================================

test_that("Predict_CRModel_DeepSurv returns correct output structure", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  preds <- Predict_CRModel_DeepSurv(model, test_data)

  # Check output structure
  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))

  # Check CIFs matrix
  expect_true(is.matrix(preds$CIFs))
  expect_equal(ncol(preds$CIFs), nrow(test_data))

  # Check Times
  expect_true(is.numeric(preds$Times))
  expect_true(all(preds$Times >= 0))

  # Check CIF values (should be between 0 and 1)
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1, na.rm = TRUE))
})

test_that("Predict_CRModel_DeepSurv CIFs are monotonically non-decreasing", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  preds <- Predict_CRModel_DeepSurv(model, test_data)

  # For each observation, CIF should be non-decreasing over time
  for (i in seq_len(ncol(preds$CIFs))) {
    cif_curve <- preds$CIFs[, i]
    diffs <- diff(cif_curve)
    expect_true(all(diffs >= -1e-10))  # Allow small numerical tolerance
  }
})

test_that("Predict_CRModel_DeepSurv handles custom time points", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  custom_times <- c(1, 5, 10, 20, 50)
  preds <- Predict_CRModel_DeepSurv(model, test_data, newtimes = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
})

test_that("Predict_CRModel_DeepSurv includes time 0", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  preds <- Predict_CRModel_DeepSurv(model, test_data)

  # Time 0 should be included
  expect_true(0 %in% preds$Times)

  # CIF at time 0 should be 0
  time_0_idx <- which(preds$Times == 0)
  expect_true(all(abs(preds$CIFs[time_0_idx, ]) < 1e-6))
})

# ==============================================================================
# Tests for Predict_CRModel_DeepSurv - Input Validation
# ==============================================================================

test_that("Predict_CRModel_DeepSurv validates inputs", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  # Missing modelout
  expect_error(Predict_CRModel_DeepSurv(newdata = test_data),
               "argument \"modelout\" is missing")

  # Missing newdata
  expect_error(Predict_CRModel_DeepSurv(modelout = model),
               "argument \"newdata\" is missing")

  # Invalid modelout
  expect_error(Predict_CRModel_DeepSurv(modelout = "not a model", newdata = test_data),
               "'modelout' must be output from CRModel_DeepSurv")

  # Invalid newdata
  expect_error(Predict_CRModel_DeepSurv(modelout = model, newdata = "not data"),
               "'newdata' must be a data frame")

  # Missing variables in newdata
  incomplete_data <- test_data[, -which(names(test_data) == "x1")]
  expect_error(Predict_CRModel_DeepSurv(modelout = model, newdata = incomplete_data),
               "missing in newdata")

  # Invalid newtimes
  expect_error(Predict_CRModel_DeepSurv(modelout = model, newdata = test_data,
                                       newtimes = c(-1, 1)),
               "'newtimes' must be a numeric vector of non-negative values")
})

# ==============================================================================
# Tests for Factor Level Handling
# ==============================================================================

test_that("Predict_CRModel_DeepSurv handles matching factor levels", {
  model <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    maxit = 50
  )

  # Test data with same factor levels
  preds <- Predict_CRModel_DeepSurv(model, test_data)

  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(all(!is.na(preds$CIFs)))
})

# ==============================================================================
# Tests for Model Parameters
# ==============================================================================

test_that("CRModel_DeepSurv accepts different network sizes", {
  # Small network
  model_small <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 2,
    maxit = 50
  )

  # Larger network
  model_large <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 5,
    maxit = 50
  )

  expect_s3_class(model_small$model, "nnet")
  expect_s3_class(model_large$model, "nnet")

  # Both should produce valid predictions
  preds_small <- Predict_CRModel_DeepSurv(model_small, test_data)
  preds_large <- Predict_CRModel_DeepSurv(model_large, test_data)

  expect_true(all(preds_small$CIFs >= 0 & preds_small$CIFs <= 1))
  expect_true(all(preds_large$CIFs >= 0 & preds_large$CIFs <= 1))
})

test_that("CRModel_DeepSurv handles regularization", {
  # No regularization
  model_no_reg <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    decay = 0,
    maxit = 50
  )

  # With regularization
  model_reg <- CRModel_DeepSurv(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    size = 3,
    decay = 0.1,
    maxit = 50
  )

  expect_s3_class(model_no_reg$model, "nnet")
  expect_s3_class(model_reg$model, "nnet")
})