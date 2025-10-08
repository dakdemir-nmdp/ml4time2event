# ==============================================================================
# Test Suite for Competing Risks GAM Model
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
expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# Tests for CRModel_GAM - Basic Functionality
# ==============================================================================

test_that("CRModel_GAM fits basic model", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  # Check output structure
  expect_type(model, "list")
  expect_named(model, c("gam_model", "times", "varprof", "model_type",
                       "expvars", "timevar", "eventvar", "failcode", "time_range"))

  # Check model type
  expect_equal(model$model_type, "cr_gam")
  expect_s3_class(model$gam_model, "gam")

  # Check time range
  expect_true(is.numeric(model$time_range))
  expect_equal(length(model$time_range), 2)
  expect_true(model$time_range[1] == 0)
  expect_true(model$time_range[2] > 0)

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

test_that("CRModel_GAM handles factor variables", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  expect_s3_class(model$gam_model, "gam")
  expect_true(is.list(model$varprof))
})

test_that("CRModel_GAM handles different failcodes", {
  skip_if_not_installed("mgcv")
  # Test failcode = 2
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 2
  )

  expect_equal(model$failcode, 2)
  expect_s3_class(model$gam_model, "gam")
})

# ==============================================================================
# Tests for CRModel_GAM - Input Validation
# ==============================================================================

test_that("CRModel_GAM validates inputs", {
  skip_if_not_installed("mgcv")
  # Missing data
  expect_error(CRModel_GAM(expvars = expvars_numeric, timevar = "time", eventvar = "event"),
               "argument \"data\" is missing")

  # Missing expvars
  expect_error(CRModel_GAM(data = train_data, timevar = "time", eventvar = "event"),
               "argument \"expvars\" is missing")

  # Missing timevar
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric, eventvar = "event"),
               "argument \"timevar\" is missing")

  # Missing eventvar
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric, timevar = "time"),
               "argument \"eventvar\" is missing")

  # Invalid data type
  expect_error(CRModel_GAM(data = "not data", expvars = expvars_numeric,
                          timevar = "time", eventvar = "event"),
               "'data' must be a data frame")

  # Empty expvars
  expect_error(CRModel_GAM(data = train_data, expvars = character(0),
                          timevar = "time", eventvar = "event"),
               "'expvars' must be a non-empty character vector")

  # Invalid failcode
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                          timevar = "time", eventvar = "event", failcode = 0),
               "'failcode' must be a positive integer")

  # Missing timevar column
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                          timevar = "missing_time", eventvar = "event"),
               "'timevar' not found in data")

  # Missing eventvar column
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                          timevar = "time", eventvar = "missing_event"),
               "'eventvar' not found in data")

  # Missing expvar column
  expect_error(CRModel_GAM(data = train_data, expvars = c(expvars_numeric, "missing_var"),
                          timevar = "time", eventvar = "event"),
               "not found in data")
})

test_that("CRModel_GAM handles insufficient data", {
  skip_if_not_installed("mgcv")
  small_data <- train_data[1:5, ]  # Very small dataset

  expect_error(CRModel_GAM(data = small_data, expvars = expvars_numeric,
                          timevar = "time", eventvar = "event"),
               "Insufficient data after removing missing values")
})

test_that("CRModel_GAM handles no events of interest", {
  skip_if_not_installed("mgcv")
  # Create data with no events of failcode = 1
  no_event_data <- train_data
  no_event_data$event[no_event_data$event == 1] <- 2  # Change all events to competing

  expect_error(CRModel_GAM(data = no_event_data, expvars = expvars_numeric,
                          timevar = "time", eventvar = "event", failcode = 1),
               "No events of type 1 in training data")
})

# ==============================================================================
# Tests for Predict_CRModel_GAM - Basic Functionality
# ==============================================================================

test_that("Predict_CRModel_GAM returns correct output structure", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  preds <- Predict_CRModel_GAM(model, test_data)

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

test_that("Predict_CRModel_GAM CIFs are monotonically non-decreasing", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  preds <- Predict_CRModel_GAM(model, test_data)

  # For each observation, CIF should be non-decreasing over time
  for (i in seq_len(ncol(preds$CIFs))) {
    cif_curve <- preds$CIFs[, i]
    diffs <- diff(cif_curve)
    expect_true(all(diffs >= -1e-10))  # Allow small numerical tolerance
  }
})

test_that("Predict_CRModel_GAM handles custom time points", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  custom_times <- c(1, 5, 10, 20, 50)
  preds <- Predict_CRModel_GAM(model, test_data, newtimes = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
})

test_that("Predict_CRModel_GAM includes time 0", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  preds <- Predict_CRModel_GAM(model, test_data)

  # Time 0 should be included
  expect_true(0 %in% preds$Times)

  # CIF at time 0 should be 0
  time_0_idx <- which(preds$Times == 0)
  expect_true(all(abs(preds$CIFs[time_0_idx, ]) < 1e-6))
})

# ==============================================================================
# Tests for Predict_CRModel_GAM - Input Validation
# ==============================================================================

test_that("Predict_CRModel_GAM validates inputs", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  # Missing modelout
  expect_error(Predict_CRModel_GAM(newdata = test_data),
               "argument \"modelout\" is missing")

  # Missing newdata
  expect_error(Predict_CRModel_GAM(modelout = model),
               "argument \"newdata\" is missing")

  # Invalid modelout
  expect_error(Predict_CRModel_GAM(modelout = "not a model", newdata = test_data),
               "'modelout' must be output from CRModel_GAM")

  # Invalid newdata
  expect_error(Predict_CRModel_GAM(modelout = model, newdata = "not data"),
               "'newdata' must be a data frame")

  # Missing variables in newdata
  incomplete_data <- test_data[, -which(names(test_data) == "x1")]
  expect_error(Predict_CRModel_GAM(modelout = model, newdata = incomplete_data),
               "missing in newdata")

  # Invalid newtimes
  expect_error(Predict_CRModel_GAM(modelout = model, newdata = test_data,
                                  newtimes = c(-1, 1)),
               "'newtimes' must be a numeric vector of non-negative values")
})

# ==============================================================================
# Tests for Factor Level Handling
# ==============================================================================

test_that("Predict_CRModel_GAM handles matching factor levels", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  # Test data with same factor levels
  preds <- Predict_CRModel_GAM(model, test_data)

  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(all(!is.na(preds$CIFs)))
})

# ==============================================================================
# Tests for Model Parameters
# ==============================================================================

test_that("CRModel_GAM handles different shrinkage thresholds", {
  skip_if_not_installed("mgcv")
  # Low threshold (more variables get smoothed)
  model_low <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    shrinkTreshold = 2
  )

  # High threshold (fewer variables get smoothed)
  model_high <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    shrinkTreshold = 20
  )

  expect_s3_class(model_low$gam_model, "gam")
  expect_s3_class(model_high$gam_model, "gam")

  # Both should produce valid predictions
  preds_low <- Predict_CRModel_GAM(model_low, test_data)
  preds_high <- Predict_CRModel_GAM(model_high, test_data)

  expect_true(all(preds_low$CIFs >= 0 & preds_low$CIFs <= 1))
  expect_true(all(preds_high$CIFs >= 0 & preds_high$CIFs <= 1))
})