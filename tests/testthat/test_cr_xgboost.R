# ==============================================================================
# Test Suite for Competing Risks XGBoost Model
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
# Tests for CRModel_xgboost - Basic Functionality
# ==============================================================================

test_that("CRModel_xgboost fits basic model", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10  # Small number for testing
  )

  # Check output structure
  expect_type(model, "list")
  expect_named(model, c("xgb_model", "times", "varprof", "model_type",
                       "expvars", "timevar", "eventvar", "failcode", "time_range", "feature_names"))

  # Check model type
  expect_equal(model$model_type, "cr_xgboost")
  expect_s3_class(model$xgb_model, "xgb.Booster")

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

  # Check feature names
  expect_true(is.character(model$feature_names))
  expect_true(length(model$feature_names) > 0)
})

test_that("CRModel_xgboost handles factor variables", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  expect_s3_class(model$xgb_model, "xgb.Booster")
  expect_true(is.list(model$varprof))
  expect_true(length(model$feature_names) > length(expvars_numeric))  # Should have dummy variables
})

test_that("CRModel_xgboost handles different failcodes", {
  skip_if_not_installed("xgboost")
  # Test failcode = 2
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 2,
    nrounds = 10
  )

  expect_equal(model$failcode, 2)
  expect_s3_class(model$xgb_model, "xgb.Booster")
})

# ==============================================================================
# Tests for CRModel_xgboost - Input Validation
# ==============================================================================

test_that("CRModel_xgboost validates inputs", {
  skip_if_not_installed("xgboost")
  # Missing data
  expect_error(CRModel_xgboost(expvars = expvars_numeric, timevar = "time", eventvar = "event"),
               "argument \"data\" is missing")

  # Missing expvars
  expect_error(CRModel_xgboost(data = train_data, timevar = "time", eventvar = "event"),
               "argument \"expvars\" is missing")

  # Missing timevar
  expect_error(CRModel_xgboost(data = train_data, expvars = expvars_numeric, eventvar = "event"),
               "argument \"timevar\" is missing")

  # Missing eventvar
  expect_error(CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time"),
               "argument \"eventvar\" is missing")

  # Invalid data type
  expect_error(CRModel_xgboost(data = "not data", expvars = expvars_numeric,
                              timevar = "time", eventvar = "event"),
               "'data' must be a data frame")

  # Empty expvars
  expect_error(CRModel_xgboost(data = train_data, expvars = character(0),
                              timevar = "time", eventvar = "event"),
               "'expvars' must be a non-empty character vector")

  # Invalid failcode
  expect_error(CRModel_xgboost(data = train_data, expvars = expvars_numeric,
                              timevar = "time", eventvar = "event", failcode = 0),
               "'failcode' must be a positive integer")

  # Non-existent timevar
  expect_error(CRModel_xgboost(data = train_data, expvars = expvars_numeric,
                              timevar = "nonexistent", eventvar = "event"),
               "'timevar' not found in data")

  # Non-existent eventvar
  expect_error(CRModel_xgboost(data = train_data, expvars = expvars_numeric,
                              timevar = "time", eventvar = "nonexistent"),
               "'eventvar' not found in data")

  # Non-existent expvars
  expect_error(CRModel_xgboost(data = train_data, expvars = c("nonexistent"),
                              timevar = "time", eventvar = "event"),
               "not found in data")
})

test_that("CRModel_xgboost handles missing data", {
  skip_if_not_installed("xgboost")
  # Create data with missing values
  train_data_na <- train_data
  train_data_na$x1[1:10] <- NA

  expect_warning(
    model <- CRModel_xgboost(
      data = train_data_na,
      expvars = expvars_numeric,
      timevar = "time",
      eventvar = "event",
      failcode = 1,
      nrounds = 10
    ),
    "Removed .* rows with missing values"
  )

  expect_s3_class(model$xgb_model, "xgb.Booster")
})

test_that("CRModel_xgboost handles insufficient data", {
  skip_if_not_installed("xgboost")
  # Create very small dataset
  small_data <- train_data[1:5, ]

  expect_error(CRModel_xgboost(
    data = small_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  ), "Insufficient data")
})

# ==============================================================================
# Tests for Predict_CRModel_xgboost - Basic Functionality
# ==============================================================================

test_that("Predict_CRModel_xgboost works with basic model", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)

  # Check output structure
  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))

  # Check CIFs
  expect_true(is.matrix(preds$CIFs))
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1))
  expect_equal(ncol(preds$CIFs), nrow(test_data))  # Columns = observations

  # Check Times
  expect_true(is.numeric(preds$Times))
  expect_true(length(preds$Times) > 0)
  expect_true(all(preds$Times >= 0))
  expect_equal(nrow(preds$CIFs), length(preds$Times))  # Rows = times
})

test_that("Predict_CRModel_xgboost handles custom times", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  custom_times <- c(1, 5, 10, 15)
  preds <- Predict_CRModel_xgboost(model, test_data, newtimes = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
})

test_that("Predict_CRModel_xgboost handles factor variables in newdata", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)

  expect_true(is.matrix(preds$CIFs))
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1))
})

# ==============================================================================
# Tests for Predict_CRModel_xgboost - Input Validation
# ==============================================================================

test_that("Predict_CRModel_xgboost validates inputs", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  # Missing modelout
  expect_error(Predict_CRModel_xgboost(newdata = test_data),
               "argument \"modelout\" is missing")

  # Missing newdata
  expect_error(Predict_CRModel_xgboost(modelout = model),
               "argument \"newdata\" is missing")

  # Invalid modelout
  expect_error(Predict_CRModel_xgboost(modelout = "not model", newdata = test_data),
               "'modelout' must be output from CRModel_xgboost")

  # Invalid newdata
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = "not data"),
               "'newdata' must be a data frame")

  # Missing variables in newdata
  test_data_missing <- test_data[, !(names(test_data) %in% c("x1"))]
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data_missing),
               "missing in newdata")

  # Invalid newtimes
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data, newtimes = "not numeric"),
               "'newtimes' must be a numeric vector")
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data, newtimes = c(-1, 1)),
               "'newtimes' must be a numeric vector of non-negative values")
})

test_that("Predict_CRModel_xgboost handles missing factor levels", {
  skip_if_not_installed("xgboost")
  # Create test data with missing factor level
  test_data_missing_level <- test_data
  test_data_missing_level$cat1 <- factor(c("A", "B"), levels = c("A", "B", "C"))  # Missing "C"

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  # Should work without error (missing levels get 0 in dummy coding)
  preds <- Predict_CRModel_xgboost(model, test_data_missing_level)
  expect_true(is.matrix(preds$CIFs))
})

# ==============================================================================
# Tests for CRModel_xgboost - Model Parameters
# ==============================================================================

test_that("CRModel_xgboost accepts XGBoost parameters", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    eta = 0.1,
    max_depth = 3,
    nrounds = 5,
    subsample = 0.8
  )

  expect_s3_class(model$xgb_model, "xgb.Booster")
})

test_that("CRModel_xgboost handles verbose parameter", {
  skip_if_not_installed("xgboost")
  # Should work with verbose = TRUE (no error expected)
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 5,
    verbose = TRUE
  )

  expect_s3_class(model$xgb_model, "xgb.Booster")
})

# ==============================================================================
# Tests for CIF Properties
# ==============================================================================

test_that("CR XGBoost CIFs have correct properties", {
  skip_if_not_installed("xgboost")
  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)

  # CIFs should be between 0 and 1
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1))

  # CIFs should be non-decreasing over time for each observation
  for (i in seq_len(ncol(preds$CIFs))) {
    cif_curve <- preds$CIFs[, i]
    expect_true(all(diff(cif_curve) >= -1e-10))  # Allow small numerical errors
  }

  # Check that times are sorted
  expect_true(all(diff(preds$Times) >= 0))
})

# ==============================================================================
# Tests for Model Comparison with Other CR Models
# ==============================================================================

test_that("CR XGBoost has similar interface to other CR models", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("mgcv")

  # Fit CR XGBoost
  xgb_model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1,
    nrounds = 10
  )

  # Fit CR GAM for comparison
  gam_model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    failcode = 1
  )

  # Both should have similar output structure
  expect_named(xgb_model, c("xgb_model", "times", "varprof", "model_type",
                           "expvars", "timevar", "eventvar", "failcode", "time_range", "feature_names"))
  expect_named(gam_model, c("gam_model", "times", "varprof", "model_type",
                           "expvars", "timevar", "eventvar", "failcode", "time_range"))

  # Both should produce valid predictions
  xgb_preds <- Predict_CRModel_xgboost(xgb_model, test_data)
  gam_preds <- Predict_CRModel_GAM(gam_model, test_data)

  expect_true(is.matrix(xgb_preds$CIFs))
  expect_true(is.matrix(gam_preds$CIFs))
  expect_true(all(xgb_preds$CIFs >= 0 & xgb_preds$CIFs <= 1))
  expect_true(all(gam_preds$CIFs >= 0 & gam_preds$CIFs <= 1))
})