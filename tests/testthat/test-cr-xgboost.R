# ==============================================================================
# Test Suite for Competing Risks XGBoost Model
# ==============================================================================

library(testthat)

# Load the package
devtools::load_all()

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(42)
n_train <- 50
n_test <- 20

train_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
  x1 = rnorm(n_train),
  x2 = rnorm(n_train),
  x3 = rnorm(n_train, mean = 1),
  x4 = rnorm(n_train, mean = -1),
  cat1 = factor(sample(c("A", "B", "C"), n_train, replace = TRUE)),
  cat2 = factor(sample(c("Low", "High"), n_train, replace = TRUE))
)

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

expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# CRModel_xgboost - Basic behaviour
# ==============================================================================

test_that("CRModel_xgboost fits basic model", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    nrounds = 10
  )

  expect_type(model, "list")
  expect_named(model, c(
    "xgb_model", "xgb_models_all_causes", "event_codes", "event_codes_numeric",
    "default_event_code", "times", "varprof", "model_type", "expvars",
    "timevar", "eventvar", "time_range", "feature_names"
  ))

  expect_equal(model$model_type, "cr_xgboost")
  expect_s3_class(model$xgb_model, "xgb.Booster")

  expect_equal(model$default_event_code, "1")
  expect_equal(model$event_codes, c("1", "2"))
  expect_equal(model$event_codes_numeric, c(1, 2))

  expect_true(is.numeric(model$time_range))
  expect_equal(length(model$time_range), 2)
  expect_true(model$time_range[1] >= 0)
  expect_true(model$time_range[2] > 0)

  expect_true(is.numeric(model$times))
  expect_true(length(model$times) > 0)
  expect_true(all(model$times >= 0))

  expect_true(is.list(model$varprof))
  expect_equal(length(model$varprof), length(expvars_numeric))

  expect_true(is.character(model$feature_names))
  expect_true(length(model$feature_names) > 0)
})

test_that("CRModel_xgboost handles factor predictors", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = NULL,
    nrounds = 10
  )

  expect_s3_class(model$xgb_model, "xgb.Booster")
  expect_true(length(model$feature_names) > length(expvars_numeric))
})

test_that("CRModel_xgboost respects event code ordering", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(2, 1),
    nrounds = 10
  )

  expect_equal(model$default_event_code, "2")
  expect_equal(model$event_codes, c("2", "1"))
})

# ==============================================================================
# CRModel_xgboost - Validation
# ==============================================================================

test_that("CRModel_xgboost validates inputs", {
  skip_if_not_installed("xgboost")

  expect_error(
    CRModel_xgboost(expvars = expvars_numeric, timevar = "time", eventvar = "event"),
    "argument \"data\" is missing"
  )

  expect_error(
    CRModel_xgboost(data = train_data, timevar = "time", eventvar = "event"),
    "argument \"expvars\" is missing"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, eventvar = "event"),
    "argument \"timevar\" is missing"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time"),
    "argument \"eventvar\" is missing"
  )

  expect_error(
    CRModel_xgboost(data = "not data", expvars = expvars_numeric, timevar = "time", eventvar = "event"),
    "'data' must be a data frame"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = character(0), timevar = "time", eventvar = "event"),
    "'expvars' must be a non-empty character vector"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time", eventvar = "event", event_codes = numeric(0)),
    "'event_codes' must be NULL or a non-empty vector"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time", eventvar = "event", event_codes = 0),
    "The following event_codes are not present in the data"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time", eventvar = "event", event_codes = "not-numeric"),
    "event_codes are not present in the data"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "missing", eventvar = "event"),
    "'timevar' not found in data"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = expvars_numeric, timevar = "time", eventvar = "missing"),
    "'eventvar' not found in data"
  )

  expect_error(
    CRModel_xgboost(data = train_data, expvars = "missing", timevar = "time", eventvar = "event"),
    "not found in data"
  )
})

test_that("CRModel_xgboost warns when removing missing cases", {
  skip_if_not_installed("xgboost")

  train_data_na <- train_data
  train_data_na$x1[1:10] <- NA

  expect_warning(
    CRModel_xgboost(
      data = train_data_na,
      expvars = expvars_numeric,
      timevar = "time",
      eventvar = "event",
      event_codes = c(1, 2),
      nrounds = 10
    ),
    "Removed .* rows with missing values"
  )
})

test_that("CRModel_xgboost requires sufficient rows", {
  skip_if_not_installed("xgboost")

  small_data <- train_data[1:5, ]

  expect_error(
    CRModel_xgboost(
      data = small_data,
      expvars = expvars_numeric,
      timevar = "time",
      eventvar = "event",
      event_codes = c(1, 2),
      nrounds = 10
    ),
    "Insufficient data"
  )
})

# ==============================================================================
# Predict_CRModel_xgboost - Behaviour
# ==============================================================================

test_that("Predict_CRModel_xgboost returns CIF matrix", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = NULL,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)

  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))
  expect_true(is.matrix(preds$CIFs))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(is.numeric(preds$Times))
  expect_equal(nrow(preds$CIFs), length(preds$Times))
  expect_true(all(diff(preds$Times) >= 0))

  expect_true(all(is.na(preds$CIFs) | (preds$CIFs >= 0 & preds$CIFs <= 1)))
})

test_that("Predict_CRModel_xgboost honours custom time grid", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = NULL,
    nrounds = 10
  )

  custom_times <- c(1, 5, 10, 15)
  preds <- Predict_CRModel_xgboost(model, test_data, new_times = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
})

test_that("Predict_CRModel_xgboost validates inputs", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    nrounds = 10
  )

  expect_error(Predict_CRModel_xgboost(newdata = test_data), "'modelout' is missing")
  expect_error(Predict_CRModel_xgboost(modelout = model), "'newdata' is missing")
  expect_error(Predict_CRModel_xgboost(modelout = "not model", newdata = test_data), "must be output from CRModel_xgboost")
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = "not data"), "Input 'newdata' must be a data frame.")

  test_data_missing <- test_data[, !(names(test_data) %in% c("x1"))]
  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data_missing), "The following 'expvars' are missing in 'newdata': x1")

  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data, new_times = "not numeric"),
               "'new_times' must be a numeric vector of non-negative values")

  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data, new_times = c(-1, 1)),
               "'new_times' must be a numeric vector of non-negative values")

  expect_error(Predict_CRModel_xgboost(modelout = model, newdata = test_data, event_of_interest = "99"),
               "Input 'event_of_interest' was not present in training data. Available event codes: 1, 2")
})

test_that("Predict_CRModel_xgboost supports factor predictors in newdata", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = NULL,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)
  expect_true(is.matrix(preds$CIFs))
})

test_that("Predict_CRModel_xgboost returns different CIFs for different events", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    nrounds = 10
  )

  preds_event1 <- Predict_CRModel_xgboost(model, test_data, event_of_interest = "1")
  preds_event2 <- Predict_CRModel_xgboost(model, test_data, event_of_interest = "2")

  expect_false(identical(preds_event1$CIFs, preds_event2$CIFs))
})

test_that("CR XGBoost CIF curves are non-decreasing", {
  skip_if_not_installed("xgboost")

  model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = NULL,
    nrounds = 10
  )

  preds <- Predict_CRModel_xgboost(model, test_data)

  for (i in seq_len(ncol(preds$CIFs))) {
    curve <- preds$CIFs[, i]
    curve <- curve[!is.na(curve)]
    if (length(curve) > 1) {
      expect_true(all(diff(curve) >= -1e-8))
    }
  }
})

# ==============================================================================
# Model comparison
# ==============================================================================

test_that("CR XGBoost interface aligns with CR GAM", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("mgcv")

  xgb_model <- CRModel_xgboost(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    nrounds = 10
  )

  gam_model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_of_interest = 1,
    event_codes = c(1, 2) # Ensure event_of_interest is part of event_codes
  )

  expect_named(xgb_model, c(
    "xgb_model", "xgb_models_all_causes", "event_codes", "event_codes_numeric",
    "default_event_code", "times", "varprof", "model_type", "expvars",
    "timevar", "eventvar", "time_range", "feature_names"
  ))

  expect_named(gam_model, c(
    "gam_model", "gam_models_all_causes", "event_codes", "event_codes_numeric",
    "default_event_code", "default_event_code_numeric", "times", "varprof",
    "model_type", "expvars", "timevar", "eventvar", "time_range"
  ))

  xgb_preds <- Predict_CRModel_xgboost(xgb_model, test_data)
  gam_preds <- Predict_CRModel_GAM(gam_model, test_data)

  expect_true(is.matrix(xgb_preds$CIFs))
  expect_true(is.matrix(gam_preds$CIFs))
  expect_true(all(is.na(xgb_preds$CIFs) | (xgb_preds$CIFs >= 0 & xgb_preds$CIFs <= 1)))
  expect_true(all(is.na(gam_preds$CIFs) | (gam_preds$CIFs >= 0 & gam_preds$CIFs <= 1)))
})