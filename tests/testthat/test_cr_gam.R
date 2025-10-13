# ================================
# Test Suite for CRModel_GAM
# ================================

library(testthat)

# -------------------------------
# Data Setup
# -------------------------------

set.seed(42)
n_train <- 60
n_test <- 25

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

# -------------------------------
# CRModel_GAM - Core Behaviour
# -------------------------------

test_that("CRModel_GAM fits basic model", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  expect_s3_class(model$gam_model, "gam")
  expect_type(model$gam_models_all_causes, "list")
  expect_named(model, c(
    "gam_model", "gam_models_all_causes", "event_codes", "event_codes_numeric",
    "default_event_code", "default_event_code_numeric", "times", "varprof",
    "model_type", "expvars", "timevar", "eventvar", "time_range"
  ))

  expect_equal(model$model_type, "cr_gam")
  expect_identical(model$event_codes, c("1", "2"))
  expect_identical(model$event_codes_numeric, c(1, 2))
  expect_identical(model$default_event_code, "1")
  expect_identical(model$default_event_code_numeric, 1)
  expect_true(is.numeric(model$times))
  expect_true(all(model$times >= 0))
  expect_length(model$varprof, length(expvars_numeric))
  expect_equal(model$time_range[1], 0)
  expect_true(model$time_range[2] > 0)
})

test_that("CRModel_GAM handles factor variables", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = 1
  )

  expect_s3_class(model$gam_model, "gam")
  expect_length(model$varprof, length(expvars_all))
})

test_that("CRModel_GAM respects event code ordering", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(2, 1)
  )

  expect_identical(model$event_codes, c("2", "1"))
  expect_identical(model$default_event_code, "2")
})

# -------------------------------
# CRModel_GAM - Validation
# -------------------------------

test_that("CRModel_GAM validates inputs", {
  skip_if_not_installed("mgcv")

  expect_error(CRModel_GAM(expvars = expvars_numeric, timevar = "time", eventvar = "event"),
               "argument \"data\" is missing")
  expect_error(CRModel_GAM(data = train_data, timevar = "time", eventvar = "event"),
               "argument \"expvars\" is missing")
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric, eventvar = "event"),
               "argument \"timevar\" is missing")
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric, timevar = "time"),
               "argument \"eventvar\" is missing")

  expect_error(CRModel_GAM(data = "not data", expvars = expvars_numeric,
                           timevar = "time", eventvar = "event"),
               "'data' must be a data frame")
  expect_error(CRModel_GAM(data = train_data, expvars = character(0),
                           timevar = "time", eventvar = "event"),
               "'expvars' must be a non-empty character vector")
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                           timevar = "missing", eventvar = "event"),
               "'timevar' not found in data")
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                           timevar = "time", eventvar = "missing"),
               "'eventvar' not found in data")
  expect_error(CRModel_GAM(data = train_data, expvars = c(expvars_numeric, "missing"),
                           timevar = "time", eventvar = "event"),
               "not found in data")

  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                           timevar = "time", eventvar = "event", event_codes = "A"),
               "not present in the data")
  expect_error(CRModel_GAM(data = train_data, expvars = expvars_numeric,
                           timevar = "time", eventvar = "event", event_codes = 5),
               "not present in the data")
})

test_that("CRModel_GAM handles insufficient data", {
  skip_if_not_installed("mgcv")
  tiny <- train_data[1:5, ]

  expect_error(
    CRModel_GAM(data = tiny, expvars = expvars_numeric,
                timevar = "time", eventvar = "event"),
    "Insufficient data"
  )
})

test_that("CRModel_GAM requires events for requested code", {
  skip_if_not_installed("mgcv")
  no_event_data <- train_data
  no_event_data$event[no_event_data$event == 1] <- 2

  expect_error(
    CRModel_GAM(data = no_event_data, expvars = expvars_numeric,
                timevar = "time", eventvar = "event", event_codes = 1),
    "not present in the data"
  )
})

# -------------------------------
# Predict_CRModel_GAM - Behaviour
# -------------------------------

test_that("Predict_CRModel_GAM returns expected structure", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  preds <- Predict_CRModel_GAM(model, test_data)

  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))
  expect_true(is.matrix(preds$CIFs))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(is.numeric(preds$Times))
  expect_true(all(preds$Times >= 0))
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1, na.rm = TRUE))
})

test_that("Predict_CRModel_GAM outputs monotone CIFs", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  preds <- Predict_CRModel_GAM(model, test_data)

  apply(preds$CIFs, 2, function(curve) {
    diffs <- diff(curve)
    expect_true(all(diffs >= -1e-8))
  })
})

test_that("Predict_CRModel_GAM handles custom time grid", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  custom_times <- c(0, 1, 5, 10, 25)
  preds <- Predict_CRModel_GAM(model, test_data, newtimes = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
})

test_that("Predict_CRModel_GAM includes time zero with zero CIF", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  preds <- Predict_CRModel_GAM(model, test_data)
  zero_idx <- which(preds$Times == 0)
  expect_equal(length(zero_idx), 1)
  expect_true(all(abs(preds$CIFs[zero_idx, ]) < 1e-6))
})

test_that("Predict_CRModel_GAM validates inputs", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  expect_error(Predict_CRModel_GAM(newdata = test_data), "argument \"modelout\" is missing")
  expect_error(Predict_CRModel_GAM(modelout = model), "argument \"newdata\" is missing")
  expect_error(Predict_CRModel_GAM(modelout = "not a model", newdata = test_data),
               "'modelout' must be output")
  expect_error(Predict_CRModel_GAM(modelout = model, newdata = "not data"),
               "'newdata' must be a data frame")

  incomplete <- test_data[, setdiff(names(test_data), "x1")]
  expect_error(Predict_CRModel_GAM(modelout = model, newdata = incomplete),
               "missing in newdata")

  expect_error(Predict_CRModel_GAM(modelout = model, newdata = test_data, newtimes = c(-1, 1)),
               "must be a numeric vector of non-negative")
})

test_that("Predict_CRModel_GAM handles factor levels", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  preds <- Predict_CRModel_GAM(model, test_data)
  expect_equal(ncol(preds$CIFs), nrow(test_data))
})

test_that("Predict_CRModel_GAM supports alternate event selection", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  preds1 <- Predict_CRModel_GAM(model, test_data)
  preds2 <- Predict_CRModel_GAM(model, test_data, event_of_interest = 2)

  expect_equal(colnames(preds1$CIFs), colnames(preds2$CIFs))
  expect_false(identical(preds1$CIFs, preds2$CIFs))
})

test_that("Predict_CRModel_GAM errors on unknown event code", {
  skip_if_not_installed("mgcv")
  model <- CRModel_GAM(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2)
  )

  expect_error(Predict_CRModel_GAM(model, test_data, event_of_interest = 4),
               "Available event codes")
})

# -------------------------------
# CRModel_GAM - Parameter Checks
# -------------------------------

test_that("CRModel_GAM handles different shrinkage thresholds", {
  skip_if_not_installed("mgcv")
  model_low <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    shrinkTreshold = 2
  )

  model_high <- CRModel_GAM(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    event_codes = c(1, 2),
    shrinkTreshold = 25
  )

  expect_s3_class(model_low$gam_model, "gam")
  expect_s3_class(model_high$gam_model, "gam")

  preds_low <- Predict_CRModel_GAM(model_low, test_data)
  preds_high <- Predict_CRModel_GAM(model_high, test_data)

  expect_true(all(preds_low$CIFs >= 0 & preds_low$CIFs <= 1))
  expect_true(all(preds_high$CIFs >= 0 & preds_high$CIFs <= 1))
})