library(testthat)
library(here)
library(survival) # For Surv object

context("Testing cr_rulefit functions")

# --- Test Data Setup ---
set.seed(123)
n_obs <- 70  # Reduced for faster testing
cr_data <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)

train_indices <- 1:50
test_indices <- 51:70
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]
expvars <- c("x1", "x2", "x3")
custom_times <- quantile(train_data$time[train_data$status != 0], probs = c(0.25, 0.5, 0.75))

# -----------------------------------------------------------------------------
# CRModel_rulefit - Core Behaviour
# -----------------------------------------------------------------------------

test_that("CRModel_rulefit fits basic model with event metadata", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = "time",
    eventvar = "status",
    event_codes = 1,
    ntree = 10,
    nsample = 50
  )

  expect_s3_class(model_rulefit, "ml4t2e_cr_rulefit")
  expect_named(model_rulefit, c(
    "rulefit_model", "times", "varprof", "model_type",
    "expvars", "timevar", "eventvar", "event_codes",
    "event_codes_numeric", "default_event_code", "default_event_code_numeric",
    "time_range"
  ))

  expect_equal(model_rulefit$model_type, "cr_rulefit")
  expect_equal(model_rulefit$expvars, expvars)
  expect_equal(model_rulefit$timevar, "time")
  expect_equal(model_rulefit$eventvar, "status")
  expect_identical(model_rulefit$event_codes, "1")
  expect_identical(model_rulefit$event_codes_numeric, 1)
  expect_identical(model_rulefit$default_event_code, "1")
  expect_identical(model_rulefit$default_event_code_numeric, 1)
  expect_true(all(model_rulefit$times >= 0))
  expect_equal(model_rulefit$time_range[1], 0)
})

test_that("CRModel_rulefit validates inputs", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  expect_error(CRModel_rulefit(expvars = expvars, timevar = "time", eventvar = "status"),
               "argument \"data\" is missing")
  expect_error(CRModel_rulefit(data = train_data, expvars = character(0), timevar = "time", eventvar = "status"),
               "`expvars` must be a non-empty character vector")
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars, timevar = "missing", eventvar = "status"),
               "`timevar` not found in data")
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars, timevar = "time", eventvar = "missing"),
               "`eventvar` not found in data")
  expect_error(CRModel_rulefit(data = train_data, expvars = c(expvars, "missing"), timevar = "time", eventvar = "status"),
               "not found in data")
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars, timevar = "time", eventvar = "status", event_codes = c(1, 2)),
               "`event_codes` must be a single value \\(one event code\\)")
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars, timevar = "time", eventvar = "status", event_codes = 5),
               "not present in the training data")
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars, timevar = "time", eventvar = "status", event_codes = "A"),
               "not present in the training data")
})

test_that("CRModel_rulefit requires events for selected code", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  no_event_data <- train_data
  no_event_data$status[no_event_data$status == 1] <- 2

  expect_error(
    CRModel_rulefit(data = no_event_data, expvars = expvars, timevar = "time", eventvar = "status", event_codes = 1),
    "not present in the training data"
  )
})

# -----------------------------------------------------------------------------
# Predict_CRModel_rulefit - Behaviour
# -----------------------------------------------------------------------------

test_that("Predict_CRModel_rulefit returns CIF matrix and times", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = "time",
    eventvar = "status",
    event_codes = 1,
    ntree = 10,
    nsample = 50
  )

  preds <- Predict_CRModel_rulefit(model_rulefit, test_data)

  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))
  expect_true(is.matrix(preds$CIFs))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1, na.rm = TRUE))
  expect_true(all(preds$Times >= 0))
})

test_that("Predict_CRModel_rulefit handles custom time grid", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = "time",
    eventvar = "status",
    event_codes = 1,
    ntree = 10,
    nsample = 50
  )

  preds <- Predict_CRModel_rulefit(model_rulefit, test_data, new_times = custom_times)
  expect_equal(length(preds$Times), length(custom_times))
  expect_equal(nrow(preds$CIFs), length(custom_times))
})

test_that("Predict_CRModel_rulefit validates inputs", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = "time",
    eventvar = "status",
    event_codes = 1,
    ntree = 10,
    nsample = 50
  )

  expect_error(Predict_CRModel_rulefit(modelout = list(), newdata = test_data),
               "'modelout' must be output from CRModel_rulefit")
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit),
               "'newdata' is missing")
  missing_var_data <- test_data[, setdiff(names(test_data), "x1"), drop = FALSE]
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit, newdata = missing_var_data),
               "missing in newdata")
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit, newdata = test_data, new_times = c(-1, 1)),
               "numeric vector of non-negative")
})

test_that("Predict_CRModel_rulefit enforces training event code", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(
    data = train_data,
    expvars = expvars,
    timevar = "time",
    eventvar = "status",
    event_codes = 1,
    ntree = 10,
    nsample = 50
  )

  expect_error(Predict_CRModel_rulefit(model_rulefit, test_data, event_of_interest = 2),
               "can only predict for the event they were trained on")
})
