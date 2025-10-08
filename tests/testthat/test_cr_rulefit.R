library(testthat)
library(here)
library(survival) # For Surv object

context("Testing cr_rulefit functions")

# --- Test Data Setup ---
set.seed(123)
n_obs <- 100
cr_data <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)

train_indices <- 1:80
test_indices <- 81:100
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]
expvars <- c("x1", "x2", "x3")
time_points <- c(quantile(train_data$time[train_data$status != 0], 0.25),
                 median(train_data$time[train_data$status != 0]),
                 quantile(train_data$time[train_data$status != 0], 0.75))

# --- Tests for CRModel_rulefit ---

test_that("CRModel_rulefit fits basic model with correct interface", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(data = train_data, expvars = expvars,
                                  timevar = "time", eventvar = "status",
                                  failcode = 1, ntree = 10, nsample = 50)

  # Check output structure
  expect_type(model_rulefit, "list")
  expect_named(model_rulefit, c("rulefit_model", "times", "varprof", "model_type",
                               "expvars", "timevar", "eventvar", "failcode",
                               "time_range"))
  expect_s3_class(model_rulefit, "ml4t2e_cr_rulefit")

  # Check model type and basic properties
  expect_equal(model_rulefit$model_type, "cr_rulefit")
  expect_equal(model_rulefit$expvars, expvars)
  expect_equal(model_rulefit$timevar, "time")
  expect_equal(model_rulefit$eventvar, "status")
  expect_equal(model_rulefit$failcode, 1)
  expect_true(is.numeric(model_rulefit$times))
  expect_true(length(model_rulefit$times) > 0)
})

test_that("CRModel_rulefit validates inputs", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  # Test missing data
  expect_error(CRModel_rulefit(expvars = expvars, timevar = "time", eventvar = "status"),
               "argument \"data\" is missing, with no default")

  # Test invalid expvars
  expect_error(CRModel_rulefit(data = train_data, expvars = character(0),
                              timevar = "time", eventvar = "status"),
               "'expvars' must be a non-empty character vector")

  # Test missing timevar
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars,
                              timevar = "nonexistent", eventvar = "status"),
               "'timevar' not found in data")

  # Test missing eventvar
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars,
                              timevar = "time", eventvar = "nonexistent"),
               "'eventvar' not found in data")

  # Test invalid failcode
  expect_error(CRModel_rulefit(data = train_data, expvars = expvars,
                              timevar = "time", eventvar = "status", failcode = 3),
               "'failcode' must be 1 or 2")
})

# --- Tests for Predict_CRModel_rulefit ---

test_that("Predict_CRModel_rulefit returns predictions in correct format", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(data = train_data, expvars = expvars,
                                  timevar = "time", eventvar = "status",
                                  failcode = 1, ntree = 10, nsample = 50)

  predictions <- Predict_CRModel_rulefit(modelout = model_rulefit, newdata = test_data)

  # Check output structure
  expect_type(predictions, "list")
  expect_named(predictions, c("CIFs", "Times"))

  # Check dimensions: CIFs should be [times, observations]
  expect_true(is.matrix(predictions$CIFs))
  expect_equal(nrow(predictions$CIFs), length(predictions$Times))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions$CIFs >= 0 & predictions$CIFs <= 1, na.rm = TRUE))

  # Check that time 0 has CIF = 0
  expect_true(all(predictions$CIFs[1, ] == 0))
})

test_that("Predict_CRModel_rulefit handles custom time points", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(data = train_data, expvars = expvars,
                                  timevar = "time", eventvar = "status",
                                  failcode = 1, ntree = 10, nsample = 50)

  predictions <- Predict_CRModel_rulefit(modelout = model_rulefit, newdata = test_data,
                                        newtimes = time_points)

  # Check dimensions match requested time points
  expect_equal(length(predictions$Times), length(time_points))
  expect_equal(nrow(predictions$CIFs), length(time_points))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))

  # Check that times match requested times
  expect_equal(as.numeric(predictions$Times), as.numeric(sort(time_points)))
})

test_that("Predict_CRModel_rulefit validates inputs", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("pseudo")
  skip_if_not_installed("fastcmprsk")

  model_rulefit <- CRModel_rulefit(data = train_data, expvars = expvars,
                                  timevar = "time", eventvar = "status",
                                  failcode = 1, ntree = 10, nsample = 50)

  # Test invalid modelout
  expect_error(Predict_CRModel_rulefit(modelout = list(), newdata = test_data),
               "'modelout' must be output from CRModel_rulefit")

  # Test missing newdata
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit),
               "argument \"newdata\" is missing, with no default")

  # Test missing variables in newdata
  test_data_missing <- test_data[, -which(names(test_data) == "x1")]
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit, newdata = test_data_missing),
               "variables missing in newdata")

  # Test invalid newtimes
  expect_error(Predict_CRModel_rulefit(modelout = model_rulefit, newdata = test_data,
                                      newtimes = c(-1, 1)),
               "'newtimes' must be a numeric vector of non-negative values")
})
