library(testthat)
library(randomForestSRC) # Required for the functions being tested
library(survival) # For Surv object

context("Testing cr_random_forest functions")

# --- Test Data Setup ---
# Create a simple competing risks dataset
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
time_points <- c(quantile(train_data$time[train_data$status != 0], 0.25, na.rm = TRUE),
                 median(train_data$time[train_data$status != 0], na.rm = TRUE),
                 quantile(train_data$time[train_data$status != 0], 0.75, na.rm = TRUE))

# --- Tests for CRModel_RF ---

test_that("CRModel_RF fits basic model with correct interface", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(data = train_data, expvars = expvars,
                        timevar = "time", eventvar = "status",
                        event_codes = "1", ntree = 50)

  expect_error(
    Predict_CRModel_RF(modelout = model_rf, newdata = test_data, event_of_interest = "999"),
    "RF models can only predict"
  )

  # Check output structure
  expect_type(model_rf, "list")
  expect_setequal(names(model_rf),
                  c("rf_model", "times", "varprof", "model_type",
                    "expvars", "timevar", "eventvar", "event_codes",
                    "event_code_numeric", "time_range"))

  # Check model type and class
  expect_equal(model_rf$model_type, "cr_rf")
  expect_s3_class(model_rf$rf_model, "rfsrc")
  expect_s3_class(model_rf, "ml4t2e_cr_rf")

  # Check basic properties
  expect_equal(model_rf$expvars, expvars)
  expect_equal(model_rf$timevar, "time")
  expect_equal(model_rf$eventvar, "status")
  expect_equal(model_rf$event_codes, "1")
  expect_equal(model_rf$event_code_numeric, 1)
  expect_true(is.numeric(model_rf$time_range))
  expect_length(model_rf$time_range, 2)
})

test_that("CRModel_RF handles different event codes", {
  skip_if_not_installed("randomForestSRC")

  # Test with event code = 2
  model_rf_2 <- CRModel_RF(data = train_data, expvars = expvars,
                          timevar = "time", eventvar = "status",
                          event_codes = "2", ntree = 50)

  expect_equal(model_rf_2$event_codes, "2")
  expect_s3_class(model_rf_2$rf_model, "rfsrc")
})

test_that("CRModel_RF handles additional parameters", {
  skip_if_not_installed("randomForestSRC")

  # Test with different parameters
  model_rf_params <- CRModel_RF(data = train_data, expvars = expvars,
                               timevar = "time", eventvar = "status",
                               event_codes = "1", ntree = 30, samplesize = 40, nsplit = 3)

  expect_s3_class(model_rf_params$rf_model, "rfsrc")
  expect_equal(model_rf_params$rf_model$ntree, 30)
})

test_that("CRModel_RF validates inputs", {
  skip_if_not_installed("randomForestSRC")

  # Test missing data
  expect_error(CRModel_RF(expvars = expvars, timevar = "time", eventvar = "status"),
               "argument \"data\" is missing, with no default")

  # Test invalid expvars
  expect_error(CRModel_RF(data = train_data, expvars = character(0),
                         timevar = "time", eventvar = "status"),
               "'expvars' must be a non-empty character vector")

  # Test missing timevar
  expect_error(CRModel_RF(data = train_data, expvars = expvars,
                         timevar = "nonexistent", eventvar = "status"),
               "'timevar' not found in data")

  # Test missing eventvar
  expect_error(CRModel_RF(data = train_data, expvars = expvars,
                         timevar = "time", eventvar = "nonexistent"),
               "'eventvar' not found in data")

  # Test missing expvars in data
  expect_error(CRModel_RF(data = train_data, expvars = c("x1", "nonexistent"),
                         timevar = "time", eventvar = "status"),
               "expvars not found in data")
})

# --- Tests for Predict_CRModel_RF ---

test_that("Predict_CRModel_RF returns predictions in correct format", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(data = train_data, expvars = expvars,
                        timevar = "time", eventvar = "status",
                        event_codes = "1", ntree = 50)

  predictions <- Predict_CRModel_RF(modelout = model_rf, newdata = test_data, event_of_interest = "1")

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

test_that("Predict_CRModel_RF handles custom time points", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(data = train_data, expvars = expvars,
                        timevar = "time", eventvar = "status",
                        event_codes = "1", ntree = 50)

  predictions <- Predict_CRModel_RF(modelout = model_rf, newdata = test_data,
                                   newtimes = time_points, event_of_interest = "1")

  # Check dimensions match requested time points
  expect_equal(length(predictions$Times), length(time_points))
  expect_equal(nrow(predictions$CIFs), length(time_points))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))

  # Check that times match requested times
  expect_equal(as.numeric(predictions$Times), as.numeric(sort(time_points)))
})

test_that("Predict_CRModel_RF validates inputs", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(data = train_data, expvars = expvars,
                        timevar = "time", eventvar = "status",
                        event_codes = "1", ntree = 50)

  expect_error(
    Predict_CRModel_RF(modelout = model_rf, newdata = test_data, event_of_interest = "999"),
    "RF models can only predict"
  )

  # Test invalid modelout
  expect_error(Predict_CRModel_RF(modelout = list(), newdata = test_data, event_of_interest = "1"),
               "'modelout' must be output from CRModel_RF")

  # Test missing newdata
  expect_error(Predict_CRModel_RF(modelout = model_rf),
               "'newdata' is missing")

  # Test missing variables in newdata
  test_data_missing <- test_data[, -which(names(test_data) == "x1")]
  expect_error(Predict_CRModel_RF(modelout = model_rf, newdata = test_data_missing, event_of_interest = "1"),
               "variables missing in newdata")

  # Test invalid newtimes
  expect_error(Predict_CRModel_RF(modelout = model_rf, newdata = test_data,
                                 newtimes = c(-1, 1), event_of_interest = "1"),
               "'newtimes' must be a numeric vector of non-negative values")
})
