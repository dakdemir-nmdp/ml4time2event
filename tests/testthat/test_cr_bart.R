library(testthat)
library(BART) # Required for the functions being tested

context("Testing cr_bart functions")

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

# --- Tests for CRModel_BART ---

test_that("CRModel_BART fits basic model with correct interface", {
  skip_if_not_installed("BART")

  model_bart <- CRModel_BART(data = train_data, expvars = expvars,
                            timevar = "time", eventvar = "status",
                            ntree = 20, ndpost = 50, nskip = 25)

  # Check output structure
  expect_type(model_bart, "list")
  expect_named(model_bart, c("bart_model", "times", "varprof", "model_type",
                            "expvars", "timevar", "eventvar", "failcode",
                            "time_range", "x_train", "times_train", "delta_train"))

  # Check model type and class
  expect_equal(model_bart$model_type, "cr_bart")
  expect_s3_class(model_bart$bart_model, "criskbart")
  expect_s3_class(model_bart, "ml4t2e_cr_bart")

  # Check basic properties
  expect_equal(model_bart$expvars, expvars)
  expect_equal(model_bart$timevar, "time")
  expect_equal(model_bart$eventvar, "status")
  expect_equal(model_bart$failcode, 1)
  expect_true(is.numeric(model_bart$time_range))
  expect_length(model_bart$time_range, 2)
})

test_that("CRModel_BART handles different failcode values", {
  skip_if_not_installed("BART")

  # Test with failcode = 2
  model_bart_2 <- CRModel_BART(data = train_data, expvars = expvars,
                              timevar = "time", eventvar = "status",
                              failcode = 2, ntree = 20, ndpost = 50, nskip = 25)

  expect_equal(model_bart_2$failcode, 2)
  expect_s3_class(model_bart_2$bart_model, "criskbart")
})

test_that("CRModel_BART handles additional parameters", {
  skip_if_not_installed("BART")

  # Test with different parameters
  model_bart_params <- CRModel_BART(data = train_data, expvars = expvars,
                                   timevar = "time", eventvar = "status",
                                   K = 5, ntree = 30, ndpost = 50, nskip = 25)

  expect_s3_class(model_bart_params$bart_model, "criskbart")
})

test_that("CRModel_BART validates inputs", {
  skip_if_not_installed("BART")

  # Test missing data
  expect_error(CRModel_BART(expvars = expvars, timevar = "time", eventvar = "status"),
               "argument \"data\" is missing, with no default")

  # Test invalid expvars
  expect_error(CRModel_BART(data = train_data, expvars = character(0),
                           timevar = "time", eventvar = "status"),
               "'expvars' must be a non-empty character vector")

  # Test missing timevar
  expect_error(CRModel_BART(data = train_data, expvars = expvars,
                           timevar = "nonexistent", eventvar = "status"),
               "'timevar' not found in data")

  # Test missing eventvar
  expect_error(CRModel_BART(data = train_data, expvars = expvars,
                           timevar = "time", eventvar = "nonexistent"),
               "'eventvar' not found in data")

  # Test missing expvars in data
  expect_error(CRModel_BART(data = train_data, expvars = c("x1", "nonexistent"),
                           timevar = "time", eventvar = "status"),
               "expvars not found in data")
})

# --- Tests for Predict_CRModel_BART ---

test_that("Predict_CRModel_BART returns predictions in correct format", {
  skip_if_not_installed("BART")

  model_bart <- CRModel_BART(data = train_data, expvars = expvars,
                            timevar = "time", eventvar = "status",
                            ntree = 20, ndpost = 50, nskip = 25)

  predictions <- Predict_CRModel_BART(modelout = model_bart, newdata = test_data)

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

test_that("Predict_CRModel_BART handles custom time points", {
  skip_if_not_installed("BART")

  model_bart <- CRModel_BART(data = train_data, expvars = expvars,
                            timevar = "time", eventvar = "status",
                            ntree = 20, ndpost = 50, nskip = 25)

  predictions <- Predict_CRModel_BART(modelout = model_bart, newdata = test_data,
                                     newtimes = time_points)

  # Check dimensions match requested time points
  expect_equal(length(predictions$Times), length(time_points))
  expect_equal(nrow(predictions$CIFs), length(time_points))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))

  # Check that times match requested times
  expect_equal(as.numeric(predictions$Times), as.numeric(sort(time_points)))
})

test_that("Predict_CRModel_BART validates inputs", {
  skip_if_not_installed("BART")

  model_bart <- CRModel_BART(data = train_data, expvars = expvars,
                            timevar = "time", eventvar = "status",
                            ntree = 20, ndpost = 50, nskip = 25)

  # Test invalid modelout
  expect_error(Predict_CRModel_BART(modelout = list(), newdata = test_data),
               "'modelout' must be output from CRModel_BART")

  # Test missing newdata
  expect_error(Predict_CRModel_BART(modelout = model_bart),
               "argument \"newdata\" is missing, with no default")

  # Test missing variables in newdata
  test_data_missing <- test_data[, -which(names(test_data) == "x1")]
  expect_error(Predict_CRModel_BART(modelout = model_bart, newdata = test_data_missing),
               "variables missing in newdata")

  # Test invalid newtimes
  expect_error(Predict_CRModel_BART(modelout = model_bart, newdata = test_data,
                                   newtimes = c(-1, 1)),
               "'newtimes' must be a numeric vector of non-negative values")
})
