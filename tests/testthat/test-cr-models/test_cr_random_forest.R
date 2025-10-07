library(testthat)
library(here)
# library(data.table) # Removed
library(randomForestSRC) # Required for the functions being tested
library(survival) # For Surv object

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_random_forest.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing cr_random_forest functions")

# --- Test Data Setup ---
# Create a simple competing risks dataset
set.seed(123)
n_obs <- 50
cr_data <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)
# Ensure status 0 has max time for simplicity in this small dataset? Not strictly necessary.
# cr_data[status == 0, time := max(time) + 1]

# Define formula and parameters
cr_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices <- 1:40
test_indices <- 41:50
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]
time_points <- c(quantile(train_data$time[train_data$status != 0], 0.25),
                 median(train_data$time[train_data$status != 0]),
                 quantile(train_data$time[train_data$status != 0], 0.75))


# --- Tests for CRModel_RF ---

test_that("CRModel_RF runs and returns a model object", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(formula = cr_formula, data = train_data)

  # Check output type
  expect_s3_class(model_rf, "rfsrc")
  expect_true(!is.null(model_rf$cif)) # Check if CIF results are present

  # Check basic model properties
  expect_equal(model_rf$call[[1]], quote(rfsrc))
  expect_equal(model_rf$n, nrow(train_data))
})

test_that("CRModel_RF handles additional parameters", {
  skip_if_not_installed("randomForestSRC")

  # Test with different ntree and nodesize
  model_rf_params <- CRModel_RF(formula = cr_formula, data = train_data, ntree = 10, nodesize = 10)

  expect_s3_class(model_rf_params, "rfsrc")
  expect_equal(model_rf_params$ntree, 10)
  # Note: nodesize might not be directly stored with that name, check rfsrc object structure if needed
  # expect_equal(model_rf_params$nodesize, 10) # This might fail depending on rfsrc object
})

test_that("CRModel_RF requires formula and data", {
  skip_if_not_installed("randomForestSRC")
  expect_error(CRModel_RF(formula = cr_formula), "argument \"data\" is missing")
  expect_error(CRModel_RF(data = train_data), "argument \"formula\" is missing")
})


# --- Tests for Predict_CRModel_RF ---

test_that("Predict_CRModel_RF returns predictions in correct format", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(formula = cr_formula, data = train_data, ntree = 10) # Faster model
  predictions <- Predict_CRModel_RF(model = model_rf, data = test_data, times = time_points)

  # Check output structure
  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]])) # Expecting matrix for CIF of event 1
  expect_true(is.matrix(predictions[[2]])) # Expecting matrix for CIF of event 2

  # Check dimensions
  # Rows = number of test observations, Cols = number of time points
  expect_equal(nrow(predictions[[1]]), nrow(test_data))
  expect_equal(ncol(predictions[[1]]), length(time_points))
  expect_equal(nrow(predictions[[2]]), nrow(test_data))
  expect_equal(ncol(predictions[[2]]), length(time_points))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions[[1]] >= 0 & predictions[[1]] <= 1, na.rm = TRUE))
  expect_true(all(predictions[[2]] >= 0 & predictions[[2]] <= 1, na.rm = TRUE))

  # Check sum of CIFs <= 1 for each time point and observation
  expect_true(all(predictions[[1]] + predictions[[2]] <= 1.0001, na.rm = TRUE)) # Allow small tolerance
})

test_that("Predict_CRModel_RF handles single time point", {
  skip_if_not_installed("randomForestSRC")

  model_rf <- CRModel_RF(formula = cr_formula, data = train_data, ntree = 10)
  single_time <- median(train_data$time[train_data$status != 0])
  predictions <- Predict_CRModel_RF(model = model_rf, data = test_data, times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]]))
  expect_equal(ncol(predictions[[1]]), 1) # Should have 1 column
  expect_equal(nrow(predictions[[1]]), nrow(test_data))
})

test_that("Predict_CRModel_RF requires model, data, and times", {
  skip_if_not_installed("randomForestSRC")
  model_rf <- CRModel_RF(formula = cr_formula, data = train_data, ntree = 10)
  expect_error(Predict_CRModel_RF(data = test_data, times = time_points), "argument \"model\" is missing")
  expect_error(Predict_CRModel_RF(model = model_rf, times = time_points), "argument \"data\" is missing")
  expect_error(Predict_CRModel_RF(model = model_rf, data = test_data), "argument \"times\" is missing")
})
