library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # Required for Surv and coxph
library(cmprsk)   # Required for predict.coxphcr (if used internally) or CIF calculation logic

# Assuming the functions are available in the environment
source(here("R/models/cr_cox.R"))

context("Testing cr_cox functions")

# --- Test Data Setup ---
# Reusing the setup from test_cr_random_forest.R for consistency
set.seed(123)
n_obs <- 50
cr_data <- data.frame( # Replaced data.table()
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)
# Create a Surv object suitable for coxph (status needs modification for cause-specific Cox)
# For cause-specific Cox, we typically model each cause separately, treating others as censored.

# Data for cause 1 model: status = 1 if event 1, 0 otherwise
train_data_c1 <- cr_data # Removed copy()
train_data_c1$status_c1 <- ifelse(train_data_c1$status == 1, 1, 0) # Base R assignment
cr_formula_c1 <- Surv(time, status_c1) ~ x1 + x2 + x3

# Data for cause 2 model: status = 1 if event 2, 0 otherwise
train_data_c2 <- cr_data # Removed copy()
train_data_c2$status_c2 <- ifelse(train_data_c2$status == 2, 1, 0) # Base R assignment
cr_formula_c2 <- Surv(time, status_c2) ~ x1 + x2 + x3

# Use original data split
train_indices <- 1:40
test_indices <- 41:50
train_data_orig <- cr_data[train_indices, ]
test_data_orig <- cr_data[test_indices, ]

# Time points for prediction
time_points <- c(quantile(train_data_orig$time[train_data_orig$status != 0], 0.25),
                 median(train_data_orig$time[train_data_orig$status != 0]),
                 quantile(train_data_orig$time[train_data_orig$status != 0], 0.75))

# Original formula for CRModel_Cox (assuming it handles the multi-state internally)
cr_formula_orig <- Surv(time, status) ~ x1 + x2 + x3


# --- Tests for CRModel_Cox ---

test_that("CRModel_Cox runs and returns a list of coxph models", {
  skip_if_not_installed("survival")
  # This test assumes CRModel_Cox fits separate cause-specific models
  # If it uses another approach (like cmprsk::crr), the expectation needs change

  # Need to understand the expected output structure of CRModel_Cox
  # Assuming it returns a list, one element per cause
  model_cox_list <- CRModel_Cox(formula = cr_formula_orig, data = train_data_orig)

  expect_type(model_cox_list, "list")
  # Assuming models for event 1 and 2 are returned
  expect_length(model_cox_list, 2)
  expect_s3_class(model_cox_list[[1]], "coxph")
  expect_s3_class(model_cox_list[[2]], "coxph")

  # Check if models used the correct data subset/status internally (hard to verify directly)
})

test_that("CRModel_Cox requires formula and data", {
  skip_if_not_installed("survival")
  expect_error(CRModel_Cox(formula = cr_formula_orig), "argument \"data\" is missing")
  expect_error(CRModel_Cox(data = train_data_orig), "argument \"formula\" is missing")
})


# --- Tests for Predict_CRModel_Cox ---

test_that("Predict_CRModel_Cox returns predictions in correct format", {
  skip_if_not_installed("survival")
  skip_if_not_installed("cmprsk") # Often needed for CIF calculation

  # Fit the model first
  model_cox_list <- CRModel_Cox(formula = cr_formula_orig, data = train_data_orig)

  # Get predictions
  predictions <- Predict_CRModel_Cox(model = model_cox_list, data = test_data_orig, times = time_points)

  # Check output structure
  expect_type(predictions, "list")
  expect_length(predictions, 2) # One element per cause
  expect_true(is.matrix(predictions[[1]])) # Expecting matrix for CIF of event 1
  expect_true(is.matrix(predictions[[2]])) # Expecting matrix for CIF of event 2

  # Check dimensions
  # Rows = number of test observations, Cols = number of time points
  expect_equal(nrow(predictions[[1]]), nrow(test_data_orig))
  expect_equal(ncol(predictions[[1]]), length(time_points))
  expect_equal(nrow(predictions[[2]]), nrow(test_data_orig))
  expect_equal(ncol(predictions[[2]]), length(time_points))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions[[1]] >= 0 & predictions[[1]] <= 1, na.rm = TRUE))
  expect_true(all(predictions[[2]] >= 0 & predictions[[2]] <= 1, na.rm = TRUE))

  # Check sum of CIFs <= 1 for each time point and observation
  expect_true(all(predictions[[1]] + predictions[[2]] <= 1.0001, na.rm = TRUE)) # Allow small tolerance
})

test_that("Predict_CRModel_Cox handles single time point", {
  skip_if_not_installed("survival")
  skip_if_not_installed("cmprsk")

  model_cox_list <- CRModel_Cox(formula = cr_formula_orig, data = train_data_orig)
  single_time <- median(train_data_orig$time[train_data_orig$status != 0])
  predictions <- Predict_CRModel_Cox(model = model_cox_list, data = test_data_orig, times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]]))
  expect_equal(ncol(predictions[[1]]), 1) # Should have 1 column
  expect_equal(nrow(predictions[[1]]), nrow(test_data_orig))
})

test_that("Predict_CRModel_Cox requires model, data, and times", {
  skip_if_not_installed("survival")
  model_cox_list <- CRModel_Cox(formula = cr_formula_orig, data = train_data_orig)
  expect_error(Predict_CRModel_Cox(data = test_data_orig, times = time_points), "argument \"model\" is missing")
  expect_error(Predict_CRModel_Cox(model = model_cox_list, times = time_points), "argument \"data\" is missing")
  expect_error(Predict_CRModel_Cox(model = model_cox_list, data = test_data_orig), "argument \"times\" is missing")
})
