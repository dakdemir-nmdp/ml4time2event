library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
# library(dbarts) # Required for the functions being tested - load if needed, or skip

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_bart.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing cr_bart functions")

# --- Test Data Setup ---
# Reusing the setup from previous CR tests
set.seed(123)
n_obs <- 50
cr_data <- data.frame( # Replaced data.table()
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)), # BART might handle factors
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)

# Define formula and parameters
cr_formula <- Surv(time, status) ~ x1 + x2 + x3
train_indices <- 1:40
test_indices <- 41:50
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]

# Time points for prediction
time_points <- c(quantile(train_data$time[train_data$status != 0], 0.25),
                 median(train_data$time[train_data$status != 0]),
                 quantile(train_data$time[train_data$status != 0], 0.75))


# --- Tests for CRModel_BART ---

test_that("CRModel_BART runs and returns a model object", {
  skip_if_not_installed("dbarts") # Or relevant BART package
  # Need to know the expected output structure. Does it fit one model per cause?
  # Assuming it fits one model per cause, similar to Cox/RuleFit wrappers.
  # Let's assume it requires a failcode argument.

  model_bart_list <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 1, ntree = 10) # Assume failcode needed

  # Check output type - what class does the BART package return?
  # dbarts might return a list or a custom object. Assuming a list for now.
  expect_type(model_bart_list, "list")
  expect_length(model_bart_list, 1) # Assuming one model per failcode run
  # Check for expected elements within the BART model object (e.g., 'fit', 'call')
  expect_true(!is.null(model_bart_list[[1]]$fit)) # Example check for dbarts object structure

})

test_that("CRModel_BART requires formula, data, and failcode", {
  skip_if_not_installed("dbarts")
  expect_error(CRModel_BART(formula = cr_formula, data = train_data), "argument \"failcode\" is missing")
  expect_error(CRModel_BART(formula = cr_formula, failcode = 1), "argument \"data\" is missing")
  expect_error(CRModel_BART(data = train_data, failcode = 1), "argument \"formula\" is missing")
})

test_that("CRModel_BART handles different failcode", {
  skip_if_not_installed("dbarts")
  model_bart_c1 <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 1)
  model_bart_c2 <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 2)
  expect_type(model_bart_c1, "list")
  expect_type(model_bart_c2, "list")
  # Check if underlying model fits differ (they should)
  # This comparison depends heavily on the BART object structure
  expect_false(identical(model_bart_c1[[1]]$fit$train.data, model_bart_c2[[1]]$fit$train.data)) # Example check
})


# --- Tests for Predict_CRModel_BART ---

test_that("Predict_CRModel_BART returns predictions in correct format", {
  skip_if_not_installed("dbarts")

  # Fit the model first (for cause 1)
  model_bart_list <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 1, ntree = 10)

  # Get predictions
  # How does prediction work? Does the BART predict method handle competing risks CIF?
  # Assuming the wrapper function handles CIF calculation.
  predictions <- Predict_CRModel_BART(model = model_bart_list, data = test_data, times = time_points)

  # Check output structure
  expect_type(predictions, "list")
  expect_length(predictions, 1) # Predicts for the failcode model was built for
  expect_true(is.matrix(predictions[[1]])) # Expecting matrix for CIF of event 1

  # Check dimensions
  # Rows = number of test observations, Cols = number of time points
  expect_equal(nrow(predictions[[1]]), nrow(test_data))
  expect_equal(ncol(predictions[[1]]), length(time_points))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions[[1]] >= 0 & predictions[[1]] <= 1, na.rm = TRUE))
})

test_that("Predict_CRModel_BART handles single time point", {
  skip_if_not_installed("dbarts")

  model_bart_list <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 1, ntree = 10)
  single_time <- median(train_data$time[train_data$status != 0])
  predictions <- Predict_CRModel_BART(model = model_bart_list, data = test_data, times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]]))
  expect_equal(ncol(predictions[[1]]), 1) # Should have 1 column
  expect_equal(nrow(predictions[[1]]), nrow(test_data))
})

test_that("Predict_CRModel_BART requires model, data, and times", {
  skip_if_not_installed("dbarts")
  model_bart_list <- CRModel_BART(formula = cr_formula, data = train_data, failcode = 1, ntree = 10)
  expect_error(Predict_CRModel_BART(data = test_data, times = time_points), "argument \"model\" is missing")
  expect_error(Predict_CRModel_BART(model = model_bart_list, times = time_points), "argument \"data\" is missing")
  expect_error(Predict_CRModel_BART(model = model_bart_list, data = test_data), "argument \"times\" is missing")
})
