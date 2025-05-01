library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
# library(pre) # Required for the functions being tested - load if needed, or skip
# library(cmprsk) # May be needed for prediction logic if not handled by pre directly

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_rulefit.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing cr_rulefit functions")

# --- Test Data Setup ---
# Reusing the setup from previous CR tests
set.seed(123)
n_obs <- 50
cr_data <- data.frame( # Replaced data.table()
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)), # pre handles factors
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


# --- Tests for CRModel_rulefit ---

test_that("CRModel_rulefit runs and returns a model object", {
  skip_if_not_installed("pre")
  # Need to know the expected output structure. Does it fit one model per cause?
  # Assuming it fits one model per cause, similar to Cox.
  # Let's assume it requires a failcode argument like Fine-Gray.

  model_rulefit_list <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 1, maxdepth = 2) # Assume failcode needed

  # Check output type - what class does pre return for competing risks?
  # It might return a list of 'pre' objects or a custom class.
  # Let's assume a list of 'pre' objects for now.
  expect_type(model_rulefit_list, "list")
  expect_length(model_rulefit_list, 1) # Assuming one model per failcode run
  expect_s3_class(model_rulefit_list[[1]], "pre")

  # Check basic model properties if possible
  # expect_equal(model_rulefit$call[[1]], quote(pre)) # Check call if available
})

test_that("CRModel_rulefit requires formula, data, and failcode", {
  skip_if_not_installed("pre")
  expect_error(CRModel_rulefit(formula = cr_formula, data = train_data), "argument \"failcode\" is missing")
  expect_error(CRModel_rulefit(formula = cr_formula, failcode = 1), "argument \"data\" is missing")
  expect_error(CRModel_rulefit(data = train_data, failcode = 1), "argument \"formula\" is missing")
})

test_that("CRModel_rulefit handles different failcode", {
  skip_if_not_installed("pre")
  model_rf_c1 <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 1)
  model_rf_c2 <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 2)
  expect_s3_class(model_rf_c1[[1]], "pre")
  expect_s3_class(model_rf_c2[[1]], "pre")
  # Check if underlying model coefficients/rules differ (they should)
  # This comparison might be complex depending on the 'pre' object structure
  expect_false(identical(model_rf_c1[[1]]$glmnet.fit, model_rf_c2[[1]]$glmnet.fit)) # Example check
})


# --- Tests for Predict_CRModel_rulefit ---

test_that("Predict_CRModel_rulefit returns predictions in correct format", {
  skip_if_not_installed("pre")
  # skip_if_not_installed("cmprsk") # May be needed

  # Fit the model first (for cause 1)
  model_rulefit_list <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 1)

  # Get predictions
  # How does prediction work? Does pre::predict handle competing risks CIF directly?
  # Assuming the wrapper function handles CIF calculation.
  predictions <- Predict_CRModel_rulefit(model = model_rulefit_list, data = test_data, times = time_points)

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

test_that("Predict_CRModel_rulefit handles single time point", {
  skip_if_not_installed("pre")

  model_rulefit_list <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 1)
  single_time <- median(train_data$time[train_data$status != 0])
  predictions <- Predict_CRModel_rulefit(model = model_rulefit_list, data = test_data, times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]]))
  expect_equal(ncol(predictions[[1]]), 1) # Should have 1 column
  expect_equal(nrow(predictions[[1]]), nrow(test_data))
})

test_that("Predict_CRModel_rulefit requires model, data, and times", {
  skip_if_not_installed("pre")
  model_rulefit_list <- CRModel_rulefit(formula = cr_formula, data = train_data, failcode = 1, maxdepth = 2)
  expect_error(Predict_CRModel_rulefit(data = test_data, times = time_points), "argument \"model\" is missing")
  expect_error(Predict_CRModel_rulefit(model = model_rulefit_list, times = time_points), "argument \"data\" is missing")
  expect_error(Predict_CRModel_rulefit(model = model_rulefit_list, data = test_data), "argument \"times\" is missing")
})
