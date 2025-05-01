library(testthat)
library(here)
# library(data.table) # Removed
library(survival) # For Surv object
library(cmprsk)   # Required for crr and predict.crr

# Assuming the functions are available in the environment
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_fine_gray.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")

context("Testing cr_fine_gray functions")

# --- Test Data Setup ---
# Reusing the setup from previous CR tests
set.seed(123)
n_obs <- 50
cr_data <- data.frame( # Replaced data.table()
  ftime = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  fstatus = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = sample(c("A", "B"), n_obs, replace = TRUE), # Keep as character for crr
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)
# crr requires numeric matrix for covariates, no factors allowed directly in formula
# Need to convert factor x2 to dummy variables
# Replace data.table assignment with base R
cr_data$x2A <- ifelse(cr_data$x2 == "A", 1, 0)
cr_data$x2B <- ifelse(cr_data$x2 == "B", 1, 0)

# Define formula and parameters for crr (using dummy vars)
# Note: crr uses cov1, cov2 arguments, not formula
# The CRModel_FineGray function likely handles the formula conversion internally
cr_formula_orig <- Surv(ftime, fstatus) ~ x1 + x2 + x3 # Original formula for wrapper

# Prepare data for crr (numeric matrix for covariates)
# Assuming CRModel_FineGray handles this conversion
cov_vars <- c("x1", "x2A", "x2B", "x3") # Example, depends on internal handling
# cov_matrix <- as.matrix(cr_data[, ..cov_vars]) # Example matrix

# Use original data split
train_indices <- 1:40
test_indices <- 41:50
train_data_orig <- cr_data[train_indices, ]
test_data_orig <- cr_data[test_indices, ]

# Time points for prediction
time_points <- c(quantile(train_data_orig$ftime[train_data_orig$fstatus != 0], 0.25),
                 median(train_data_orig$ftime[train_data_orig$fstatus != 0]),
                 quantile(train_data_orig$ftime[train_data_orig$fstatus != 0], 0.75))


# --- Tests for CRModel_FineGray ---

test_that("CRModel_FineGray runs and returns a crr model object", {
  skip_if_not_installed("cmprsk")

  # Assuming the function fits a model for cause 1 by default or specified
  # Need to know how the function selects the cause of interest
  # Let's assume it fits for cause 1 if not specified
  model_fg <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 1, iter.max = 10) # Specify failcode

  # Check output type
  expect_s3_class(model_fg, "crr")

  # Check basic model properties (if available in crr object)
  # expect_equal(model_fg$n, nrow(train_data_orig)) # crr object structure might differ
})

test_that("CRModel_FineGray requires formula, data, and failcode", {
  skip_if_not_installed("cmprsk")
  expect_error(CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig), "argument \"failcode\" is missing")
  expect_error(CRModel_FineGray(formula = cr_formula_orig, failcode = 1), "argument \"data\" is missing")
  expect_error(CRModel_FineGray(data = train_data_orig, failcode = 1), "argument \"formula\" is missing")
})

test_that("CRModel_FineGray handles different failcode", {
  skip_if_not_installed("cmprsk")
  # Fit for cause 2
  model_fg_c2 <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 2)
  expect_s3_class(model_fg_c2, "crr")
  # Check if coefficients differ from cause 1 model (they should)
  model_fg_c1 <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 1)
  expect_false(identical(coef(model_fg_c1), coef(model_fg_c2)))
})


# --- Tests for Predict_CRModel_FG ---

test_that("Predict_CRModel_FG returns predictions in correct format", {
  skip_if_not_installed("cmprsk")

  # Fit the model first (for cause 1)
  model_fg <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 1)

  # Get predictions
  # Predict_CRModel_FG likely calls predict.crr internally
  # Need to ensure test_data covariates are formatted correctly (numeric matrix)
  # Assuming Predict_CRModel_FG handles this conversion based on the formula
  predictions <- Predict_CRModel_FG(model = model_fg, data = test_data_orig, times = time_points)

  # Check output structure (predict.crr output is slightly different)
  # It returns a matrix where rows are time points, columns are subjects
  # The wrapper function should reformat this to subjects x time points
  expect_type(predictions, "list") # Should still return a list per cause
  expect_length(predictions, 1) # Only predicts for the failcode the model was built for
  expect_true(is.matrix(predictions[[1]])) # Expecting matrix for CIF of event 1

  # Check dimensions
  # Rows = number of test observations, Cols = number of time points
  expect_equal(nrow(predictions[[1]]), nrow(test_data_orig))
  expect_equal(ncol(predictions[[1]]), length(time_points))

  # Check values are probabilities (between 0 and 1)
  expect_true(all(predictions[[1]] >= 0 & predictions[[1]] <= 1, na.rm = TRUE))
})

test_that("Predict_CRModel_FG handles single time point", {
  skip_if_not_installed("cmprsk")

  model_fg <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 1)
  single_time <- median(train_data_orig$ftime[train_data_orig$fstatus != 0])
  predictions <- Predict_CRModel_FG(model = model_fg, data = test_data_orig, times = single_time)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions[[1]]))
  expect_equal(ncol(predictions[[1]]), 1) # Should have 1 column
  expect_equal(nrow(predictions[[1]]), nrow(test_data_orig))
})

test_that("Predict_CRModel_FG requires model, data, and times", {
  skip_if_not_installed("cmprsk")
  model_fg <- CRModel_FineGray(formula = cr_formula_orig, data = train_data_orig, failcode = 1, iter.max = 10)
  expect_error(Predict_CRModel_FG(data = test_data_orig, times = time_points), "argument \"model\" is missing")
  expect_error(Predict_CRModel_FG(model = model_fg, times = time_points), "argument \"data\" is missing")
  expect_error(Predict_CRModel_FG(model = model_fg, data = test_data_orig), "argument \"times\" is missing")
})
