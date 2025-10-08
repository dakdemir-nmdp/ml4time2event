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

test_that("CRModel_FineGray runs and returns expected model structure", {
  skip_if_not_installed("fastcmprsk")

  # Fit model using the correct interface
  model_fg <- CRModel_FineGray(data = train_data_orig,
                              expvars = c("x1", "x2A", "x2B", "x3"),
                              timevar = "ftime",
                              eventvar = "fstatus",
                              failcode = 1)

  # Check output structure
  expect_type(model_fg, "list")
  expect_s3_class(model_fg, "ml4t2e_cr_finegray")
  expect_equal(model_fg$model_type, "fine_gray")
  expect_equal(model_fg$failcode, 1)
  expect_true("fg_model" %in% names(model_fg))
  expect_true("time_range" %in% names(model_fg))
  expect_true("varprof" %in% names(model_fg))
})

test_that("CRModel_FineGray requires data, expvars, timevar, eventvar", {
  skip_if_not_installed("fastcmprsk")
  expect_error(CRModel_FineGray(expvars = c("x1", "x2A", "x2B", "x3"), timevar = "ftime", eventvar = "fstatus"), "argument \"data\" is missing")
  expect_error(CRModel_FineGray(data = train_data_orig, timevar = "ftime", eventvar = "fstatus"), "argument \"expvars\" is missing")
  expect_error(CRModel_FineGray(data = train_data_orig, expvars = c("x1", "x2A", "x2B", "x3"), eventvar = "fstatus"), "argument \"timevar\" is missing")
  expect_error(CRModel_FineGray(data = train_data_orig, expvars = c("x1", "x2A", "x2B", "x3"), timevar = "ftime"), "argument \"eventvar\" is missing")
})

test_that("CRModel_FineGray handles different failcode", {
  skip_if_not_installed("fastcmprsk")
  # Fit for cause 2
  model_fg_c2 <- CRModel_FineGray(data = train_data_orig,
                                 expvars = c("x1", "x2A", "x2B", "x3"),
                                 timevar = "ftime",
                                 eventvar = "fstatus",
                                 failcode = 2)
  expect_s3_class(model_fg_c2, "ml4t2e_cr_finegray")
  # Check if models differ for different failcodes
  model_fg_c1 <- CRModel_FineGray(data = train_data_orig,
                                 expvars = c("x1", "x2A", "x2B", "x3"),
                                 timevar = "ftime",
                                 eventvar = "fstatus",
                                 failcode = 1)
  expect_equal(model_fg_c1$failcode, 1)
  expect_equal(model_fg_c2$failcode, 2)
})


# --- Tests for Predict_CRModel_FineGray ---

test_that("Predict_CRModel_FineGray returns predictions in correct format", {
  skip_if_not_installed("fastcmprsk")

  # Fit the model first (for cause 1)
  model_fg <- CRModel_FineGray(data = train_data_orig,
                              expvars = c("x1", "x2A", "x2B", "x3"),
                              timevar = "ftime",
                              eventvar = "fstatus",
                              failcode = 1)

  # Get predictions
  predictions <- Predict_CRModel_FineGray(modelout = model_fg, newdata = test_data_orig)

  # Check output structure
  expect_type(predictions, "list")
  expect_true("CIFs" %in% names(predictions))
  expect_true("Times" %in% names(predictions))
  expect_true(is.matrix(predictions$CIFs))
  expect_true(is.vector(predictions$Times))

  # CIFs should be [times, observations]
  expect_equal(ncol(predictions$CIFs), nrow(test_data_orig))
  expect_equal(length(predictions$Times), nrow(predictions$CIFs))

  # CIFs should start at 0 and be non-decreasing
  expect_true(all(predictions$CIFs[1, ] == 0))
  expect_true(all(diff(predictions$CIFs[, 1]) >= 0))  # Check monotonicity for first observation
})

test_that("Predict_CRModel_FineGray handles custom time points", {
  skip_if_not_installed("fastcmprsk")

  model_fg <- CRModel_FineGray(data = train_data_orig,
                              expvars = c("x1", "x2A", "x2B", "x3"),
                              timevar = "ftime",
                              eventvar = "fstatus",
                              failcode = 1)
  custom_times <- c(1, 5, 10)
  predictions <- Predict_CRModel_FineGray(modelout = model_fg, newdata = test_data_orig, newtimes = custom_times)

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions$CIFs))
  expect_equal(length(predictions$Times), length(custom_times))
  expect_equal(nrow(predictions$CIFs), length(custom_times))
  expect_equal(ncol(predictions$CIFs), nrow(test_data_orig))
})

test_that("Predict_CRModel_FineGray requires model and data", {
  skip_if_not_installed("fastcmprsk")

  model_fg <- CRModel_FineGray(data = train_data_orig,
                              expvars = c("x1", "x2A", "x2B", "x3"),
                              timevar = "ftime",
                              eventvar = "fstatus",
                              failcode = 1)

  expect_error(Predict_CRModel_FineGray(modelout = model_fg), "argument \"newdata\" is missing")
  expect_error(Predict_CRModel_FineGray(newdata = test_data_orig), "argument \"modelout\" is missing")
})
