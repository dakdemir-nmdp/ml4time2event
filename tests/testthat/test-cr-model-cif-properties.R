library(survival)
library(testthat)

# Load the package to make all CR model functions available
devtools::load_all()

# Helper function to check for non-monotonicity in a CIF matrix
# A valid CIF must be non-decreasing.
# This function returns TRUE if any CIF curve is non-monotonic.
is_cif_non_monotonic <- function(cif_matrix) {
  # Input: A matrix where columns are observations and rows are times
  if (!is.matrix(cif_matrix)) {
    cif_matrix <- as.matrix(cif_matrix)
  }
  
  # Apply check to each column (each observation's CIF)
  non_monotonic_cols <- apply(cif_matrix, 2, function(cif_curve) {
    # Check for NAs first
    if (any(is.na(cif_curve))) return(FALSE) # Can't check monotonicity with NAs
    # Check if any difference is negative (i.e., the CIF decreases)
    any(diff(cif_curve) < -1e-9) # Use a small tolerance
  })
  
  any(non_monotonic_cols)
}

# Helper function to check for near-constancy
# This checks if the CIF jumps at the first step and then stays flat.
is_cif_near_constant <- function(cif_matrix) {
  if (!is.matrix(cif_matrix)) {
    cif_matrix <- as.matrix(cif_matrix)
  }
  
  near_constant_cols <- apply(cif_matrix, 2, function(cif_curve) {
    if (any(is.na(cif_curve)) || length(cif_curve) < 3) return(FALSE)
    
    # Check if the first increment is large and subsequent increments are tiny
    increments <- diff(cif_curve)
    first_increment <- increments[1]
    rest_increments <- increments[-1]
    
    # Heuristic: first increment is > 0.01 and max of the rest is < 1% of the first
    first_increment > 0.01 && all(abs(rest_increments) < first_increment * 0.01)
  })
  
  any(near_constant_cols)
}


# 1. Prepare the PBC dataset for competing risks
data(pbc, package = "survival")
pbc$status[pbc$status == 2] <- 3 # Use 3 for death to avoid confusion with event 2
pbc$status[pbc$status == 1] <- 2 # Use 2 for transplant
pbc$status[pbc$status == 0] <- 1 # Use 1 for censored
pbc <- na.omit(pbc) # Simple case for testing, remove NAs
pbc$status <- as.numeric(pbc$status) # Ensure event variable is numeric

# Define variables
time_var <- "time"
event_var <- "status"
exp_vars <- c("age", "sex", "bili", "chol", "albumin")
event_codes <- c("2", "3") # Transplant and Death

# Split data
set.seed(123)
train_indices <- sample(seq_len(nrow(pbc)), size = 0.7 * nrow(pbc))
train_data <- pbc[train_indices, ]
test_data <- pbc[-train_indices, ]

# Define evaluation times
eval_times <- seq(min(train_data$time), max(train_data$time), length.out = 20)

# --- Test for CRModel_FineGray ---
test_that("CRModel_FineGray produces monotonic CIFs", {
  
  # Train model
  fg_model_out <- CRModel_FineGray(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = "2" # Event of interest: transplant (coded as 2)
  )
  
  # Predict CIFs
  fg_preds <- Predict_CRModel_FineGray(
    modelout = fg_model_out,
    newdata = test_data,
    new_times = eval_times
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(fg_preds$CIFs), 
               "Fine-Gray model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(fg_preds$CIFs),
               "Fine-Gray model produced near-constant CIFs.")
})

# --- Test for CRModel_Cox ---
test_that("CRModel_Cox produces monotonic CIFs", {
  
  # Train model
  cox_model_out <- CRModel_Cox(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = event_codes
  )
  
  # Predict CIFs
  cox_preds <- Predict_CRModel_Cox(
    modelout = cox_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(cox_preds$CIFs), 
               "Cox model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(cox_preds$CIFs),
               "Cox model produced near-constant CIFs.")
})

# --- Test for CRModel_GAM ---
test_that("CRModel_GAM produces monotonic CIFs", {
  
  # Train model
  gam_model_out <- CRModel_GAM(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_of_interest = 2
  )

  # Predict CIFs
  gam_preds <- Predict_CRModel_GAM(
    modelout = gam_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = 2
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(gam_preds$CIFs), 
               "GAM model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(gam_preds$CIFs),
               "GAM model produced near-constant CIFs.")
})

# --- Test for CRModel_RF ---
test_that("CRModel_RF produces monotonic CIFs", {
  
  # Train model
  rf_model_out <- CRModel_RF(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = "2"
  )
  
  # Predict CIFs
  rf_preds <- Predict_CRModel_RF(
    modelout = rf_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(rf_preds$CIFs), 
               "Random Forest model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(rf_preds$CIFs),
               "Random Forest model produced near-constant CIFs.")
})

# --- Test for CRModel_xgboost ---
test_that("CRModel_xgboost produces monotonic CIFs", {
  
  # Train model
  xgb_model_out <- CRModel_xgboost(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = event_codes
  )
  
  # Predict CIFs
  xgb_preds <- Predict_CRModel_xgboost(
    modelout = xgb_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(xgb_preds$CIFs), 
               "XGBoost model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(xgb_preds$CIFs),
               "XGBoost model produced near-constant CIFs.")
})

# --- Test for CRModel_FineGray ---
test_that("CRModel_FineGray produces monotonic CIFs", {
  
  # Train model
  fg_model_out <- CRModel_FineGray(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = "2" # Event of interest: transplant (coded as 2)
  )
  
  # Predict CIFs
  fg_preds <- Predict_CRModel_FineGray(
    modelout = fg_model_out,
    newdata = test_data,
    new_times = eval_times
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(fg_preds$CIFs), 
               "Fine-Gray model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(fg_preds$CIFs),
               "Fine-Gray model produced near-constant CIFs.")
})

# --- Test for CRModel_Cox ---
test_that("CRModel_Cox produces monotonic CIFs", {
  
  # Train model
  cox_model_out <- CRModel_Cox(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = event_codes
  )
  
  # Predict CIFs
  cox_preds <- Predict_CRModel_Cox(
    modelout = cox_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(cox_preds$CIFs), 
               "Cox model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(cox_preds$CIFs),
               "Cox model produced near-constant CIFs.")
})

# --- Test for CRModel_GAM ---
test_that("CRModel_GAM produces monotonic CIFs", {
  
  # Train model
  gam_model_out <- CRModel_GAM(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_of_interest = 2
  )

  # Predict CIFs
  gam_preds <- Predict_CRModel_GAM(
    modelout = gam_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = 2
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(gam_preds$CIFs), 
               "GAM model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(gam_preds$CIFs),
               "GAM model produced near-constant CIFs.")
})

# --- Test for CRModel_RF ---
test_that("CRModel_RF produces monotonic CIFs", {
  
  # Train model
  rf_model_out <- CRModel_RF(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = "2"
  )
  
  # Predict CIFs
  rf_preds <- Predict_CRModel_RF(
    modelout = rf_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(rf_preds$CIFs), 
               "Random Forest model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(rf_preds$CIFs),
               "Random Forest model produced near-constant CIFs.")
})

# --- Test for CRModel_xgboost ---
test_that("CRModel_xgboost produces monotonic CIFs", {
  
  # Train model
  xgb_model_out <- CRModel_xgboost(
    data = train_data,
    timevar = time_var,
    eventvar = event_var,
    expvars = exp_vars,
    event_codes = event_codes
  )
  
  # Predict CIFs
  xgb_preds <- Predict_CRModel_xgboost(
    modelout = xgb_model_out,
    newdata = test_data,
    new_times = eval_times,
    event_of_interest = "2"
  )
  
  # Check for monotonicity
  expect_false(is_cif_non_monotonic(xgb_preds$CIFs), 
               "XGBoost model produced non-monotonic CIFs.")
               
  # Check for near-constancy
  expect_false(is_cif_near_constant(xgb_preds$CIFs),
               "XGBoost model produced near-constant CIFs.")
})

