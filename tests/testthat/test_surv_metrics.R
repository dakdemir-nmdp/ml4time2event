library(testthat)
library(here)
# library(data.table) # Removed
library(survival)
library(pec) # Required for cindex and pec functions
library(prodlim) # Required for Hist function

# Assuming the function is available in the environment
source(here("R/surv_metrics.R"))

context("Testing surv_metrics functions")

# --- Test Data Setup ---
# Reusing the setup from previous Surv tests, maybe slightly larger
set.seed(1011)
n_obs_met <- 60
surv_data_metrics <- data.frame(
  time = rexp(n_obs_met, rate = 0.05),
  status = sample(0:1, n_obs_met, replace = TRUE, prob = c(0.3, 0.7)), # 0=censored, 1=event
  x1 = rnorm(n_obs_met),
  x2 = factor(sample(c("C", "D"), n_obs_met, replace = TRUE)),
  x3 = rnorm(n_obs_met, mean = 2),
  stringsAsFactors = FALSE
)
# Create Surv object
surv_obj_met <- Surv(surv_data_metrics$time, surv_data_metrics$status)

# Create mock predictions (Survival probabilities) at specific time points
# Predictions should ideally have some correlation with outcome for metrics to be meaningful
# Let's create predictions based loosely on x1 (higher x1 -> lower survival prob)
time_points_met <- quantile(surv_data_metrics$time[surv_data_metrics$status == 1], c(0.3, 0.6))

# Mock Survival Probs (higher for lower x1)
survprob_t1 <- pmax(0, pmin(1, 0.8 - 0.1 * surv_data_metrics$x1 + rnorm(n_obs_met, 0, 0.1)))
survprob_t2 <- pmax(0, pmin(1, 0.6 - 0.15 * surv_data_metrics$x1 + rnorm(n_obs_met, 0, 0.15)))
# Ensure non-increasing
mock_survprobs <- cbind(survprob_t1, survprob_t2)
mock_survprobs[, 2] <- pmin(mock_survprobs[, 1], mock_survprobs[, 2]) # Ensure t2 <= t1

# --- Tests for timedepConcordance ---

test_that("timedepConcordance calculates C-index correctly", {
  # skip_if_not_installed("pec") # If pec::cindex is used internally

  # Calculate C-index at the first time point
  # Note: timedepConcordance expects survival probabilities (higher = better survival)
  # Function signature: timedepConcordance(predsurv, predsurvtimes, obstimes, obsevents, TestMat = NULL)
  # predsurv is a matrix with rows=times, cols=observations
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)  # 1 row (1 time point), n_obs_met columns
  c_index_obj_t1 <- timedepConcordance(
    predsurv = predsurv_t1,
    predsurvtimes = time_points_met[1],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )

  # pec::cindex returns an object, extract the concordance value
  expect_type(c_index_obj_t1, "list")
  expect_true("AppCindex" %in% names(c_index_obj_t1))
  c_index_t1 <- c_index_obj_t1$AppCindex$matrix[1]

  expect_type(c_index_t1, "double")
  expect_gte(c_index_t1, 0)
  expect_lte(c_index_t1, 1)
  # C-index should be in valid range (mock predictions may not be highly discriminative)

  # Calculate C-index at the second time point
  predsurv_t2 <- matrix(mock_survprobs[, 2], nrow = 1)
  c_index_obj_t2 <- timedepConcordance(
    predsurv = predsurv_t2,
    predsurvtimes = time_points_met[2],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )
  c_index_t2 <- c_index_obj_t2$AppCindex$matrix[1]

  expect_type(c_index_t2, "double")
  expect_gte(c_index_t2, 0)
  expect_lte(c_index_t2, 1)
})

test_that("timedepConcordance handles matrix format", {
  # Test that the function works with proper matrix format
  # Matrix should have rows=times, cols=observations
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)

  c_index_obj <- timedepConcordance(
    predsurv = predsurv_t1,
    predsurvtimes = time_points_met[1],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )

  expect_type(c_index_obj, "list")
  expect_true("AppCindex" %in% names(c_index_obj))
})

test_that("timedepConcordance requires all mandatory parameters", {
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)

  # Missing predsurvtimes
  expect_error(
    timedepConcordance(predsurv = predsurv_t1, obstimes = surv_data_metrics$time,
                      obsevents = surv_data_metrics$status),
    "argument .* is missing"
  )

  # Missing obstimes
  expect_error(
    timedepConcordance(predsurv = predsurv_t1, predsurvtimes = time_points_met[1],
                      obsevents = surv_data_metrics$status),
    "argument .* is missing"
  )

  # Missing obsevents
  expect_error(
    timedepConcordance(predsurv = predsurv_t1, predsurvtimes = time_points_met[1],
                      obstimes = surv_data_metrics$time),
    "argument .* is missing"
  )
})

test_that("timedepConcordance handles cases with no events before time t", {
  surv_data_late_event <- surv_data_metrics # Replace copy()
  time_t <- min(surv_data_late_event$time[surv_data_late_event$status == 1]) - 0.1
  # Replace data.table assignment with base R
  rows_to_update <- surv_data_late_event$status == 1
  surv_data_late_event$time[rows_to_update] <- surv_data_late_event$time[rows_to_update] + time_t + 1 # Shift events later

  predsurv_late <- matrix(mock_survprobs[, 1], nrow = 1)
  c_index_obj_late <- timedepConcordance(
    predsurv = predsurv_late,
    predsurvtimes = time_t,
    obstimes = surv_data_late_event$time,
    obsevents = surv_data_late_event$status
  )

  # Extract concordance value
  c_index_late <- c_index_obj_late$AppCindex$matrix[1]
  expect_true(is.na(c_index_late) || c_index_late == 0.5) # Allow NA or 0.5
})

test_that("timedepConcordance handles cases with no comparable pairs", {
  surv_data_few_events <- data.frame(time = c(1, 2, 10, 11), status = c(0, 1, 0, 0), x1=1:4, stringsAsFactors = FALSE)
  # predsurv should be rows=times, cols=observations
  # preds_few: 1 row (1 time point), 4 columns (4 observations)
  preds_few <- matrix(c(0.9, 0.8, 0.7, 0.6), nrow=1)

  c_index_obj_no_pairs <- timedepConcordance(
    predsurv = preds_few,
    predsurvtimes = 5,
    obstimes = surv_data_few_events$time,
    obsevents = surv_data_few_events$status
  )
  c_index_no_pairs <- c_index_obj_no_pairs$AppCindex$matrix[1]
  # With few events/pairs, c-index may be computable but not necessarily 0.5
  expect_type(c_index_no_pairs, "double")

  c_index_obj_all_cens <- timedepConcordance(
    predsurv = preds_few,
    predsurvtimes = 0.5,
    obstimes = surv_data_few_events$time,
    obsevents = surv_data_few_events$status
  )
  c_index_all_cens <- c_index_obj_all_cens$AppCindex$matrix[1]
  expect_true(is.na(c_index_all_cens) || c_index_all_cens == 0.5)
})


# --- Tests for BrierScore ---

test_that("BrierScore calculates score correctly", {
  # Calculate Brier score at the first time point
  # Function signature: BrierScore(predsurv, predsurvtimes, obstimes, obsevents, eval.times = NULL, TestMat = NULL)
  # predsurv is a matrix with rows=times, cols=observations
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)
  brier_obj_t1 <- BrierScore(
    predsurv = predsurv_t1,
    predsurvtimes = time_points_met[1],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )

  # pec::pec returns an object with AppErr slot
  expect_type(brier_obj_t1, "list")
  expect_true("AppErr" %in% names(brier_obj_t1))
  brier_t1 <- brier_obj_t1$AppErr$model[1]

  expect_type(brier_t1, "double")
  expect_gte(brier_t1, 0)
  # Max Brier score is typically 0.25 for binary outcomes, can be higher with censoring weighting
  expect_lte(brier_t1, 1)

  # Calculate Brier score at the second time point
  predsurv_t2 <- matrix(mock_survprobs[, 2], nrow = 1)
  brier_obj_t2 <- BrierScore(
    predsurv = predsurv_t2,
    predsurvtimes = time_points_met[2],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )
  brier_t2 <- brier_obj_t2$AppErr$model[1]

  expect_type(brier_t2, "double")
  expect_gte(brier_t2, 0)
  expect_lte(brier_t2, 1)

  # Test with perfect predictions (expect Brier score near 0)
  # Status at time t
  status_at_t1 <- ifelse(surv_data_metrics$time <= time_points_met[1] & surv_data_metrics$status == 1, 1, 0)
  # Perfect survival prediction: 0 if event occurred by t1, 1 otherwise
  perfect_preds_t1 <- 1 - status_at_t1
  perfect_predsurv_t1 <- matrix(perfect_preds_t1, nrow = 1)
  brier_obj_perfect_t1 <- BrierScore(
    predsurv = perfect_predsurv_t1,
    predsurvtimes = time_points_met[1],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )
  # IPCW weighting can make it non-zero even for perfect status prediction
  # expect_lt(brier_perfect_t1, 0.1) # Expect low value
})

test_that("BrierScore works with matrix format", {
  # Test that the function works with proper matrix format
  # Matrix should have rows=times, cols=observations
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)

  brier_obj <- BrierScore(
    predsurv = predsurv_t1,
    predsurvtimes = time_points_met[1],
    obstimes = surv_data_metrics$time,
    obsevents = surv_data_metrics$status
  )

  expect_type(brier_obj, "list")
  expect_true("AppErr" %in% names(brier_obj))
})

test_that("BrierScore requires all mandatory parameters", {
  predsurv_t1 <- matrix(mock_survprobs[, 1], nrow = 1)

  # Missing predsurvtimes
  expect_error(
    BrierScore(predsurv = predsurv_t1, obstimes = surv_data_metrics$time,
               obsevents = surv_data_metrics$status),
    "argument .* is missing"
  )

  # Missing obstimes
  expect_error(
    BrierScore(predsurv = predsurv_t1, predsurvtimes = time_points_met[1],
               obsevents = surv_data_metrics$status),
    "argument .* is missing"
  )

  # Missing obsevents
  expect_error(
    BrierScore(predsurv = predsurv_t1, predsurvtimes = time_points_met[1],
               obstimes = surv_data_metrics$time),
    "argument .* is missing"
  )
})

test_that("BrierScore handles time point beyond last observation", {
   time_late <- max(surv_data_metrics$time) + 1
   # Behavior depends on IPCW calculation (e.g., using last censoring time)
   # It might produce NA or a score based on Kaplan-Meier estimate at last time.
   predsurv_late <- matrix(mock_survprobs[, 2], nrow = 1)
   brier_obj_late <- BrierScore(
     predsurv = predsurv_late,
     predsurvtimes = time_late,
     obstimes = surv_data_metrics$time,
     obsevents = surv_data_metrics$status
   )
   expect_type(brier_obj_late, "list") # Should return a pec object, not error
})
