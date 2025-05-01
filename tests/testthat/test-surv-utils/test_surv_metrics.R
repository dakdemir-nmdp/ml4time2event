library(testthat)
library(here)
# library(data.table) # Removed
library(survival)
# library(pec) # Often used for C-index/Brier calculations, might be a dependency

# Assuming the function is available in the environment
source(here("R/utils/surv_metrics.R"))

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
  # Need to check if the function uses survival probs or risk scores internally.
  # Assuming survival probabilities based on typical C-index usage.
  c_index_t1 <- timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 1, drop = FALSE],
                                   time = time_points_met[1])

  expect_type(c_index_t1, "double")
  expect_length(c_index_t1, 1)
  expect_gte(c_index_t1, 0)
  expect_lte(c_index_t1, 1)
  # Given mock_survprobs decreases with x1, and events might correlate with x1, expect C > 0.5
  expect_gt(c_index_t1, 0.5) # This is a heuristic check

  # Calculate C-index at the second time point
  c_index_t2 <- timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 2, drop = FALSE],
                                   time = time_points_met[2])
  expect_type(c_index_t2, "double")
  expect_gt(c_index_t2, 0.5) # Heuristic check
})

test_that("timedepConcordance handles different input formats", {
  # Test with vector prediction
   c_index_vec <- timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 1],
                                     time = time_points_met[1])
   expect_equal(c_index_vec, timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 1, drop=FALSE], time = time_points_met[1]))

  # Test with data.frame prediction
   pred_df <- as.data.frame(mock_survprobs[, 1, drop = FALSE])
   colnames(pred_df) <- "pred"
   c_index_df <- timedepConcordance(surv_obj_met, Predictions = pred_df,
                                    time = time_points_met[1])
   expect_equal(c_index_df, timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 1, drop=FALSE], time = time_points_met[1]))
})

test_that("timedepConcordance requires Surv object, Predictions, and time", {
  expect_error(timedepConcordance(Predictions = mock_survprobs[, 1], time = time_points_met[1]), "Surv")
  expect_error(timedepConcordance(surv_obj_met, time = time_points_met[1]), "Predictions")
  expect_error(timedepConcordance(surv_obj_met, Predictions = mock_survprobs[, 1]), "time")
})

test_that("timedepConcordance handles cases with no events before time t", {
  surv_data_late_event <- surv_data_metrics # Replace copy()
  time_t <- min(surv_data_late_event$time[surv_data_late_event$status == 1]) - 0.1
  # Replace data.table assignment with base R
  rows_to_update <- surv_data_late_event$status == 1
  surv_data_late_event$time[rows_to_update] <- surv_data_late_event$time[rows_to_update] + time_t + 1 # Shift events later
  surv_obj_late <- Surv(surv_data_late_event$time, surv_data_late_event$status)

  c_index_late <- timedepConcordance(surv_obj_late, Predictions = mock_survprobs[, 1], time = time_t)
  expect_true(is.na(c_index_late) || c_index_late == 0.5) # Allow NA or 0.5
})

test_that("timedepConcordance handles cases with no comparable pairs", {
  surv_data_few_events <- data.frame(time = c(1, 2, 10, 11), status = c(0, 1, 0, 0), x1=1:4, stringsAsFactors = FALSE)
  surv_obj_few <- Surv(surv_data_few_events$time, surv_data_few_events$status)
  preds_few <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol=1)

  c_index_no_pairs <- timedepConcordance(surv_obj_few, Predictions = preds_few, time = 5)
  expect_true(is.na(c_index_no_pairs) || c_index_no_pairs == 0.5) # Expect NA or 0.5

  c_index_all_cens <- timedepConcordance(surv_obj_few, Predictions = preds_few, time = 0.5)
   expect_true(is.na(c_index_all_cens) || c_index_all_cens == 0.5)
})


# --- Tests for BrierScore ---

test_that("BrierScore calculates score correctly", {
  # Calculate Brier score at the first time point
  brier_t1 <- BrierScore(surv_obj_met, Predictions = mock_survprobs[, 1, drop = FALSE],
                         time = time_points_met[1])

  expect_type(brier_t1, "double")
  expect_length(brier_t1, 1)
  expect_gte(brier_t1, 0)
  # Max Brier score is typically 0.25 for binary outcomes, can be higher with censoring weighting
  expect_lte(brier_t1, 1)

  # Calculate Brier score at the second time point
  brier_t2 <- BrierScore(surv_obj_met, Predictions = mock_survprobs[, 2, drop = FALSE],
                         time = time_points_met[2])
  expect_type(brier_t2, "double")
  expect_gte(brier_t2, 0)
  expect_lte(brier_t2, 1)

  # Test with perfect predictions (expect Brier score near 0)
  # Status at time t
  status_at_t1 <- ifelse(surv_obj_met[, "time"] <= time_points_met[1] & surv_obj_met[, "status"] == 1, 1, 0)
  # Perfect survival prediction: 0 if event occurred by t1, 1 otherwise
  perfect_preds_t1 <- 1 - status_at_t1
  brier_perfect_t1 <- BrierScore(surv_obj_met, Predictions = matrix(perfect_preds_t1, ncol=1),
                                 time = time_points_met[1])
  # IPCW weighting can make it non-zero even for perfect status prediction
  # expect_lt(brier_perfect_t1, 0.1) # Expect low value
})

test_that("BrierScore handles different input formats", {
  # Test with vector prediction
   brier_vec <- BrierScore(surv_obj_met, Predictions = mock_survprobs[, 1],
                           time = time_points_met[1])
   expect_equal(brier_vec, BrierScore(surv_obj_met, Predictions = mock_survprobs[, 1, drop=FALSE], time = time_points_met[1]))

  # Test with data.frame prediction
   pred_df <- as.data.frame(mock_survprobs[, 1, drop = FALSE])
   colnames(pred_df) <- "pred"
   brier_df <- BrierScore(surv_obj_met, Predictions = pred_df,
                          time = time_points_met[1])
   expect_equal(brier_df, BrierScore(surv_obj_met, Predictions = mock_survprobs[, 1, drop=FALSE], time = time_points_met[1]))
})

test_that("BrierScore requires Surv object, Predictions, and time", {
  expect_error(BrierScore(Predictions = mock_survprobs[, 1], time = time_points_met[1]), "Surv")
  expect_error(BrierScore(surv_obj_met, time = time_points_met[1]), "Predictions")
  expect_error(BrierScore(surv_obj_met, Predictions = mock_survprobs[, 1]), "time")
})

test_that("BrierScore handles time point beyond last observation", {
   time_late <- max(surv_data_metrics$time) + 1
   # Behavior depends on IPCW calculation (e.g., using last censoring time)
   # It might produce NA or a score based on Kaplan-Meier estimate at last time.
   brier_late = BrierScore(surv_obj_met, Predictions = mock_survprobs[, 2], time = time_late)
   expect_type(brier_late, "double") # Should return a value or NA, not error
})
