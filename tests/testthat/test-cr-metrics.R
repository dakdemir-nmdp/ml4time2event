library(testthat)
library(here)
# library(data.table) # Removed
library(survival)
# library(pec) # Often used for C-index calculations, might be a dependency

# Assuming the function is available in the environment
source(here("R/cr_metrics.R"))  # Ensure timedepConcordanceCR is available when running tests standalone

context("Testing cr_metrics functions")

# --- Test Data Setup ---
# Reusing the setup from previous CR tests
set.seed(456) # Use different seed for potentially different outcomes
n_obs <- 60 # Slightly larger dataset for metrics
cr_data_metrics <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)
# Create Surv object
surv_obj <- Surv(cr_data_metrics$time, cr_data_metrics$status, type = "mstate")

# Create mock predictions (CIFs) for two causes at specific time points
# Predictions should ideally have some correlation with outcome for C-index to be meaningful
# Let's create predictions based loosely on x1
time_points_metrics <- quantile(cr_data_metrics$time[cr_data_metrics$status != 0], c(0.3, 0.6))

# Mock CIFs for cause 1 (higher for higher x1)
cif1_t1 <- pmax(0, pmin(1, 0.1 + 0.1 * cr_data_metrics$x1 + rnorm(n_obs, 0, 0.05)))
cif1_t2 <- pmax(0, pmin(1, 0.2 + 0.15 * cr_data_metrics$x1 + rnorm(n_obs, 0, 0.08)))
mock_cif1 <- rbind(cif1_t1, cif1_t2)
attr(mock_cif1, "Times") <- time_points_metrics

# Mock CIFs for cause 2 (higher for lower x1)
cif2_t1 <- pmax(0, pmin(1, 0.15 - 0.05 * cr_data_metrics$x1 + rnorm(n_obs, 0, 0.05)))
cif2_t2 <- pmax(0, pmin(1, 0.25 - 0.08 * cr_data_metrics$x1 + rnorm(n_obs, 0, 0.08)))
mock_cif2 <- rbind(cif2_t1, cif2_t2)
attr(mock_cif2, "Times") <- time_points_metrics

# Ensure CIF1 + CIF2 <= 1 (crude adjustment)
total_cif <- mock_cif1 + mock_cif2
scale_factor <- ifelse(total_cif > 1, 1 / (total_cif + 0.01), 1) # Add small epsilon
mock_cif1 <- mock_cif1 * scale_factor
mock_cif2 <- mock_cif2 * scale_factor

mock_predictions <- list(mock_cif1, mock_cif2)


# --- Tests for timedepConcordanceCR ---

test_that("timedepConcordanceCR calculates C-index for specified cause", {
  # skip_if_not_installed("pec") # If pec::cindex is used internally

  # Calculate C-index for cause 1 at the first time point
  # Note: timedepConcordanceCR expects predictions for *one* cause at a time
  # The input 'Predictions' should be the matrix for the cause of interest
  c_index_cause1_t1 <- timedepConcordanceCR(
    surv_obj,
    Predictions = mock_predictions[[1]][1, , drop = FALSE],
    time = time_points_metrics[1],
    cause = 1
  )

  expect_type(c_index_cause1_t1, "double")
  expect_length(c_index_cause1_t1, 1)
  expect_gte(c_index_cause1_t1, 0)
  expect_lte(c_index_cause1_t1, 1)
  # Note: C-index may be < 0.5 due to random data correlation

  # Calculate C-index for cause 2 at the second time point
  c_index_cause2_t2 <- timedepConcordanceCR(
    surv_obj,
    Predictions = mock_predictions[[2]][2, , drop = FALSE],
    time = time_points_metrics[2],
    cause = 2
  )

  expect_type(c_index_cause2_t2, "double")
  expect_length(c_index_cause2_t2, 1)
  expect_gte(c_index_cause2_t2, 0)
  expect_lte(c_index_cause2_t2, 1)
   # Note: C-index may be < 0.5 due to random data correlation
})

test_that("timedepConcordanceCR handles different input formats", {
  # Test with vector prediction (should work if only one time point implicitly)
   c_index_vec <- timedepConcordanceCR(
     surv_obj,
     Predictions = mock_predictions[[1]][1, ],
     time = time_points_metrics[1],
     cause = 1,
     pred_times = time_points_metrics[1]
   )
   expect_equal(
     c_index_vec,
     timedepConcordanceCR(
       surv_obj,
       Predictions = mock_predictions[[1]][1, , drop = FALSE],
       time = time_points_metrics[1],
       cause = 1,
       pred_times = time_points_metrics[1]
     )
   )

  # Test with data.frame prediction
   pred_df <- as.data.frame(matrix(mock_predictions[[1]][1, ], nrow = 1))
   colnames(pred_df) <- paste0("pred_", seq_len(ncol(pred_df)))
   c_index_df <- timedepConcordanceCR(
     surv_obj,
     Predictions = pred_df,
     time = time_points_metrics[1],
     cause = 1,
     pred_times = time_points_metrics[1]
   )
   expect_equal(
     c_index_df,
     timedepConcordanceCR(
       surv_obj,
       Predictions = mock_predictions[[1]][1, , drop = FALSE],
       time = time_points_metrics[1],
       cause = 1,
       pred_times = time_points_metrics[1]
     )
   )
})


test_that("timedepConcordanceCR requires Surv object, Predictions, time, and cause", {
  expect_error(timedepConcordanceCR(Predictions = mock_predictions[[1]][1, ], time = time_points_metrics[1], cause = 1), "Surv")
  expect_error(timedepConcordanceCR(surv_obj, time = time_points_metrics[1], cause = 1), "Predictions")
  expect_error(timedepConcordanceCR(surv_obj, Predictions = mock_predictions[[1]][1, ], cause = 1), "time")
  expect_error(timedepConcordanceCR(surv_obj, Predictions = mock_predictions[[1]][1, ], time = time_points_metrics[1]), "cause")
})

test_that("timedepConcordanceCR handles cases with no events of interest before time t", {
  # Create data where cause 1 only happens after time t
  cr_data_late_event <- cr_data_metrics # Use standard assignment for copy
  time_t <- min(cr_data_late_event$time[cr_data_late_event$status == 1]) - 0.1
  # Replace data.table assignment with base R subsetting and assignment
  rows_to_update <- cr_data_late_event$status == 1
  cr_data_late_event$time[rows_to_update] <- cr_data_late_event$time[rows_to_update] + time_t + 1 # Shift cause 1 events later
  surv_obj_late <- Surv(cr_data_late_event$time, cr_data_late_event$status, type = "mstate")

  # Predictions don't matter much here, expect NA or 0.5? Depends on implementation.
  # pec::cindex returns NA in this case.
  c_index_late <- timedepConcordanceCR(surv_obj_late, Predictions = mock_predictions[[1]][1, ],
                                       time = time_t, cause = 1, pred_times = time_points_metrics[1])
  expect_true(is.na(c_index_late) || c_index_late == 0.5) # Allow NA or 0.5
})

test_that("timedepConcordanceCR handles cases with no comparable pairs", {
  # Scenario 1: Only one subject with event - no comparable pairs
  cr_data_single <- data.frame(time = c(5), status = c(1), x1 = 1, stringsAsFactors = FALSE)
  surv_obj_single <- Surv(cr_data_single$time, cr_data_single$status)
  preds_single <- matrix(c(0.5), nrow=1)
  
  c_index_single <- timedepConcordanceCR(surv_obj_single, Predictions = preds_single, time = 10, cause = 1, pred_times = 5)
  expect_true(is.na(c_index_single)) # Should return NA when no comparable pairs
  
  # Scenario 2: Event happens after all other observations - no comparable pairs
  cr_data_late_event <- data.frame(time = c(1, 2, 3, 10), status = c(0, 0, 0, 1), x1=1:4, stringsAsFactors = FALSE)
  surv_obj_late <- Surv(cr_data_late_event$time, cr_data_late_event$status)
  preds_late <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=1)
  
  c_index_late <- timedepConcordanceCR(surv_obj_late, Predictions = preds_late, time = 15, cause = 1, pred_times = 5)
  expect_true(is.na(c_index_late)) # Should return NA when no comparable pairs

  # All censored before time t - no events to evaluate
  cr_data_early_cens <- data.frame(time = c(1, 2, 10, 11), status = c(0, 1, 0, 0), x1=1:4, stringsAsFactors = FALSE)
  surv_obj_early <- Surv(cr_data_early_cens$time, cr_data_early_cens$status)
  preds_early <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=1)
  
  c_index_all_cens <- timedepConcordanceCR(surv_obj_early, Predictions = preds_early, time = 0.5, cause = 1, pred_times = 0.3)
  expect_true(is.na(c_index_all_cens)) # Should return NA when no events before time t
})

test_that("timedepConcordanceCR handles cases with zero concordance", {
  # Case where event has lower prediction than all comparable subjects (all discordant)
  cr_data_discordant <- data.frame(time = c(1, 2, 10, 11), status = c(0, 1, 0, 0), x1=1:4, stringsAsFactors = FALSE)
  surv_obj_discordant <- Surv(cr_data_discordant$time, cr_data_discordant$status)
  preds_discordant <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=1)  # Event at t=2 has pred=0.2, others at t=10,11 have pred=0.3,0.4
  
  c_index_zero <- timedepConcordanceCR(surv_obj_discordant, Predictions = preds_discordant, time = 5, cause = 1, pred_times = 2)
  expect_equal(c_index_zero, 0) # Should return 0 when all pairs are discordant
  expect_false(is.na(c_index_zero)) # Should not be NA since there are comparable pairs
})
