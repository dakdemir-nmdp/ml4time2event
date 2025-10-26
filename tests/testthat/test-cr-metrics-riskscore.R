# Source the concordance function robustly
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("R", "cr_metrics.R"))

# Test for concordance function with risk scores

test_that("timedepConcordanceCR handles risk scores correctly (perfect order)", {
  library(survival)
  n <- 10
  times <- seq(1, n) # perfectly ordered event times
  status <- rep(1, n) # all events
  SurvObj <- Surv(times, status)
  # Risk scores: higher = earlier event
  risk_scores <- rev(seq(1, n)) # highest risk for earliest event
  Predictions <- matrix(risk_scores, nrow = 1)
  pred_time <- n
  cidx <- timedepConcordanceCR(SurvObj, Predictions, pred_time, cause = 1, pred_times = pred_time)
  expect_true(cidx > 0.99)

  # Reverse risk scores: lowest risk for earliest event
  Predictions_rev <- matrix(seq(1, n), nrow = 1)
  cidx_rev <- timedepConcordanceCR(SurvObj, Predictions_rev, pred_time, cause = 1, pred_times = pred_time)
  expect_true(cidx_rev < 0.01)
})
