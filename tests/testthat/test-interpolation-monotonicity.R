test_that("survprobMatInterpolator preserves non-monotonicity without forcing it", {
  # Create a non-monotonic survival curve to test the removal of the `cummin` fix
  times <- c(0, 10, 20, 30)
  # Non-monotonic probabilities for the first observation (0.8 -> 0.85)
  probs_mat <- matrix(c(1, 0.8, 0.85, 0.7, 1, 0.9, 0.8, 0.7), nrow = 4, ncol = 2)
  new_times <- c(5, 15, 25)

  # Expected output from pure linear interpolation of the non-monotonic curve
  # For the first observation, the value at time 15 should be 0.825, which is > 0.8
  expected_interp <- matrix(c(0.9, 0.825, 0.775, 0.95, 0.85, 0.75), nrow = 3, ncol = 2)

  # The `cummin` fix would change the first column to c(0.9, 0.8, 0.775)
  # We want to ensure the raw, non-monotonic interpolation is returned.
  
  # Get result from function
  result_interp <- survprobMatInterpolator(probs_mat, times, new_times)

  # This test will fail if cummin() is still active, as it would alter the first column's result
  expect_equal(result_interp, expected_interp)
})

