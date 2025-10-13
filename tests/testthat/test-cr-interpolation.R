test_that("cifMatInterpolaltor preserves non-monotonicity and propagates NAs", {
  # Define times and a matrix with non-monotonic and NA values
  times <- c(10, 20, 30)
  newtimes <- c(15, 25)
  
  # Observation 1 is non-monotonic (0.15 -> 0.12)
  # Observation 2 contains an NA
  probs_mat <- matrix(c(
    0.1, 0.15, 0.12,
    0.2, NA,   0.3
  ), nrow = 2, byrow = TRUE)
  
  # Expected result from pure linear interpolation:
  # Obs 1, time 15: 0.125
  # Obs 1, time 25: 0.135 (non-monotonic result)
  # Obs 2 should be all NA because of the NA in the input
  
  # The function `apply` will transpose the matrix, so we expect a 2x2 matrix
  # with columns for observations and rows for newtimes.
  expected_output <- matrix(c(
    0.125, 0.135, # obs 1
    NA,    NA     # obs 2
  ), nrow = 2, ncol = 2, byrow = FALSE)
  
  # Get the actual output from the function
  actual_output <- cifMatInterpolaltor(probs_mat, times, newtimes)
  
  # The test will fail if the function forces monotonicity or does not propagate NAs.
  expect_equal(actual_output, expected_output)
})
