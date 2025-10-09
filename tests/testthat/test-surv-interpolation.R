test_that("survprobMatListAveraging handles NA values correctly", {
  m1 <- matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2, ncol = 2)
  m2 <- matrix(c(0.8, 0.7, NA, 0.5), nrow = 2, ncol = 2)
  m3 <- matrix(c(0.85, 0.75, 0.65, 0.55), nrow = 2, ncol = 2)
  
  list_probs_mat <- list(m1, m2, m3)
  
  # The current implementation with na.omit will produce a result, hiding the NA.
  # The corrected implementation should propagate the NA.
  
  # Calculate expected result without na.omit
  h1 <- -log(m1 + 1e-10)
  h2 <- -log(m2 + 1e-10) # This will have an NA
  h3 <- -log(m3 + 1e-10)
  
  h_avg <- array(c(h1, h2, h3), dim = c(2, 2, 3))
  mean_h <- apply(h_avg, c(1, 2), mean, na.rm = FALSE)
  expected_probs <- exp(-mean_h)
  
  # Get actual result from the function
  result_probs <- survprobMatListAveraging(list_probs_mat)
  
  # The test will fail if the function does not propagate NAs.
  # We expect an NA where m2 has an NA.
  expect_true(is.na(result_probs[1, 2]))
  
  # For the non-NA value, it should be close to the expected value
  expect_equal(result_probs[1, 1], expected_probs[1, 1])
})
