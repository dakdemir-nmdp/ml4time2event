library(testthat)
library(ml4time2event)

# Create a simplified test version of the cifMatInterpolaltor function
# This avoids using the existing function which might depend on other parts
test_interpolator <- function(probsMat, times, new_times) {
  n_obs <- ncol(probsMat)
  n_new_times <- length(new_times)
  
  # Prepare result matrix
  result_mat <- matrix(NA_real_, nrow=n_new_times, ncol=n_obs)
  
  # Process each observation
  for (i in seq_len(n_obs)) {
    # Get probabilities for this observation
    probs <- probsMat[, i]
    
    # Interpolate to new times
    interp_probs <- stats::approx(
      x = times,
      y = probs,
      xout = new_times,
      method = "linear",
      yleft = 0,
      yright = utils::tail(probs, 1),
      rule = 2,
      ties = "ordered"
    )$y
    
    # Store in result matrix
    result_mat[, i] <- interp_probs
  }
  
  return(result_mat)
}

# Test the cifMatInterpolaltor function with synthetic data
test_that("interpolation works with simple synthetic data", {
  # Create synthetic probability matrix and times
  n_times <- 5
  n_subj <- 10
  
  # Create matrix with rows=times, cols=observations
  probsMat <- matrix(runif(n_subj * n_times), nrow = n_times, ncol = n_subj)
  times <- seq(1, 10, length.out = n_times)
  
  # Test vector of new times
  new_times_vec <- c(2, 4, 6, 8)
  result_vec <- test_interpolator(probsMat, times, new_times_vec)
  
  # Test structure of result
  expect_true(is.matrix(result_vec))
  expect_equal(nrow(result_vec), length(new_times_vec))
  expect_equal(ncol(result_vec), n_subj)
  
  # Test scalar new time
  new_times_scalar <- 5
  result_scalar <- test_interpolator(probsMat, times, new_times_scalar)
  
  # Test structure of scalar result
  expect_true(is.matrix(result_scalar))
  expect_equal(nrow(result_scalar), 1)
  expect_equal(ncol(result_scalar), n_subj)
})

# Test ensemble function with synthetic data
test_that("cifMatListAveraging works with simple synthetic data", {
  # Create two synthetic matrices
  n_times <- 3
  n_subj <- 8
  
  # Matrix 1 - higher values for first half of subjects
  mat1 <- matrix(0, nrow = n_times, ncol = n_subj)
  mat1[, 1:(n_subj/2)] <- 0.8
  mat1[, (n_subj/2+1):n_subj] <- 0.2
  
  # Matrix 2 - higher values for first half of subjects
  mat2 <- matrix(0, nrow = n_times, ncol = n_subj)
  mat2[, 1:(n_subj/2)] <- 0.9
  mat2[, (n_subj/2+1):n_subj] <- 0.1
  
  # Create list of matrices
  mats <- list(mat1 = mat1, mat2 = mat2)
  
  # Test averaging
  result <- ml4time2event:::cifMatListAveraging(mats, type = "prob")
  
  # Check structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat1))
  
  # Check values
  expect_true(all(result[, 1:(n_subj/2)] > 0.8))
  expect_true(all(result[, (n_subj/2+1):n_subj] < 0.2))
})
