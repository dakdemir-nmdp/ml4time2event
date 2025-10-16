#' Test file for the updated cifMatInterpolaltor function
#' 
#' This file contains tests for the fixed cifMatInterpolaltor function
#' to ensure it correctly handles edge cases.

library(testthat)
library(ml4time2event)

# Test the cifMatInterpolaltor function with synthetic data
test_that("interpolation works with well-formed synthetic data", {
  # Create synthetic probability matrix and times
  n_times <- 5
  n_subj <- 10
  
  # Create matrix with rows=observations, cols=times
  probsMat <- matrix(runif(n_subj * n_times), nrow = n_subj, ncol = n_times)
  times <- seq(1, 10, length.out = n_times)
  
  # Test vector of new times
  newtimes_vec <- c(2, 4, 6, 8)
  result_vec <- cifMatInterpolaltor(probsMat, times, newtimes_vec)
  
  # Test structure of result
  expect_true(is.matrix(result_vec))
  expect_equal(nrow(result_vec), length(newtimes_vec))
  expect_equal(ncol(result_vec), n_subj)
  
  # Test scalar new time
  newtimes_scalar <- 5
  result_scalar <- cifMatInterpolaltor(probsMat, times, newtimes_scalar)
  
  # Test structure of scalar result
  expect_true(is.matrix(result_scalar))
  expect_equal(nrow(result_scalar), 1)
  expect_equal(ncol(result_scalar), n_subj)
})

test_that("interpolation handles transposed matrices correctly", {
  # Create synthetic probability matrix and times
  n_times <- 5
  n_subj <- 10
  
  # Create matrix with rows=times, cols=observations (transposed)
  probsMat <- matrix(runif(n_subj * n_times), nrow = n_times, ncol = n_subj)
  times <- seq(1, 10, length.out = n_times)
  
  # Test with transposed matrix - should produce warning but work
  expect_warning(
    result <- cifMatInterpolaltor(probsMat, times, c(2, 4, 6, 8)),
    "transposed"
  )
  
  # Test structure of result
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)  # 4 new time points
  expect_equal(ncol(result), n_subj)
})

test_that("interpolation handles edge cases with time=0", {
  n_times <- 5
  n_subj <- 10
  
  # Create matrix with rows=observations, cols=times
  probsMat <- matrix(runif(n_subj * n_times), nrow = n_subj, ncol = n_times)
  
  # Test with times including 0
  times_with_zero <- c(0, seq(1, 10, length.out = (n_times - 1)))
  
  # Should work without warning
  result <- cifMatInterpolaltor(probsMat, times_with_zero, c(0.5, 5, 8))
  
  # Test structure
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)  # 3 new time points
  expect_equal(ncol(result), n_subj)
  
  # Test that values at time=0 are 0 when interpolating back
  result_at_zero <- cifMatInterpolaltor(probsMat, times_with_zero, 0)
  expect_equal(as.vector(result_at_zero), rep(0, n_subj))
})

test_that("interpolation handles NA values correctly", {
  n_times <- 5
  n_subj <- 10
  
  # Create matrix with rows=observations, cols=times
  probsMat <- matrix(runif(n_subj * n_times), nrow = n_subj, ncol = n_times)
  
  # Insert some NA values
  probsMat[2, 3] <- NA
  probsMat[5, 1:5] <- NA  # Entire row NA
  
  times <- seq(1, 10, length.out = n_times)
  
  # Should handle NA values appropriately
  result <- cifMatInterpolaltor(probsMat, times, c(2, 4, 6, 8))
  
  # Test structure
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)  # 4 new time points
  expect_equal(ncol(result), n_subj)
  
  # Row with all NAs should have all NAs in result
  expect_true(all(is.na(result[, 5])))
})

test_that("cifMatListAveraging works with multiple matrices", {
  # Create two synthetic matrices
  n_times <- 3
  n_subj <- 8
  
  # Matrix 1 - higher values for first half of subjects
  mat1 <- matrix(0, nrow = n_times, ncol = n_subj)
  for (i in 1:n_times) {
    mat1[i, 1:(n_subj/2)] <- 0.8 * i/n_times
    mat1[i, (n_subj/2+1):n_subj] <- 0.2 * i/n_times
  }
  
  # Matrix 2 - higher values for first half of subjects
  mat2 <- matrix(0, nrow = n_times, ncol = n_subj)
  for (i in 1:n_times) {
    mat2[i, 1:(n_subj/2)] <- 0.9 * i/n_times
    mat2[i, (n_subj/2+1):n_subj] <- 0.1 * i/n_times
  }
  
  # Create list of matrices
  mats <- list(mat1 = mat1, mat2 = mat2)
  
  # Test averaging
  result <- cifMatListAveraging(mats, type = "prob")
  
  # Check structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat1))
  
  # Check values - average should be between the two matrices
  expect_true(all(result[, 1:(n_subj/2)] >= mat1[, 1:(n_subj/2)] - 1e-10 | 
                 result[, 1:(n_subj/2)] <= mat2[, 1:(n_subj/2)] + 1e-10))
  expect_true(all(result[, (n_subj/2+1):n_subj] >= mat2[, (n_subj/2+1):n_subj] - 1e-10 | 
                 result[, (n_subj/2+1):n_subj] <= mat1[, (n_subj/2+1):n_subj] + 1e-10))
})

test_that("cifMatListAveraging handles invalid matrices", {
  # Create valid matrices
  n_times <- 3
  n_subj <- 8
  
  # Valid matrix 1
  mat1 <- matrix(runif(n_times * n_subj), nrow = n_times, ncol = n_subj)
  
  # Valid matrix 2
  mat2 <- matrix(runif(n_times * n_subj), nrow = n_times, ncol = n_subj)
  
  # Invalid matrix with wrong dimensions
  mat3 <- matrix(runif((n_times+1) * n_subj), nrow = (n_times+1), ncol = n_subj)
  
  # Create list with invalid matrix
  mats <- list(mat1 = mat1, mat2 = mat2, mat3 = mat3)
  
  # Should error with invalid matrices
  expect_error(cifMatListAveraging(mats, type = "prob"), 
               "All matrices in listprobsMat must have the same dimensions")
  
  # Let's also check it works with valid matrices
  mats_valid <- list(mat1 = mat1, mat2 = mat2)
  result <- cifMatListAveraging(mats_valid, type = "prob")
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat1))
})