library(testthat)
library(here)
# library(data.table) # Removed
library(stats) # For approxfun

# Assuming the functions are available in the environment
source(here("R/utils/surv_interpolation.R"))

context("Testing surv_interpolation functions")

# --- Test Data Setup ---
# For survivalProbsInterpolator
times_single <- c(1, 3, 5, 8)
probs_single <- c(0.9, 0.7, 0.6, 0.5)
new_times_single <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

# For survprobMatInterpolator & survprobMatListAveraging
times_mat <- c(2, 5, 10)
probs_mat1 <- matrix(c(
  0.9, 0.7, 0.5, # Subject 1
  0.8, 0.6, 0.3  # Subject 2
), nrow = 2, byrow = TRUE)
probs_mat2 <- matrix(c(
  0.85, 0.75, 0.55, # Subject 1
  0.82, 0.55, 0.35  # Subject 2
), nrow = 2, byrow = TRUE)

new_times_mat <- c(0, 1, 2, 3, 4, 5, 8, 10, 12)


# --- Tests for survivalProbsInterpolator ---

test_that("survivalProbsInterpolator interpolates correctly (constant, right-continuous)", {
  interpolated_probs <- survivalProbsInterpolator(new_times_single, probs_single, times_single)
  # Expected: yleft=1, step down at times_single, yright=min(probs)=0.5
  expected_probs <- c(1, 0.9, 0.9, 0.7, 0.7, 0.6, 0.6, 0.6, 0.5, 0.5)
  expect_equal(interpolated_probs, expected_probs)
})

test_that("survivalProbsInterpolator handles yleft and yright", {
  # Test time before first time point
  expect_equal(survivalProbsInterpolator(0.5, probs_single, times_single), 1)
  # Test time after last time point
  expect_equal(survivalProbsInterpolator(10, probs_single, times_single), min(probs_single))
})

test_that("survivalProbsInterpolator handles unsorted input", {
  times_unsorted <- c(5, 1, 8, 3)
  probs_unsorted <- c(0.6, 0.9, 0.5, 0.7) # Corresponds to 1=0.9, 3=0.7, 5=0.6, 8=0.5
  interpolated_probs <- survivalProbsInterpolator(new_times_single, probs_unsorted, times_unsorted)
  expected_probs <- c(1, 0.9, 0.9, 0.7, 0.7, 0.6, 0.6, 0.6, 0.5, 0.5)
  expect_equal(interpolated_probs, expected_probs)
})

test_that("survivalProbsInterpolator handles single time/prob input", {
   interpolated <- survivalProbsInterpolator(c(0, 5, 10), 0.8, 5)
   expect_equal(interpolated, c(1, 0.8, 0.8))
})

test_that("survivalProbsInterpolator handles NA in probs (uses min of non-NA for yright)", {
   probs_na <- c(0.9, NA, 0.6, 0.5)
   interpolated <- survivalProbsInterpolator(c(9, 10), probs_na, times_single)
   expect_equal(interpolated, c(0.5, 0.5)) # yright should be min(0.9, 0.6, 0.5) = 0.5
})


# --- Tests for survprobMatInterpolator ---

test_that("survprobMatInterpolator interpolates matrix correctly", {
  interpolated_mat <- survprobMatInterpolator(probs_mat1, times_mat, new_times_mat)

  # Check dimensions: rows = newtimes, cols = observations
  expect_equal(nrow(interpolated_mat), length(new_times_mat))
  expect_equal(ncol(interpolated_mat), nrow(probs_mat1))

  # Check interpolation for subject 1 (probs: 0.9@t=2, 0.7@t=5, 0.5@t=10)
  # Expected at new_times_mat = c(0, 1, 2, 3, 4, 5, 8, 10, 12)
  # Should be: c(1, 1, 0.9, 0.9, 0.9, 0.7, 0.7, 0.5, 0.5)
  expect_equal(interpolated_mat[, 1], c(1, 1, 0.9, 0.9, 0.9, 0.7, 0.7, 0.5, 0.5))

  # Check interpolation for subject 2 (probs: 0.8@t=2, 0.6@t=5, 0.3@t=10)
  # Expected at new_times_mat = c(0, 1, 2, 3, 4, 5, 8, 10, 12)
  # Should be: c(1, 1, 0.8, 0.8, 0.8, 0.6, 0.6, 0.3, 0.3)
  expect_equal(interpolated_mat[, 2], c(1, 1, 0.8, 0.8, 0.8, 0.6, 0.6, 0.3, 0.3))
})

test_that("survprobMatInterpolator handles time 0 correctly", {
  # Case 1: time 0 not in input times (tested above)
  interpolated_mat <- survprobMatInterpolator(probs_mat1, times_mat, new_times_mat)
  expect_equal(interpolated_mat[1, ], c(1, 1)) # First row (time 0) should be 1

  # Case 2: time 0 is in input times
  times_with_zero <- c(0, times_mat)
  probs_with_zero <- cbind(rep(1, nrow(probs_mat1)), probs_mat1)
  interpolated_mat_zero <- survprobMatInterpolator(probs_with_zero, times_with_zero, new_times_mat)
  expect_equal(interpolated_mat_zero[, 1], c(1, 1, 0.9, 0.9, 0.9, 0.7, 0.7, 0.5, 0.5))
  expect_equal(interpolated_mat_zero[, 2], c(1, 1, 0.8, 0.8, 0.8, 0.6, 0.6, 0.3, 0.3))
})

test_that("survprobMatInterpolator enforces monotonicity", {
  # Create a matrix where interpolation might initially increase
  probs_nonmono <- matrix(c(0.9, 0.7, 0.8), nrow = 1) # Increases from 0.7 to 0.8
  times_nonmono <- c(2, 5, 10)
  new_times_nonmono <- c(1, 3, 6, 11)
  # Initial interpolation: 1, 0.9, 0.7, 0.8
  # Monotonicity correction should make it: 1, 0.9, 0.7, 0.7
  interpolated_mat <- survprobMatInterpolator(probs_nonmono, times_nonmono, new_times_nonmono)
  expect_equal(as.vector(interpolated_mat), c(1, 0.9, 0.7, 0.7))
})

test_that("survprobMatInterpolator handles single new time", {
   interpolated_mat <- survprobMatInterpolator(probs_mat1, times_mat, newtimes = 4)
   expect_true(is.matrix(interpolated_mat))
   expect_equal(nrow(interpolated_mat), 1)
   expect_equal(ncol(interpolated_mat), nrow(probs_mat1))
   # Expected at t=4: Subj1=0.9, Subj2=0.8
   expect_equal(as.vector(interpolated_mat), c(0.9, 0.8))
})


# --- Tests for survprobMatListAveraging ---

test_that("survprobMatListAveraging averages correctly on cumulative hazard scale", {
  list_mats <- list(probs_mat1, probs_mat2)
  # Transpose because survprobMatInterpolator output has times as rows
  list_mats_transposed <- list(t(survprobMatInterpolator(probs_mat1, times_mat, new_times_mat)),
                               t(survprobMatInterpolator(probs_mat2, times_mat, new_times_mat)))

  averaged_mat <- survprobMatListAveraging(list_mats_transposed)

  # Check dimensions (should match input matrix dims: newtimes x observations)
  expect_equal(dim(averaged_mat), dim(list_mats_transposed[[1]]))

  # Check a specific value (e.g., Subject 1 at newtime = 8)
  # Mat1: S(8) = 0.7
  # Mat2: S(8) = 0.55 (interpolated between 0.75@t=5 and 0.55@t=10 -> step down at 10 -> 0.75)
  # Let's re-run interpolation for mat2, subj 1: times=c(0,2,5,10), probs=c(1,0.85,0.75,0.55) -> S(8)=0.75
  # H1 = -log(0.7) ~= 0.3567
  # H2 = -log(0.75) ~= 0.2877
  # Mean H = (0.3567 + 0.2877) / 2 ~= 0.3222
  # Averaged S = exp(-0.3222) ~= 0.7246
  subj1_idx <- 1
  time8_idx <- which(new_times_mat == 8)
  expect_equal(averaged_mat[time8_idx, subj1_idx], exp(-( -log(0.7) + (-log(0.75)) ) / 2), tolerance = 1e-6)

  # Check another value (Subject 2 at newtime = 12)
  # Mat1: S(12) = 0.3
  # Mat2: S(12) = 0.35
  # H1 = -log(0.3) ~= 1.2040
  # H2 = -log(0.35) ~= 1.0498
  # Mean H = (1.2040 + 1.0498) / 2 ~= 1.1269
  # Averaged S = exp(-1.1269) ~= 0.3240
  subj2_idx <- 2
  time12_idx <- which(new_times_mat == 12)
   expect_equal(averaged_mat[time12_idx, subj2_idx], exp(-( -log(0.3) + (-log(0.35)) ) / 2), tolerance = 1e-6)
})

test_that("survprobMatListAveraging handles list with one matrix", {
  list_one <- list(t(survprobMatInterpolator(probs_mat1, times_mat, new_times_mat)))
  averaged_mat <- survprobMatListAveraging(list_one)
  expect_equal(averaged_mat, list_one[[1]])
})

test_that("survprobMatListAveraging handles empty list", {
  expect_null(survprobMatListAveraging(list()))
})

test_that("survprobMatListAveraging handles inconsistent dimensions", {
  mat_wrong_dim <- matrix(1:6, nrow=3) # Different rows
  list_inconsistent <- list(t(survprobMatInterpolator(probs_mat1, times_mat, new_times_mat)),
                            mat_wrong_dim)
  expect_error(survprobMatListAveraging(list_inconsistent), "must have the same dimensions")
})

test_that("survprobMatListAveraging handles NAs (na.omit behavior)", {
   mat1 <- matrix(c(0.8, 0.6, 0.7, 0.5), nrow=2) # times x obs
   mat2 <- matrix(c(0.7, NA, 0.6, 0.4), nrow=2)
   mat3 <- matrix(c(NA, NA, 0.5, 0.3), nrow=2)
   list_na <- list(mat1, mat2, mat3)

   averaged <- survprobMatListAveraging(list_na)

   # Check element [1, 1]: mean(-log(0.8), -log(0.7)) = exp(-(-0.2231 - 0.3567)/2) = exp(0.2899) = 0.748
   expect_equal(averaged[1,1], exp(-mean(c(-log(0.8), -log(0.7)))))
   # Check element [1, 2]: only mat1 has value -> should be mat1 value
   expect_equal(averaged[1,2], mat1[1,2])
   # Check element [2, 1]: mean(-log(0.7), -log(0.6), -log(0.5))
   expect_equal(averaged[2,1], exp(-mean(c(-log(0.7), -log(0.6), -log(0.5)))))
   # Check element [2, 2]: mean(-log(0.5), -log(0.4), -log(0.3))
   expect_equal(averaged[2,2], exp(-mean(c(-log(0.5), -log(0.4), -log(0.3)))))
})
