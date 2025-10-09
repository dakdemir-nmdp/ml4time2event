library(testthat)
library(here)
# library(data.table) # Removed
library(stats) # For approxfun

# Assuming the functions are available in the environment
source(here("R/cr_interpolation.R"))

context("Testing cr_interpolation functions")

# --- Test Data Setup ---
# For cifInterpolator
times_single_cif <- c(1, 3, 5, 8)
probs_single_cif <- c(0.1, 0.3, 0.4, 0.5) # Non-decreasing CIF
new_times_single_cif <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

# For cifMatInterpolaltor & cifMatListAveraging
times_mat_cif <- c(2, 5, 10)
# Rows = observations, Cols = times
cif_mat1 <- matrix(c(
  0.1, 0.3, 0.5, # Subject 1
  0.2, 0.4, 0.7  # Subject 2
), nrow = 2, byrow = TRUE)
cif_mat2 <- matrix(c(
  0.15, 0.25, 0.45, # Subject 1
  0.18, 0.45, 0.65  # Subject 2
), nrow = 2, byrow = TRUE)

new_times_mat_cif <- c(0, 1, 2, 3, 4, 5, 8, 10, 12)


# --- Tests for cifInterpolator ---

test_that("cifInterpolator interpolates correctly (linear)", {
  interpolated_probs <- cifInterpolator(new_times_single_cif, probs_single_cif, times_single_cif)
  # Expected: yleft=0, linear between points, yright=max(probs)=0.5
  # t=0 -> 0
  # t=1 -> 0.1
  # t=2 -> linear(1,0.1) to (3,0.3) -> 0.1 + (0.3-0.1)/(3-1)*(2-1) = 0.1 + 0.1 = 0.2
  # t=3 -> 0.3
  # t=4 -> linear(3,0.3) to (5,0.4) -> 0.3 + (0.4-0.3)/(5-3)*(4-3) = 0.3 + 0.05 = 0.35
  # t=5 -> 0.4
  # t=6 -> linear(5,0.4) to (8,0.5) -> 0.4 + (0.5-0.4)/(8-5)*(6-5) = 0.4 + 0.0333 = 0.4333
  # t=7 -> linear(5,0.4) to (8,0.5) -> 0.4 + (0.5-0.4)/(8-5)*(7-5) = 0.4 + 0.0667 = 0.4667
  # t=8 -> 0.5
  # t=9 -> 0.5 (yright)
  expected_probs <- c(0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.4 + 1/30, 0.4 + 2/30, 0.5, 0.5)
  expect_equal(interpolated_probs, expected_probs, tolerance = 1e-6)
})

test_that("cifInterpolator handles yleft and yright", {
  # Test time before first time point
  expect_equal(cifInterpolator(0.5, probs_single_cif, times_single_cif), 0) # yleft = 0
  # Test time after last time point
  expect_equal(cifInterpolator(10, probs_single_cif, times_single_cif), max(probs_single_cif)) # yright = max
})

test_that("cifInterpolator handles unsorted input", {
  times_unsorted <- c(5, 1, 8, 3)
  probs_unsorted <- c(0.4, 0.1, 0.5, 0.3) # Corresponds to 1=0.1, 3=0.3, 5=0.4, 8=0.5
  interpolated_probs <- cifInterpolator(new_times_single_cif, probs_unsorted, times_unsorted)
  expected_probs <- c(0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.4 + 1/30, 0.4 + 2/30, 0.5, 0.5)
  expect_equal(interpolated_probs, expected_probs, tolerance = 1e-6)
})

test_that("cifInterpolator handles single time/prob input", {
   # Linear interpolation needs at least two points implicitly (t=0, p=0)
   interpolated <- cifInterpolator(c(0, 2.5, 5, 7.5), 0.6, 5)
   # Interpolates between (0,0) and (5, 0.6). yright=0.6
   # t=0 -> 0
   # t=2.5 -> 0.3
   # t=5 -> 0.6
   # t=7.5 -> 0.6
   expect_equal(interpolated, c(0, 0.3, 0.6, 0.6))
})

test_that("cifInterpolator handles NA in probs (uses max of non-NA for yright)", {
   probs_na <- c(0.1, NA, 0.4, 0.5)
   interpolated <- cifInterpolator(c(9, 10), probs_na, times_single_cif)
   expect_equal(interpolated, c(0.5, 0.5)) # yright should be max(0.1, 0.4, 0.5) = 0.5
})


# --- Tests for cifMatInterpolaltor ---

test_that("cifMatInterpolaltor interpolates matrix correctly", {
  interpolated_mat <- cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif)

  # Check dimensions: rows = newtimes, cols = observations
  expect_equal(nrow(interpolated_mat), length(new_times_mat_cif))
  expect_equal(ncol(interpolated_mat), nrow(cif_mat1))

  # Check interpolation for subject 1 (probs: 0.1@t=2, 0.3@t=5, 0.5@t=10)
  # Expected at new_times_mat_cif = c(0, 1, 2, 3, 4, 5, 8, 10, 12)
  # t=0->0; t=1->linear(0,0 to 2,0.1)=0.05; t=2->0.1; t=3->linear(2,0.1 to 5,0.3)=0.1+0.2/3*1=0.1667;
  # t=4->0.1+0.2/3*2=0.2333; t=5->0.3; t=8->linear(5,0.3 to 10,0.5)=0.3+0.2/5*3=0.42;
  # t=10->0.5; t=12->0.5(yright)
  expected_subj1 <- c(0, 0.05, 0.1, 0.1 + 1/15, 0.1 + 2/15, 0.3, 0.3 + 3/25, 0.5, 0.5)
  expect_equal(interpolated_mat[, 1], expected_subj1, tolerance = 1e-6)

  # Check interpolation for subject 2 (probs: 0.2@t=2, 0.4@t=5, 0.7@t=10)
  # t=0->0; t=1->linear(0,0 to 2,0.2)=0.1; t=2->0.2; t=3->linear(2,0.2 to 5,0.4)=0.2+0.2/3*1=0.2667;
  # t=4->0.2+0.2/3*2=0.3333; t=5->0.4; t=8->linear(5,0.4 to 10,0.7)=0.4+0.3/5*3=0.58;
  # t=10->0.7; t=12->0.7(yright)
  expected_subj2 <- c(0, 0.1, 0.2, 0.2 + 1/15, 0.2 + 2/15, 0.4, 0.4 + 9/50, 0.7, 0.7)
  expect_equal(interpolated_mat[, 2], expected_subj2, tolerance = 1e-6)
})

test_that("cifMatInterpolaltor handles time 0 correctly", {
  # Case 1: time 0 not in input times (tested above)
  interpolated_mat <- cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif)
  expect_equal(interpolated_mat[1, ], c(0, 0)) # First row (time 0) should be 0

  # Case 2: time 0 is in input times
  times_with_zero <- c(0, times_mat_cif)
  cif_with_zero <- cbind(rep(0, nrow(cif_mat1)), cif_mat1)
  interpolated_mat_zero <- cifMatInterpolaltor(cif_with_zero, times_with_zero, new_times_mat_cif)
  expected_subj1 <- c(0, 0.05, 0.1, 0.1 + 1/15, 0.1 + 2/15, 0.3, 0.3 + 3/25, 0.5, 0.5)
  expected_subj2 <- c(0, 0.1, 0.2, 0.2 + 1/15, 0.2 + 2/15, 0.4, 0.4 + 9/50, 0.7, 0.7)
  expect_equal(interpolated_mat_zero[, 1], expected_subj1, tolerance = 1e-6)
  expect_equal(interpolated_mat_zero[, 2], expected_subj2, tolerance = 1e-6)
})

test_that("cifMatInterpolaltor enforces monotonicity", {
  # Create a matrix where interpolation might initially decrease
  probs_nonmono <- matrix(c(0.1, 0.3, 0.25), nrow = 1) # Decreases from 0.3 to 0.25
  times_nonmono <- c(2, 5, 10)
  new_times_nonmono <- c(1, 3, 6, 11)
  # Initial interpolation: 0.05, 0.1667, 0.3 + (0.25-0.3)/5*1 = 0.29, 0.25(yright)
  # Monotonicity correction should make it: 0.05, 0.1667, 0.29, 0.29 (using cummax)
  interpolated_mat <- cifMatInterpolaltor(probs_nonmono, times_nonmono, new_times_nonmono)
  # Re-calc: t=1->0.05; t=3->0.1667; t=6->0.29; t=11->0.25(yright)
  # cummax: 0.05, 0.1667, 0.29, 0.29
  # Corrected: 0.05, 0.1667, 0.29, 0.29
  expect_equal(as.vector(interpolated_mat), c(0.05, 0.1666667, 0.29, 0.29), tolerance=1e-6)
})

test_that("cifMatInterpolaltor handles single new time", {
   interpolated_mat <- cifMatInterpolaltor(cif_mat1, times_mat_cif, newtimes = 4)
   expect_true(is.matrix(interpolated_mat))
   expect_equal(nrow(interpolated_mat), 1)
   expect_equal(ncol(interpolated_mat), nrow(cif_mat1))
   # Expected at t=4: Subj1=0.2333, Subj2=0.3333
   expect_equal(as.vector(interpolated_mat), c(0.1 + 2/15, 0.2 + 2/15), tolerance=1e-6)
})


# --- Tests for cifMatListAveraging ---

test_that("cifMatListAveraging averages correctly on CumHaz scale", {
  # Use interpolated matrices (newtimes x observations)
  list_mats <- list(cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif),
                    cifMatInterpolaltor(cif_mat2, times_mat_cif, new_times_mat_cif))

  averaged_mat <- cifMatListAveraging(list_mats, type = "CumHaz")

  # Check dimensions
  expect_equal(dim(averaged_mat), dim(list_mats[[1]]))

  # Check a specific value (e.g., Subject 1 at newtime = 8)
  # Mat1: S(8) = 0.42
  # Mat2: S(8) = linear(5,0.25 to 10,0.45) = 0.25 + 0.2/5*3 = 0.37
  # H1 = -log(1 - 0.42) = -log(0.58) ~= 0.5447
  # H2 = -log(1 - 0.37) = -log(0.63) ~= 0.4620
  # Mean H = (0.5447 + 0.4620) / 2 ~= 0.50335
  # Averaged P = 1 - exp(-0.50335) ~= 1 - 0.6045 = 0.3955
  subj1_idx <- 1
  time8_idx <- which(new_times_mat_cif == 8)
  val1 <- list_mats[[1]][time8_idx, subj1_idx]
  val2 <- list_mats[[2]][time8_idx, subj1_idx]
  expect_equal(averaged_mat[time8_idx, subj1_idx], 1 - exp(-mean(c(-log(1-val1), -log(1-val2)))), tolerance = 1e-6)

})

test_that("cifMatListAveraging averages correctly on prob scale", {
  list_mats <- list(cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif),
                    cifMatInterpolaltor(cif_mat2, times_mat_cif, new_times_mat_cif))

  averaged_mat <- cifMatListAveraging(list_mats, type = "prob")

  # Check dimensions
  expect_equal(dim(averaged_mat), dim(list_mats[[1]]))

  # Check a specific value (e.g., Subject 1 at newtime = 8)
  # Mat1: S(8) = 0.42
  # Mat2: S(8) = 0.37
  # Mean P = (0.42 + 0.37) / 2 = 0.395
  subj1_idx <- 1
  time8_idx <- which(new_times_mat_cif == 8)
  val1 <- list_mats[[1]][time8_idx, subj1_idx]
  val2 <- list_mats[[2]][time8_idx, subj1_idx]
  expect_equal(averaged_mat[time8_idx, subj1_idx], mean(c(val1, val2)), tolerance = 1e-6)
})


test_that("cifMatListAveraging handles list with one matrix", {
  list_one <- list(cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif))
  averaged_mat_cumhaz <- cifMatListAveraging(list_one, type = "CumHaz")
  averaged_mat_prob <- cifMatListAveraging(list_one, type = "prob")
  expect_equal(averaged_mat_cumhaz, list_one[[1]])
  expect_equal(averaged_mat_prob, list_one[[1]])
})

test_that("cifMatListAveraging handles empty list", {
  expect_null(cifMatListAveraging(list()))
})

test_that("cifMatListAveraging handles inconsistent dimensions", {
  mat_wrong_dim <- matrix(1:6, nrow=3) # Different rows
  list_inconsistent <- list(cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif),
                            mat_wrong_dim)
  expect_error(cifMatListAveraging(list_inconsistent), "must have the same dimensions")
})

test_that("cifMatListAveraging handles NAs (na.rm=TRUE behavior)", {
   mat1 <- matrix(c(0.1, 0.3, NA, 0.4), nrow=2) # times x obs, [1,2] = NA
   mat2 <- matrix(c(0.2, NA, NA, 0.5), nrow=2)  # [1,2] = NA, [2,1] = NA
   mat3 <- matrix(c(NA, NA, NA, 0.6), nrow=2)   # [1,1] = NA, [1,2] = NA, [2,1] = NA
   list_na <- list(mat1, mat2, mat3)

   # CumHaz
   averaged_ch <- cifMatListAveraging(list_na, type="CumHaz")
   # [1, 1]: 1-exp(-mean(c(-log(1-0.1), -log(1-0.2)))) = 0.1514719
   expect_equal(averaged_ch[1,1], 1 - exp(-mean(c(-log(0.9), -log(0.8)))))
   # [1, 2]: All inputs are NA (NA, NA, NA) -> should be NA
   expect_true(is.na(averaged_ch[1,2]))
   # [2, 1]: 1-exp(-mean(c(-log(1-0.3)))) = 0.3 (only mat1[2,1]=0.3 is non-NA)
   expect_equal(averaged_ch[2,1], 1 - exp(-mean(c(-log(0.7)))))
   # [2, 2]: 1-exp(-mean(c(-log(1-0.4), -log(1-0.5), -log(1-0.6)))) = 0.5067576
   expect_equal(averaged_ch[2,2], 1 - exp(-mean(c(-log(0.6), -log(0.5), -log(0.4)))))

   # Prob
   averaged_p <- cifMatListAveraging(list_na, type="prob")
   expect_equal(averaged_p[1,1], mean(c(0.1, 0.2), na.rm=TRUE))  # mean(0.1, 0.2) = 0.15
   expect_true(is.na(averaged_p[1,2])) # mean(NA, NA, NA) with na.rm=TRUE gives NaN, should be NA
   expect_equal(averaged_p[2,1], mean(c(0.3), na.rm=TRUE))       # only 0.3 is non-NA = 0.3
   expect_equal(averaged_p[2,2], mean(c(0.4, 0.5, 0.6), na.rm=TRUE)) # mean(0.4, 0.5, 0.6) = 0.5
})

test_that("cifMatListAveraging requires valid type", {
    list_mats <- list(cifMatInterpolaltor(cif_mat1, times_mat_cif, new_times_mat_cif))
    expect_error(cifMatListAveraging(list_mats, type="invalid"), "Type must be either 'CumHaz' or 'prob'")
})
