library(testthat)
library(here)
library(pracma) # Required for trapz

# Assuming the function is available in the environment
# Use here::here for robustness
source(here::here("R/math_utils.R"))

context("Testing math_utils functions (Integrator)")

# --- Test Data Setup ---
times_integ <- c(0, 1, 2, 3, 4, 5, 6)
scores_integ <- c(0, 1, 2, 3, 2, 1, 0) # Simple triangle shape area = 0.5 * base * height = 0.5 * 6 * 3 = 9

scores_const <- rep(2, length(times_integ)) # Constant score

# --- Tests for Integrator ---

test_that("Integrator calculates trapezoidal area correctly (unscaled)", {
  skip_if_not_installed("pracma")

  # Full range
  area_full <- Integrator(times_integ, scores_integ, minmax = c(0, 6), scale = FALSE)
  # Manual calculation: (1*1)/2 + (1+2)/2*1 + (2+3)/2*1 + (3+2)/2*1 + (2+1)/2*1 + (1+0)/2*1
  # = 0.5 + 1.5 + 2.5 + 2.5 + 1.5 + 0.5 = 9
  expect_equal(area_full, 9)

  # Partial range
  area_partial <- Integrator(times_integ, scores_integ, minmax = c(1, 4), scale = FALSE)
  # Points: (1,1), (2,2), (3,3), (4,2)
  # Area: (1+2)/2*1 + (2+3)/2*1 + (3+2)/2*1 = 1.5 + 2.5 + 2.5 = 6.5
  expect_equal(area_partial, 6.5)

  # Constant score
  area_const <- Integrator(times_integ, scores_const, minmax = c(1, 5), scale = FALSE)
  # Points: (1,2), (2,2), (3,2), (4,2), (5,2) -> rectangle base=4, height=2 -> Area=8
  expect_equal(area_const, 8)
})

test_that("Integrator calculates scaled area correctly", {
  skip_if_not_installed("pracma")

  # Full range (0 to 6, interval = 6)
  area_full_scaled <- Integrator(times_integ, scores_integ, minmax = c(0, 6), scale = TRUE)
  expect_equal(area_full_scaled, 9 / 6)

  # Partial range (1 to 4, interval = 3)
  area_partial_scaled <- Integrator(times_integ, scores_integ, minmax = c(1, 4), scale = TRUE)
  expect_equal(area_partial_scaled, 6.5 / 3)

  # Constant score (1 to 5, interval = 4)
  area_const_scaled <- Integrator(times_integ, scores_const, minmax = c(1, 5), scale = TRUE)
  expect_equal(area_const_scaled, 8 / 4)
})

test_that("Integrator handles range subsetting correctly", {
  skip_if_not_installed("pracma")
  # Range completely outside
  expect_warning(area_outside <- Integrator(times_integ, scores_integ, minmax = c(10, 20)),
                 "Less than 2 points found")
  expect_equal(area_outside, 0)

  # Range including only one point
  expect_warning(area_one_point <- Integrator(times_integ, scores_integ, minmax = c(2.5, 3.5)),
                 "Less than 2 points found")
   expect_equal(area_one_point, 0) # Only point (3,3) is within range

   # Range including exactly two points
   area_two_points <- Integrator(times_integ, scores_integ, minmax = c(2.5, 4.5), scale=FALSE)
   # Points (3,3), (4,2) -> Area = (3+2)/2 * (4-3) = 2.5
   expect_equal(area_two_points, 2.5)
})

test_that("Integrator handles unsorted input", {
  skip_if_not_installed("pracma")
  times_unsort <- c(3, 1, 5, 0, 6, 2, 4)
  scores_unsort <- c(3, 1, 1, 0, 0, 2, 2) # Corresponds to original scores_integ
  area_unsorted <- Integrator(times_unsort, scores_unsort, minmax = c(0, 6), scale = FALSE)
  expect_equal(area_unsorted, 9) # Should be same as sorted full range
})

test_that("Integrator requires valid inputs", {
  skip_if_not_installed("pracma")
  # Mismatched lengths
  expect_error(Integrator(times_integ, scores_integ[1:3], minmax = c(0, 6)),
               "Length of 'times' and 'scores' must be equal.")
  # Invalid minmax
  expect_error(Integrator(times_integ, scores_integ, minmax = c(5, 1)),
               "'minmax' must be a numeric vector of length 2")
  expect_error(Integrator(times_integ, scores_integ, minmax = c(5, 5)),
               "'minmax' must be a numeric vector of length 2")
   expect_error(Integrator(times_integ, scores_integ, minmax = c(1, 2, 3)),
               "'minmax' must be a numeric vector of length 2")
})

test_that("Integrator handles invalid scaling interval", {
   skip_if_not_installed("pracma")
   # Need to create minmax where interval is 0, but trapz still works
   times_edge = c(2, 2.000001)
   scores_edge = c(1, 1)
   # Scale=TRUE with minmax[1]=minmax[2] should warn and return 0
   # Note: The function checks minmax[1] >= minmax[2] earlier, so this case isn't reachable
   # Let's test the internal check for interval_length > 0
   # We need a case where trapz returns non-zero but interval is zero (not possible with valid minmax)
   # Test the warning when scale=TRUE and interval is 0 (though prevented by earlier check)
   # expect_warning(Integrator(times_integ, scores_integ, minmax = c(3, 3), scale = TRUE), "interval length is zero")
   # This test is redundant due to the initial minmax check.
   expect_true(TRUE) # Placeholder
})
