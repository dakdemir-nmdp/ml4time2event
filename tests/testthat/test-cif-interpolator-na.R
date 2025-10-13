test_that("cifInterpolator propagates NAs correctly", {
  times <- c(10, 20, 30, 40)
  probs <- c(0.1, NA, 0.3, 0.4)
  newtimes <- c(15, 25, 35)
  
  # Expected behavior:
  # 15 is between 10 and 20 (NA), so result should be NA.
  # 25 is between 20 (NA) and 30, so result should be NA.
  # 35 is between 30 and 40, so result should be 0.35.
  
  expected_result <- c(NA, NA, 0.35)
  
  actual_result <- cifInterpolator(newtimes, probs, times)
  
  expect_equal(actual_result, expected_result)
})

test_that("cifInterpolator handles NA at the end", {
  times <- c(10, 20, 30, 40)
  probs <- c(0.1, 0.2, 0.3, NA)
  newtimes <- c(15, 25, 35, 45)
  
  # Expected behavior:
  # 15 -> 0.15
  # 25 -> 0.25
  # 35 is between 30 and 40 (NA), so result should be NA.
  # 45 is > 40 (NA), and yright will be NA, so result should be NA.
  
  expected_result <- c(0.15, 0.25, NA, NA)
  
  actual_result <- cifInterpolator(newtimes, probs, times)
  
  expect_equal(actual_result, expected_result)
})
