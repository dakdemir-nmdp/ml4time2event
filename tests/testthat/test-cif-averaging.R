test_that("cifMatListAveraging propagates NAs for 'prob' type", {
  m1 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  m2 <- matrix(c(0.2, NA, 0.4, 0.5), nrow = 2)
  m3 <- matrix(c(0.15, 0.25, 0.35, 0.45), nrow = 2)
  
  list_mats <- list(m1, m2, m3)
  
  # With na.rm = TRUE, the NA is ignored. The result for that element would be (0.2 + 0.25) / 2 = 0.225
  # With na.rm = FALSE, the result for that element should be NA.
  
  # Calculate expected result with NA propagation
  expected_avg <- (m1 + m2 + m3) / 3
  
  # Run the function
  actual_avg <- cifMatListAveraging(list_mats, type = "prob")
  
  # The test will fail if the NA is ignored
  expect_equal(actual_avg, expected_avg)
})

test_that("cifMatListAveraging propagates NAs for 'CumHaz' type", {
  m1 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  m2 <- matrix(c(0.2, NA, 0.4, 0.5), nrow = 2)
  m3 <- matrix(c(0.15, 0.25, 0.35, 0.45), nrow = 2)
  
  list_mats <- list(m1, m2, m3)
  
  # With na.rm = TRUE, the NA is ignored.
  # With na.rm = FALSE, the result should be NA where m2 has an NA.
  
  # Calculate expected result with NA propagation
  h1 <- -log(1 - m1 + 1e-10)
  h2 <- -log(1 - m2 + 1e-10) # will have NA
  h3 <- -log(1 - m3 + 1e-10)
  
  h_avg_array <- array(c(h1, h2, h3), dim = c(2, 2, 3))
  mean_h <- apply(h_avg_array, c(1, 2), mean, na.rm = FALSE)
  expected_cif <- 1 - exp(-mean_h)
  
  # Run the function
  actual_cif <- cifMatListAveraging(list_mats, type = "CumHaz")
  
  # The test will fail if the NA is ignored
  expect_equal(actual_cif, expected_cif)
})
