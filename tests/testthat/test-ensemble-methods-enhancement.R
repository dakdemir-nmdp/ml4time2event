# Test-Driven Development for Additional Ensemble Methods
# Testing new ensemble methods: "median", "min", "max", "stacking"

test_that("median ensemble method works for survival models", {
  skip_if_not_installed("survival")

  # Create mock prediction data
  pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  pred2 <- matrix(c(0.85, 0.75, 0.65, 0.55), ncol = 4)
  pred3 <- matrix(c(0.95, 0.85, 0.75, 0.65), ncol = 4)
  times <- c(1, 2, 3, 4)

  # Test median ensemble
  result <- EnsemblePredictions(
    list(pred1, pred2, pred3),
    times = times,
    ensemble_method = "median"
  )

  expect_equal(nrow(result$Probs), 1)
  expect_equal(ncol(result$Probs), 4)
  expect_equal(result$Times, times)

  # Check that median is calculated correctly
  expected_median <- apply(array(c(pred1, pred2, pred3), dim = c(1, 4, 3)), 
                          c(1, 2), median)
  expect_equal(as.numeric(result$Probs), as.numeric(expected_median))
})

test_that("min ensemble method works for survival models", {
  skip_if_not_installed("survival")

  # Create mock prediction data
  pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  pred2 <- matrix(c(0.85, 0.75, 0.65, 0.55), ncol = 4)
  pred3 <- matrix(c(0.95, 0.85, 0.75, 0.65), ncol = 4)
  times <- c(1, 2, 3, 4)

  # Test min ensemble (conservative predictions)
  result <- EnsemblePredictions(
    list(pred1, pred2, pred3),
    times = times,
    ensemble_method = "min"
  )

  expect_equal(nrow(result$Probs), 1)
  expect_equal(ncol(result$Probs), 4)
  expect_equal(result$Times, times)

  # Check that min is calculated correctly
  expected_min <- apply(array(c(pred1, pred2, pred3), dim = c(1, 4, 3)), 
                        c(1, 2), min)
  expect_equal(as.numeric(result$Probs), as.numeric(expected_min))
})

test_that("max ensemble method works for survival models", {
  skip_if_not_installed("survival")

  # Create mock prediction data
  pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  pred2 <- matrix(c(0.85, 0.75, 0.65, 0.55), ncol = 4)
  pred3 <- matrix(c(0.95, 0.85, 0.75, 0.65), ncol = 4)
  times <- c(1, 2, 3, 4)

  # Test max ensemble (optimistic predictions)
  result <- EnsemblePredictions(
    list(pred1, pred2, pred3),
    times = times,
    ensemble_method = "max"
  )

  expect_equal(nrow(result$Probs), 1)
  expect_equal(ncol(result$Probs), 4)
  expect_equal(result$Times, times)

  # Check that max is calculated correctly
  expected_max <- apply(array(c(pred1, pred2, pred3), dim = c(1, 4, 3)), 
                        c(1, 2), max)
  expect_equal(as.numeric(result$Probs), as.numeric(expected_max))
})

test_that("stacking ensemble method works for survival models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")

  # Create mock prediction data with training and validation sets
  # Training predictions
  train_pred1 <- matrix(runif(20, 0.3, 0.9), ncol = 4)
  train_pred2 <- matrix(runif(20, 0.2, 0.8), ncol = 4)
  train_pred3 <- matrix(runif(20, 0.4, 0.95), ncol = 4)

  # Validation predictions
  val_pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  val_pred2 <- matrix(c(0.85, 0.75, 0.65, 0.55), ncol = 4)
  val_pred3 <- matrix(c(0.95, 0.85, 0.75, 0.65), ncol = 4)

  times <- c(1, 2, 3, 4)

  # Mock training outcomes (survival probabilities at evaluation times)
  train_outcomes <- matrix(runif(20, 0.2, 0.9), ncol = 4)

  # Ensure training predictions are named
  sl_training_predictions <- list(
    pred1 = train_pred1,
    pred2 = train_pred2,
    pred3 = train_pred3
  )

  # Test stacking ensemble with meta-learner
  result <- EnsemblePredictions(
    list(val_pred1, val_pred2, val_pred3),
    times = times,
    ensemble_method = "stacking",
    sl_training_predictions = sl_training_predictions,
    sl_actual = train_outcomes
  )

  expect_equal(nrow(result$Probs), 1)
  expect_equal(ncol(result$Probs), 4)
  expect_equal(result$Times, times)

  # Stacking results should be different from simple averaging
  avg_result <- EnsemblePredictions(
    list(val_pred1, val_pred2, val_pred3),
    times = times,
    ensemble_method = "average"
  )

  expect_false(identical(result$Probs, avg_result$Probs))
})

test_that("geometric_mean ensemble method works for survival models", {
  skip_if_not_installed("survival")

  # Create mock prediction data (avoiding zeros for geometric mean)
  pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  pred2 <- matrix(c(0.85, 0.75, 0.65, 0.55), ncol = 4)
  pred3 <- matrix(c(0.95, 0.85, 0.75, 0.65), ncol = 4)
  times <- c(1, 2, 3, 4)

  # Test geometric mean ensemble
  result <- EnsemblePredictions(
    list(pred1, pred2, pred3),
    times = times,
    ensemble_method = "geometric_mean"
  )

  expect_equal(nrow(result$Probs), 1)
  expect_equal(ncol(result$Probs), 4)
  expect_equal(result$Times, times)

  # Check that geometric mean is calculated correctly
  expected_geom <- apply(array(c(pred1, pred2, pred3), dim = c(1, 4, 3)), 
                        c(1, 2), function(x) exp(mean(log(x))))
  expect_equal(as.numeric(result$Probs), as.numeric(expected_geom), tolerance = 1e-10)
})

test_that("ensemble method validation works correctly", {
  pred1 <- matrix(c(0.9, 0.8, 0.7, 0.6), ncol = 4)
  times <- c(1, 2, 3, 4)

  # Test invalid ensemble method
  expect_error(
    EnsemblePredictions(list(pred1), times = times, ensemble_method = "invalid_method"),
    "ensemble_method must be one of: 'average', 'weighted', 'super_learner', 'median', 'min', 'max', 'geometric_mean', 'stacking'"
  )

  # Test that all new methods are recognized
  valid_methods <- c("average", "weighted", "super_learner", "median", "min", "max", 
                    "stacking", "geometric_mean")

  for (method in valid_methods) {
    if (method %in% c("weighted", "stacking")) {
      # Skip methods that need additional parameters for basic validation
      next
    }

    expect_no_error({
      result <- EnsemblePredictions(list(pred1), times = times, ensemble_method = method)
    })
  }
})