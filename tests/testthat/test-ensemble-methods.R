# ============================================================================
# Test Suite for Advanced Ensemble Methods (Phase 3)
# Focus on weighted averaging helpers and EnsemblePredictions interface
# ============================================================================

library(testthat)

context("Advanced ensemble methods")

# ----------------------------------------------------------------------------
# Helper for constructing simple prediction matrices
# ----------------------------------------------------------------------------
make_surv_matrix <- function(values) {
  matrix(values, nrow = length(values), ncol = 1)
}

make_cif_matrix <- function(values) {
  matrix(values, nrow = length(values), ncol = 1)
}

# ----------------------------------------------------------------------------
# Weighted averaging helpers
# ----------------------------------------------------------------------------

test_that("survprobMatWeightedAveraging matches manual hazard weighting", {
  preds <- list(
    ModelA = make_surv_matrix(c(0.9, 0.8)),
    ModelB = make_surv_matrix(c(0.5, 0.4))
  )

  weights <- c(ModelA = 0.7, ModelB = 0.3)

  result <- ml4time2event:::survprobMatWeightedAveraging(preds, weights)

  manual <- exp(-(
    (-log(preds$ModelA + 1e-10)) * weights["ModelA"] +
      (-log(preds$ModelB + 1e-10)) * weights["ModelB"]
  ))
  manual <- matrix(manual, nrow = nrow(preds$ModelA), ncol = ncol(preds$ModelA))

  expect_true(is.matrix(result))
  expect_equal(result, manual, tolerance = 1e-8)
})

test_that("cifMatWeightedAveraging drops mismatched dimensions with warning", {
  preds <- list(
    ModelA = make_cif_matrix(c(0.3, 0.5)),
    ModelB = make_cif_matrix(c(0.4, 0.6)),
    ModelOdd = matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 2)
  )

  weights <- c(ModelA = 0.4, ModelB = 0.6, ModelOdd = 0.2)

  expect_warning(
    result <- ml4time2event:::cifMatWeightedAveraging(preds, weights, type = "prob"),
    "dimension mismatch"
  )

  result <- matrix(result, nrow = nrow(preds$ModelA), ncol = ncol(preds$ModelA))

  expect_true(all(result >= 0 & result <= 1))
  expect_equal(dim(result), dim(preds$ModelA))
})

# ----------------------------------------------------------------------------
# EnsemblePredictions interface
# ----------------------------------------------------------------------------

test_that("EnsemblePredictions weighted requires named weights", {
  preds <- list(
    ModelA = make_surv_matrix(c(0.9, 0.8)),
    ModelB = make_surv_matrix(c(0.7, 0.6))
  )

  expect_error(
    EnsemblePredictions(preds, ensemble_method = "weighted", model_weights = c(0.6, 0.4), type = "survival", times = c(1, 2)),
    "named"
  )
})

test_that("EnsemblePredictions super learner requires training data or weights", {
  preds <- list(
    ModelA = make_surv_matrix(c(0.9, 0.8)),
    ModelB = make_surv_matrix(c(0.7, 0.6))
  )

  expect_error(
    EnsemblePredictions(preds, ensemble_method = "super_learner", type = "survival", times = c(1, 2)),
    "requires either"
  )
})

test_that("EnsemblePredictions super learner computes optimized weights", {
  train_preds <- list(
    ModelA = make_surv_matrix(c(0.95, 0.85)),
    ModelB = make_surv_matrix(c(0.6, 0.5))
  )

  actual <- make_surv_matrix(c(0.92, 0.82))

  new_preds <- list(
    ModelA = make_surv_matrix(c(0.93, 0.83)),
    ModelB = make_surv_matrix(c(0.65, 0.55))
  )

  result <- EnsemblePredictions(
    new_preds,
    ensemble_method = "super_learner",
    type = "survival",
    sl_training_predictions = train_preds,
    sl_actual = actual,
    sl_loss = "mse",
    times = c(1, 2)
  )

  weights_attr <- attr(result, "sl_weights")
  expect_true(is.numeric(weights_attr))
  expect_equal(names(weights_attr), names(new_preds))
  expect_equal(sum(weights_attr), 1, tolerance = 1e-6)

  expected_weights <- optimizeSuperLearnerWeights(train_preds, actual, loss_type = "mse")
  expect_equal(weights_attr, expected_weights[names(new_preds)], tolerance = 1e-6)

  expect_true(all(result >= 0 & result <= 1))
  expect_equal(dim(result), dim(new_preds$ModelA))
})

test_that("EnsemblePredictions super learner accepts pre-computed weights", {
  preds <- list(
    ModelA = make_surv_matrix(c(0.9, 0.8)),
    ModelB = make_surv_matrix(c(0.7, 0.6))
  )

  weights <- c(ModelA = 0.2, ModelB = 0.8)

  result <- EnsemblePredictions(
    preds,
    ensemble_method = "super_learner",
    type = "survival",
    sl_weights = weights,
    times = c(1, 2)
  )

  weights_attr <- attr(result, "sl_weights")
  expect_equal(weights_attr, weights / sum(weights))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(dim(result), dim(preds$ModelA))
})

# ----------------------------------------------------------------------------
# CR Ensemble Fix Tests
# ----------------------------------------------------------------------------

test_that("ensemble predictions work with synthetic data (CR fix)", {
  # Create synthetic test data
  set.seed(123)
  n_obs <- 20
  n_times <- 5
  eval_times <- seq(0.1, 1, by = 0.2)

  # Create test data frame
  test_data <- data.frame(
    x1 = rnorm(n_obs),
    x2 = rnorm(n_obs),
    group = sample(c("low", "high"), n_obs, replace = TRUE)
  )

  # Create synthetic prediction matrices for two models
  # Model 1: High risk for 'low' group, low risk for 'high' group
  pred_mat1 <- matrix(0, nrow = n_times, ncol = n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i / n_times
    pred_mat1[i, test_data$group == "low"] <- 0.8 * time_factor
    pred_mat1[i, test_data$group == "high"] <- 0.2 * time_factor
  }

  # Model 2: Similar pattern but more extreme
  pred_mat2 <- matrix(0, nrow = n_times, ncol = n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i / n_times
    pred_mat2[i, test_data$group == "low"] <- 0.9 * time_factor
    pred_mat2[i, test_data$group == "high"] <- 0.1 * time_factor
  }

  # Create prediction list
  model_predictions <- list(
    "Model1" = pred_mat1,
    "Model2" = pred_mat2
  )

  # Test simple averaging
  ensemble_pred <- cifMatListAveraging(model_predictions, type = "prob")

  # Test ensemble structure
  expect_true(is.matrix(ensemble_pred))
  expect_equal(nrow(ensemble_pred), n_times)
  expect_equal(ncol(ensemble_pred), n_obs)

  # Verify ensemble predictions average the two models
  # For low group, should be around 0.85
  # For high group, should be around 0.15
  low_indices <- which(test_data$group == "low")
  high_indices <- which(test_data$group == "high")

  # Check final time point predictions
  expect_true(all(abs(ensemble_pred[n_times, low_indices] - 0.85) < 0.01))
  expect_true(all(abs(ensemble_pred[n_times, high_indices] - 0.15) < 0.01))

  # Test weighted averaging
  weights <- c(Model1 = 0.75, Model2 = 0.25)
  ensemble_weighted <- cifMatWeightedAveraging(model_predictions, weights, type = "prob")

  # Test weighted ensemble structure
  expect_true(is.matrix(ensemble_weighted))
  expect_equal(dim(ensemble_weighted), dim(pred_mat1))

  # Verify weighted predictions match the weights
  # For low group: 0.75*0.8 + 0.25*0.9 = 0.6 + 0.225 = 0.825
  # For high group: 0.75*0.2 + 0.25*0.1 = 0.15 + 0.025 = 0.175
  expect_true(all(abs(ensemble_weighted[n_times, low_indices] - 0.825) < 0.01))
  expect_true(all(abs(ensemble_weighted[n_times, high_indices] - 0.175) < 0.01))
})

# ----------------------------------------------------------------------------
# Enhanced Ensemble Methods Tests (Median, Min, Max, Stacking, Geometric)
# ----------------------------------------------------------------------------

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
  expected_median <- apply(
    array(c(pred1, pred2, pred3), dim = c(1, 4, 3)),
    c(1, 2), median
  )
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
  expected_min <- apply(
    array(c(pred1, pred2, pred3), dim = c(1, 4, 3)),
    c(1, 2), min
  )
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
  expected_max <- apply(
    array(c(pred1, pred2, pred3), dim = c(1, 4, 3)),
    c(1, 2), max
  )
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
  expected_geom <- apply(
    array(c(pred1, pred2, pred3), dim = c(1, 4, 3)),
    c(1, 2), function(x) exp(mean(log(x)))
  )
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
  valid_methods <- c(
    "average", "weighted", "super_learner", "median", "min", "max",
    "stacking", "geometric_mean"
  )

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
