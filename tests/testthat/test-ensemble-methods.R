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
