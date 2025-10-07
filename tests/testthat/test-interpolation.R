# ==============================================================================
# Test Suite for Survival Interpolation Utilities
# ==============================================================================

library(testthat)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(123)
n_train <- 100

train_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = rbinom(n_train, 1, 0.7),
  x1 = rnorm(n_train),
  x2 = rnorm(n_train)
)

test_data <- data.frame(
  time = rexp(20, rate = 0.1),
  event = rbinom(20, 1, 0.7),
  x1 = rnorm(20),
  x2 = rnorm(20)
)

# ==============================================================================
# Tests for survprobMatInterpolator
# ==============================================================================

test_that("survprobMatInterpolator handles standard input", {
  # Fit a simple model
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Get predictions at model's natural times
  preds1 <- Predict_SurvModel_Cox(model, test_data)

  # Interpolate to different times
  custom_times <- c(5, 10, 15, 20, 25)
  preds2 <- Predict_SurvModel_Cox(model, test_data, newtimes = custom_times)

  # Check dimensions
  expect_equal(nrow(preds2$Probs), length(custom_times))
  expect_equal(ncol(preds2$Probs), nrow(test_data))
  expect_equal(preds2$Times, custom_times)
})

test_that("survprobMatInterpolator maintains monotonicity", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Use many time points including some very close together
  dense_times <- sort(c(seq(0, 50, by = 0.5), seq(10, 15, by = 0.1)))
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = dense_times)

  # Check each observation's survival curve is non-increasing
  for (i in seq_len(ncol(preds$Probs))) {
    surv_curve <- preds$Probs[, i]
    diffs <- diff(surv_curve)
    expect_true(all(diffs <= 1e-10),
                info = paste("Observation", i, "not monotonic"))
  }
})

test_that("survprobMatInterpolator handles time 0", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Request times starting from 0
  times_with_0 <- c(0, 5, 10, 20)
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = times_with_0)

  # Time 0 should be included
  expect_true(0 %in% preds$Times)

  # Survival at time 0 should be 1
  time_0_idx <- which(preds$Times == 0)
  expect_true(all(abs(preds$Probs[time_0_idx, ] - 1) < 1e-6))
})

test_that("survprobMatInterpolator handles extrapolation", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Request times beyond observed data
  max_observed <- model$time_range[2]
  extrap_times <- c(0, 5, 10, max_observed * 1.5, max_observed * 2)

  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = extrap_times)

  # Should not produce NAs
  expect_true(all(!is.na(preds$Probs)))

  # Extrapolated values should be constant (last observed survival)
  last_obs_idx <- which(extrap_times == max_observed * 1.5)
  beyond_idx <- which(extrap_times == max_observed * 2)

  # Values beyond max should be equal (flat extrapolation)
  for (i in seq_len(ncol(preds$Probs))) {
    expect_equal(preds$Probs[last_obs_idx, i], preds$Probs[beyond_idx, i])
  }
})

test_that("survprobMatInterpolator handles single time point", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Single time point
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = 10)

  expect_true(is.matrix(preds$Probs))
  expect_equal(nrow(preds$Probs), 1)
  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_equal(preds$Times, 10)
})

test_that("survprobMatInterpolator handles many time points", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Many time points (stress test)
  many_times <- seq(0, model$time_range[2], length.out = 500)
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = many_times)

  expect_equal(nrow(preds$Probs), 500)
  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_equal(length(preds$Times), 500)

  # All values should be valid probabilities
  expect_true(all(preds$Probs >= 0 & preds$Probs <= 1))
})

test_that("survprobMatInterpolator handles irregular time grid", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Irregular grid with clusters and gaps
  irregular_times <- c(0, 0.1, 0.2, 0.3, 5, 10, 10.1, 10.2, 20, 30, 40)
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = irregular_times)

  expect_equal(nrow(preds$Probs), length(irregular_times))
  expect_equal(preds$Times, irregular_times)

  # Should still be monotonic
  for (i in seq_len(ncol(preds$Probs))) {
    surv_curve <- preds$Probs[, i]
    expect_true(all(diff(surv_curve) <= 1e-10))
  }
})

test_that("Default time grid is generated when newtimes=NULL", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Don't specify times
  preds <- Predict_SurvModel_Cox(model, test_data)

  # Should generate default grid
  expect_equal(length(preds$Times), 50)  # Default is 50 points
  expect_equal(min(preds$Times), 0)
  expect_true(max(preds$Times) <= model$time_range[2])

  # Should be evenly spaced
  time_diffs <- diff(preds$Times)
  expect_true(all(abs(time_diffs - mean(time_diffs)) < 1e-10))
})

test_that("Interpolation works across different models", {
  # Fit different types of models
  model1 <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event",
                         varsel = "none")
  model2 <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event",
                         varsel = "backward", penalty = "AIC")

  skip_if_not_installed("glmnet")
  model3 <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event",
                         varsel = "penalized", nfolds = 3)

  # Use same time grid for all
  common_times <- c(5, 10, 15, 20)

  preds1 <- Predict_SurvModel_Cox(model1, test_data, newtimes = common_times)
  preds3 <- Predict_SurvModel_Cox(model3, test_data, newtimes = common_times)

  # Compare models 1 and 3
  expect_equal(dim(preds1$Probs), dim(preds3$Probs))
  expect_equal(preds1$Times, common_times)
  expect_equal(preds3$Times, common_times)

  # Model 2 only if it has variables
  if (length(coef(model2$cph_model)) > 0) {
    preds2 <- Predict_SurvModel_Cox(model2, test_data, newtimes = common_times)
    expect_equal(dim(preds1$Probs), dim(preds2$Probs))
    expect_equal(preds2$Times, common_times)
  }
})

# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("Interpolation enables ensemble averaging", {
  skip_if_not_installed("glmnet")

  # Fit multiple models (use penalized to ensure both work)
  model1 <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event",
                         varsel = "none")
  model2 <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event",
                         varsel = "penalized", alpha = 1, nfolds = 3)

  # Get predictions on common grid
  common_times <- seq(0, 30, by = 5)
  preds1 <- Predict_SurvModel_Cox(model1, test_data, newtimes = common_times)
  preds2 <- Predict_SurvModel_Cox(model2, test_data, newtimes = common_times)

  # Average predictions (simple average on probability scale)
  avg_probs <- (preds1$Probs + preds2$Probs) / 2

  # Averaged predictions should also be valid
  expect_true(all(avg_probs >= 0 & avg_probs <= 1))

  # Should maintain monotonicity
  for (i in seq_len(ncol(avg_probs))) {
    expect_true(all(diff(avg_probs[, i]) <= 1e-10))
  }
})

test_that("Can extract risk scores from interpolated predictions", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Get predictions on fine grid
  times <- seq(0, 30, by = 0.5)
  preds <- Predict_SurvModel_Cox(model, test_data, newtimes = times)

  # Calculate AUC (area under 1-S(t) curve) as risk score
  risk_scores <- numeric(ncol(preds$Probs))
  for (i in seq_len(ncol(preds$Probs))) {
    surv_curve <- preds$Probs[, i]
    event_prob <- 1 - surv_curve
    # Trapezoidal integration
    risk_scores[i] <- sum(diff(times) *
                         (event_prob[-1] + event_prob[-length(event_prob)]) / 2)
  }

  # Risk scores should be positive and reasonable
  expect_true(all(risk_scores > 0))
  expect_true(all(risk_scores < max(times)))  # Can't exceed time horizon

  # Higher risk should correspond to lower survival
  avg_surv_at_15 <- preds$Probs[times == 15, ]
  expect_true(cor(risk_scores, -avg_surv_at_15) > 0.5)
})

test_that("Interpolation preserves probability at observed times", {
  model <- SurvModel_Cox(train_data, c("x1", "x2"), "time", "event")

  # Get predictions at default times (includes actual event times)
  preds_default <- Predict_SurvModel_Cox(model, test_data[1:3, ])

  # Now request interpolation that includes some of those exact times
  times_subset <- preds_default$Times[c(1, 5, 10, 20, 30)]
  preds_interp <- Predict_SurvModel_Cox(model, test_data[1:3, ],
                                        newtimes = times_subset)

  # Extract probabilities at matching times
  for (t in times_subset) {
    idx_default <- which(preds_default$Times == t)
    idx_interp <- which(preds_interp$Times == t)

    if (length(idx_default) > 0 && length(idx_interp) > 0) {
      expect_equal(preds_default$Probs[idx_default, ],
                   preds_interp$Probs[idx_interp, ],
                   tolerance = 1e-6)
    }
  }
})
