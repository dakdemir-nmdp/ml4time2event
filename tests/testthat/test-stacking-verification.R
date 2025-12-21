# Verification Test for Super Learner Stacking Optimization
# Checks if weights are optimized (non-uniform) when models perform differently.

library(testthat)
library(ml4time2event)
library(survival)

test_that("Super Learner Stacking optimizes weights", {
    skip_if_not_installed("survival")

    set.seed(123)
    n <- 200
    # Simulate data where "model1" is good and "model2" is bad
    # Model 1 features: Strong signal
    # Model 2 features: Noise

    # Truth: hazard depends on x1
    x1 <- rnorm(n)
    h <- exp(2 * x1)
    time <- rexp(n, rate = h)
    event <- rep(1, n) # No censoring for simplicity of signal

    data <- data.frame(
        id = 1:n,
        time = time,
        event = event,
        x1 = x1, # Good feature
        noise = rnorm(n) # Bad feature
    )

    # We will manually create predictions to simulate models to avoid training complexity
    # or rely on simple models.
    # Let's use ml4t2e_fit with "cox" on meaningful features and maybe "cox" on noise?
    # But ml4t2e_fit takes one dataset.

    # Alternative: Use "random_forest" vs "cox" or just check if weights are non-uniform.
    # Let's fit separate models manually and then use optimizeSuperLearnerWeights directly.
    # This avoids integration complexity and tests the optimizer itself.

    time_grid <- sort(unique(time[time < quantile(time, 0.9)])) # 90% times
    if (length(time_grid) > 20) time_grid <- seq(min(time_grid), max(time_grid), length.out = 20)

    # Perfect Prediction (Model A)
    # Surv(t) = exp(-H(t)) = exp(-CumHaz(t))
    # True lambda = exp(2 * x1)
    # True CumHaz(t) = lambda * t
    # True Surv(t) = exp(-lambda * t)
    lambda <- exp(2 * x1)
    pred_A <- outer(lambda, time_grid, function(l, t) exp(-l * t))

    # Terrible Prediction (Model B) - Reverse/Random
    # Say Model B predicts constant hazard lambda=1
    pred_B <- outer(rep(1, n), time_grid, function(l, t) exp(-l * t))

    predictions_list <- list(
        ModelA = pred_A,
        ModelB = pred_B
    )

    # Actual matrix (0/1 for alive/dead at t)
    # For optimization, we need 'actual' status.
    # At time t:
    # If T > t: Alive (1)
    # If T <= t: Dead (0)
    actual_mat <- outer(time, time_grid, function(T, t) as.numeric(T > t))

    # Run Optimization
    # Need to export optimizeSuperLearnerWeights or use :::
    weights_mse <- ml4time2event:::optimizeSuperLearnerWeights(
        predictions_list,
        actual_mat,
        loss_type = "mse"
    )

    message("MSE Weights: ", paste(names(weights_mse), round(weights_mse, 3), sep = "=", collapse = ", "))

    # Expect ModelA to have much higher weight
    expect_gt(weights_mse["ModelA"], 0.8)
    expect_lt(weights_mse["ModelB"], 0.2)

    # Check LogLik Optimization
    weights_ll <- ml4time2event:::optimizeSuperLearnerWeights(
        predictions_list,
        actual_mat,
        loss_type = "loglik"
    )

    message("LogLik Weights: ", paste(names(weights_ll), round(weights_ll, 3), sep = "=", collapse = ", "))

    expect_gt(weights_ll["ModelA"], 0.8)
    expect_lt(weights_ll["ModelB"], 0.2)
})

test_that("Stacking Integration in Pipeline", {
    # Test end-to-end integration
    skip_if_not_installed("survival")

    set.seed(300)
    n <- 100
    data <- data.frame(
        id = 1:n,
        time = rexp(n),
        event = rbinom(n, 1, 0.8),
        x = rnorm(n)
    )

    task <- ml4t2e_task_surv(data, "time", "event", "x", id = "id")

    # Enable stacking
    fit <- ml4t2e_fit(
        task = task,
        models = c("cox"), # Only 1 model, so weights should be 1.0 (trivial case)
        ensemble = "stack",
        conformal_calibration = 0.2, # Needed for 3-way split to have calibration data for stacking
        keep_data = TRUE
    )

    expect_equal(fit$ensemble$strategy, "stack")
    expect_true(!is.null(fit$ensemble$weights))
    expect_equal(sum(fit$ensemble$weights), 1)
})
