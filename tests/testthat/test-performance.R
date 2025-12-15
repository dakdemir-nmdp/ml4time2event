# Performance Tests for ml4time2event
#
# Tests that models fit and predict within reasonable time limits.
# Helps catch performance regressions early.

library(testthat)
library(ml4time2event)

test_that("ShallowNN survival fits 1000 obs in <10 seconds", {
    skip_on_cran()
    skip_if_not_installed("survival")

    # Generate test data
    set.seed(123)
    n <- 1000
    p <- 10

    test_data <- data.frame(
        id = 1:n,
        time = rexp(n, rate = 0.01),
        event = rbinom(n, 1, 0.7),
        matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    # Time the fit
    elapsed <- system.time({
        fit <- ml4t2e_fit(
            task = task,
            models = "shallownn",
            keep_data = TRUE,
            controls = list(shallownn = list(maxit = 100)) # Limit iterations for speed
        )
    })["elapsed"]

    expect_lt(elapsed, 10,
        label = sprintf("ShallowNN took %.2f seconds (limit: 10s)", elapsed)
    )
})

test_that("Cox model fits 1000 obs in <5 seconds", {
    skip_on_cran()
    skip_if_not_installed("survival")

    set.seed(124)
    n <- 1000
    p <- 10

    test_data <- data.frame(
        id = 1:n,
        time = rexp(n, rate = 0.01),
        event = rbinom(n, 1, 0.7),
        matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    elapsed <- system.time({
        fit <- ml4t2e_fit(
            task = task,
            models = "cox",
            keep_data = TRUE
        )
    })["elapsed"]

    expect_lt(elapsed, 5,
        label = sprintf("Cox took %.2f seconds (limit: 5s)", elapsed)
    )
})

test_that("Random Forest fits 500 obs in <15 seconds", {
    skip_on_cran()
    skip_if_not_installed("randomForestSRC")

    set.seed(125)
    n <- 500 # Smaller for RF
    p <- 10

    test_data <- data.frame(
        id = 1:n,
        time = rexp(n, rate = 0.01),
        event = rbinom(n, 1, 0.7),
        matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    elapsed <- system.time({
        fit <- ml4t2e_fit(
            task = task,
            models = "random_forest",
            keep_data = TRUE,
            controls = list(random_forest = list(ntree = 100)) # Limit trees
        )
    })["elapsed"]

    expect_lt(elapsed, 15,
        label = sprintf("Random Forest took %.2f seconds (limit: 15s)", elapsed)
    )
})

test_that("Predictions are fast for all models", {
    skip_on_cran()
    skip_if_not_installed("survival")

    set.seed(126)
    n_train <- 200
    n_test <- 100
    p <- 5

    # Training data
    test_data <- data.frame(
        id = 1:n_train,
        time = rexp(n_train, rate = 0.01),
        event = rbinom(n_train, 1, 0.7),
        matrix(rnorm(n_train * p), nrow = n_train, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    # Test data
    newdata <- data.frame(
        id = (n_train + 1):(n_train + n_test),
        matrix(rnorm(n_test * p), nrow = n_test, ncol = p)
    )
    names(newdata)[2:(1 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    fit <- ml4t2e_fit(
        task = task,
        models = c("cox", "shallownn"),
        keep_data = TRUE,
        controls = list(shallownn = list(maxit = 50))
    )

    # Prediction should be fast (<2s for 100 obs x 2 models)
    elapsed <- system.time({
        preds <- predict(fit, newdata = newdata, times = seq(0, 500, length.out = 50))
    })["elapsed"]

    expect_lt(elapsed, 2,
        label = sprintf("Prediction took %.2f seconds (limit: 2s)", elapsed)
    )
})

test_that("Ensemble with stacking overhead is reasonable", {
    skip_on_cran()
    skip_if_not_installed("survival")

    set.seed(127)
    n <- 300
    p <- 5

    test_data <- data.frame(
        id = 1:n,
        time = rexp(n, rate = 0.01),
        event = rbinom(n, 1, 0.7),
        matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    # Fit without ensemble
    time_no_ensemble <- system.time({
        fit1 <- ml4t2e_fit(
            task = task,
            models = c("cox", "shallownn"),
            ensemble = FALSE,
            keep_data = TRUE,
            controls = list(shallownn = list(maxit = 50))
        )
    })["elapsed"]

    # Fit with stacking
    time_with_stack <- system.time({
        fit2 <- ml4t2e_fit(
            task = task,
            models = c("cox", "shallownn"),
            ensemble = "stack",
            keep_data = TRUE,
            controls = list(shallownn = list(maxit = 50))
        )
    })["elapsed"]

    # Stacking overhead should be <100% of base time (relatively relaxed for CI)
    # BUT: Skip if base time is too short for reliable percentage measurement
    overhead_ratio <- (time_with_stack - time_no_ensemble) / time_no_ensemble

    # Skip test if base time is <0.1s (overhead percentages become unreliable)
    skip_if(
        time_no_ensemble < 0.1,
        sprintf("Base time %.3fs is too short for reliable overhead measurement", time_no_ensemble)
    )

    expect_lt(overhead_ratio, 1.0,
        label = sprintf(
            "Stacking overhead %.1f%% exceeds 100%% limit (%.2fs base, %.2fs with stacking)",
            overhead_ratio * 100, time_no_ensemble, time_with_stack
        )
    )
})

test_that("Conformal calibration overhead is minimal", {
    skip_on_cran()
    skip_if_not_installed("survival")

    set.seed(128)
    n <- 300
    p <- 5

    test_data <- data.frame(
        id = 1:n,
        time = rexp(n, rate = 0.01),
        event = rbinom(n, 1, 0.7),
        matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    names(test_data)[4:(3 + p)] <- paste0("x", 1:p)

    task <- ml4t2e_task_surv(
        data = test_data,
        time = "time",
        event = "event",
        features = paste0("x", 1:p),
        id = "id"
    )

    # Fit without conformal
    time_no_conformal <- system.time({
        fit1 <- ml4t2e_fit(
            task = task,
            models = "cox",
            keep_data = TRUE
        )
    })["elapsed"]

    # Fit with conformal
    time_with_conformal <- system.time({
        fit2 <- ml4t2e_fit(
            task = task,
            models = "cox",
            conformal_calibration = 0.2,
            keep_data = TRUE
        )
    })["elapsed"]

    # Conformal overhead should be <50% of base time (relaxed for CI variability)
    # BUT: Skip if base time is too short for reliable percentage measurement
    overhead_ratio <- (time_with_conformal - time_no_conformal) / time_no_conformal

    # Skip test if base time is <0.05s (overhead percentages become unreliable)
    skip_if(
        time_no_conformal < 0.05,
        sprintf("Base time %.3fs is too short for reliable overhead measurement", time_no_conformal)
    )

    expect_lt(overhead_ratio, 0.5,
        label = sprintf(
            "Conformal overhead %.1f%% exceeds 50%% limit (%.2fs base, %.2fs with conformal)",
            overhead_ratio * 100, time_no_conformal, time_with_conformal
        )
    )
})

test_that("Baseline hazard computation is vectorized (performance check)", {
    # This tests that the vectorized baseline hazard is actually fast
    skip_on_cran()

    set.seed(129)
    n <- 2000 # Large dataset

    # Create data with many unique event times
    times <- rexp(n, rate = 0.01)
    events <- rbinom(n, 1, 0.7)
    risks <- exp(rnorm(n))

    # Sort descending (as the function expects)
    ord <- order(-times, -events)
    time_sorted <- times[ord]
    event_sorted <- events[ord]
    risks_sorted <- risks[ord]

    # Time the vectorized version
    elapsed <- system.time({
        result <- .nn_compute_baseline_hazard_vectorized(
            risks = risks_sorted,
            time_sorted = time_sorted,
            event_sorted = event_sorted
        )
    })["elapsed"]

    # Should complete in <0.5 seconds even for 2000 obs
    expect_lt(elapsed, 0.5,
        label = sprintf("Baseline hazard took %.3f seconds for n=2000 (limit: 0.5s)", elapsed)
    )

    # Verify it returns something sensible
    expect_s3_class(result, "data.frame")
    expect_true(all(c("time", "hazard") %in% names(result)))
    expect_true(all(result$hazard >= 0))
    expect_true(all(diff(result$hazard) >= 0)) # Cumulative should be monotone increasing
})
