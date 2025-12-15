# Unit tests for Conformal Calibration logic
library(testthat)
library(survival)

test_that("Weighted Quantile works correctly", {
    # Equal weights
    x <- c(10, 20, 30)
    w <- c(1, 1, 1)
    # Total weight 3.
    # Alpha 0.1 -> 90% quantile. Target cumulative weight = 0.9 * 3 = 2.7.
    # Cum weights: 1, 2, 3.
    # First index >= 2.7 is 3. Value 30.
    expect_equal(ml4t2e_weighted_quantile(x, w, alpha = 0.1), 30)

    # Skewed weights
    x <- c(10, 20, 30)
    w <- c(1, 10, 1)
    # Total 12. Target 0.9 * 12 = 10.8.
    # Cum weights: 1, 11, 12.
    # First index >= 10.8 is 2. Value 20.
    expect_equal(ml4t2e_weighted_quantile(x, w, alpha = 0.1), 20)

    # Small coverage (alpha 0.9 -> 10% quantile)
    # Target 0.1 * 12 = 1.2.
    # First index >= 1.2 is 2. Value 20.
    expect_equal(ml4t2e_weighted_quantile(x, w, alpha = 0.9), 20)

    # Very small coverage
    # Target 0.05 * 12 = 0.6.
    # First index >= 0.6 is 1. Value 10.
    expect_equal(ml4t2e_weighted_quantile(x, w, alpha = 0.95), 10)
})

test_that("Internal Score Computation Logic (Survival)", {
    # Setup synthetic data
    # S1: T=10, E=1 (Event)
    # S2: T=5, E=0 (Censored)
    # S3: T=20, E=1 (Event)
    data <- data.frame(
        time = c(10, 5, 20),
        event = c(1, 0, 1),
        id = 1:3
    )

    # Mock task
    task <- list(time_col = "time", event_col = "event")

    # Mock Censoring Model: G(t) = 1 for all t (No censoring effect)
    # We construct a survfit that returns 1
    # To force ml4t2e_predict_censoring to returned 1, we can just mock the function or pass a trivial constant model
    # Let's mock ml4t2e_predict_censoring via with_mock (if available) or just construct a trivial survfit
    # Trivial survfit: 1 observation censored at very late time.
    t_cens_data <- data.frame(t = 1000, e = 0)
    cens_model <- survival::survfit(Surv(t, e == 0) ~ 1, data = t_cens_data)

    time_grid <- c(8, 15)

    # Predictions (Arbitrary)
    # S1: p=0.9 at t=8, p=0.8 at t=15
    # S2: p=0.9 at t=8, p=0.8 at t=15
    # S3: p=0.9 at t=8, p=0.8 at t=15
    pred_matrix <- matrix(
        c(
            0.9, 0.9, 0.9, # t=8
            0.8, 0.8, 0.8
        ), # t=15
        nrow = 3, ncol = 2
    )

    res <- ml4time2event:::.compute_scores_core(pred_matrix, data, time_grid, cens_model, task, event_of_interest = NULL)

    scores <- res$scores
    weights <- res$weights

    # Check Weights
    # t=8
    # S1: T=10 > 8. Alive. W=1.
    # S2: T=5 < 8. Censored. W=0.
    # S3: T=20 > 8. Alive. W=1.
    expect_equal(weights[, 1], c(1, 0, 1), tolerance = 1e-5)

    # t=15
    # S1: T=10 <= 15. Event=1. Dead. W=1.
    # S2: T=5. Censored. W=0.
    # S3: T=20 > 15. Alive. W=1.
    expect_equal(weights[, 2], c(1, 0, 1), tolerance = 1e-5)

    # Check Scores (Abs Error)
    # t=8
    # S1: Alive. Target=1. Err = |1 - 0.9| = 0.1.
    # S2: Censored. NA.
    # S3: Alive. Target=1. Err = |1 - 0.9| = 0.1.
    expect_equal(scores[1, 1], 0.1)
    expect_true(is.na(scores[2, 1]))
    expect_equal(scores[3, 1], 0.1)

    # t=15
    # S1: Dead. Target=0. Err = |0 - 0.8| = 0.8.
    # S2: Censored. NA.
    # S3: Alive. Target=1. Err = |1 - 0.8| = 0.2.
    expect_equal(scores[1, 2], 0.8)
    expect_true(is.na(scores[2, 2]))
    expect_equal(scores[3, 2], 0.2)
})

test_that("Integration: ml4t2e_calibrate attaches scores", {
    # Create simple task
    set.seed(42)
    n <- 50
    df <- data.frame(
        time = rexp(n),
        status = sample(0:1, n, replace = TRUE),
        x = rnorm(n),
        id = 1:n
    )
    task <- ml4t2e_task_surv(df, "time", "status", "x", id = "id")

    # Fit model
    fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cox", ensemble = "none")

    # Calibrate
    # Use subset as calibration data
    cal_df <- df[1:20, ]

    fit_cal <- ml4t2e_calibrate(fit, cal_df)

    expect_true("conformal" %in% names(fit_cal))
    expect_true("scores" %in% names(fit_cal$conformal))
    expect_true("cox" %in% names(fit_cal$conformal$scores))

    # Check predict bounds
    preds <- predict(fit_cal, newdata = df[1:5, ], type = "survival", conformal_alpha = 0.1)
    expect_true("lower" %in% colnames(preds))
    expect_true("upper" %in% colnames(preds))

    # Bounds should bracket prediction (unless clamped)
    # lower <= surv <= upper
    expect_true(all(preds$lower <= preds$surv + 1e-5))
    expect_true(all(preds$surv <= preds$upper + 1e-5))
})
