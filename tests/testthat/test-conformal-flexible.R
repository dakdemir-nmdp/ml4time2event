# Verification Test for Conformal Prediction Flexible Time Interpolation
# Checks if predictions at non-grid times work correctly via interpolation.

library(testthat)
library(ml4time2event)
library(survival)

test_that("Conformal prediction interpolates for non-grid times", {
    skip_if_not_installed("survival")

    set.seed(42)
    n <- 50
    df <- data.frame(
        time = rexp(n),
        status = sample(0:1, n, replace = TRUE),
        x = rnorm(n),
        id = 1:n
    )
    task <- ml4t2e_task_surv(df, "time", "status", "x", id = "id")

    # Train model with specific grid
    grid <- c(1, 2, 3)
    fit <- ml4t2e_fit(
        task = task,
        models = "cox",
        ensemble = "none",
        conformal_calibration = 0.2,
        keep_data = TRUE,
        controls = list(times = grid)
    )

    # Check that grid contains our requested points (might include 0)
    expect_true(all(grid %in% fit$time_grid))

    # Predict at NON-GRID times
    new_times <- c(1.5, 2.5) # Midpoints

    # This should NOT error or warn
    preds <- predict(
        fit,
        newdata = df[1:5, ],
        times = new_times,
        type = "survival",
        conformal_alpha = 0.1
    )

    expect_s3_class(preds, "t2e_pred")
    expect_equal(sort(unique(preds$time)), sort(new_times))
    expect_true("lower" %in% names(preds))
    expect_true("upper" %in% names(preds))

    # Check bounds logic
    expect_true(all(preds$lower <= preds$surv + 1e-5))
    expect_true(all(preds$surv <= preds$upper + 1e-5))

    # Predict OUTSIDE grid range (extrapolation behavior)
    # R's approx with rule=2 should clamp to nearest value
    out_times <- c(0.5, 3.5)
    preds_out <- predict(
        fit,
        newdata = df[1:5, ],
        times = out_times,
        type = "survival",
        conformal_alpha = 0.1
    )

    expect_equal(sort(unique(preds_out$time)), sort(out_times))
    # Should not be all NA
    expect_false(all(is.na(preds_out$lower)))
})

test_that("Conformal prediction interpolates for Competing Risks", {
    skip_if_not_installed("survival")

    set.seed(42)
    n <- 100
    # Correct Data Structure for ml4t2e_task_cr
    # status: 0 (censored), 1 (event)
    # cause: "risk1" or "risk2"

    time <- rexp(n)
    event_raw <- sample(0:2, n, replace = TRUE) # 0=cens, 1=risk1, 2=risk2

    status <- as.numeric(event_raw > 0)
    cause <- rep(NA_character_, n)
    cause[event_raw == 1] <- "risk1"
    cause[event_raw == 2] <- "risk2"

    df <- data.frame(
        time = time,
        status = status,
        cause = cause,
        x = rnorm(n),
        id = 1:n
    )

    # ml4t2e_task_cr expects 'cause' arg to be the column name
    task <- ml4t2e_task_cr(df, "time", "status", "cause", features = "x", id = "id")

    grid <- c(1, 2, 3)
    fit <- ml4t2e_fit(
        task = task,
        models = "cox",
        ensemble = "none",
        conformal_calibration = 0.2,
        keep_data = TRUE,
        controls = list(times = grid)
    )

    # Check if grid matches (approx)
    expect_true(all(grid %in% fit$time_grid))

    new_times <- c(1.5, 2.5)

    preds <- predict(
        fit,
        newdata = df[1:5, ],
        times = new_times,
        type = "cif",
        conformal_alpha = 0.1
    )

    expect_s3_class(preds, "t2e_pred")
    expect_equal(sort(unique(preds$time)), sort(new_times))
    expect_true("lower" %in% names(preds))
    expect_true("upper" %in% names(preds))

    # Check bounds
    expect_true(all(preds$lower <= preds$cif + 1e-5))
    expect_true(all(preds$cif <= preds$upper + 1e-5))
})
