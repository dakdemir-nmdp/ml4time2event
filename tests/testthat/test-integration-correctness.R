# Integration Tests for Model Learnability
# Ensures that models can actually learn signal from data (C-index > 0.5)

library(testthat)

test_that("Integration: Models learn signal on synthetic data", {
    set.seed(42)
    n <- 300
    # Strong signal simulation
    # Hazard depends strongly on x
    x <- rnorm(n)
    beta <- 1.5
    # t ~ Exponential(exp(beta*x))
    # High x -> High hazard -> Short time
    time <- rexp(n, rate = exp(beta * x))

    # Cancellations (censoring)
    cens <- rexp(n, rate = 0.5)
    event <- as.integer(time <= cens)
    time <- pmin(time, cens)

    df <- data.frame(x = x, time = time, event = event)

    task <- ml4t2e_task_surv(df, "time", "event")

    # Fit Cox
    fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cox", ensemble = "none")
    res <- ml4t2e_evaluate(fit, metrics = "c_index")

    # C-index must be high (0.5 is random)
    expect_gt(res$value[1], 0.65)

    # Fit Random Forest (non-linear capable but should pick up linear too)
    # Limiting trees for speed
    fit_rf <- ml4t2e_fit(keep_data = TRUE, task,
        models = "random_forest", ensemble = "none",
        controls = list(random_forest = list(num.trees = 50))
    )
    res_rf <- ml4t2e_evaluate(fit_rf, metrics = "c_index")
    expect_gt(res_rf$value[1], 0.60)
})
