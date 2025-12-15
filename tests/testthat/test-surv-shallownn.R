library(testthat)

test_that("ml4t2e_fit fits ShallowNN Survival model", {
    set.seed(42)
    n_obs <- 100
    surv_data <- data.frame(
        time = rexp(n_obs, rate = 0.5),
        status = sample(0:1, n_obs, replace = TRUE),
        x1 = rnorm(n_obs),
        x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE))
    )

    task <- ml4t2e_task_surv(
        surv_data,
        time = "time",
        event = "status"
    )

    # Fit
    fit <- ml4t2e_fit(keep_data = TRUE, task, models = "shallownn")

    expect_s3_class(fit, "t2e_fit")
    expect_true("shallownn" %in% fit$model_names)

    wrapper <- fit$models$shallownn
    expect_true(inherits(wrapper, "ShallowNNSurvivalModel"))
    expect_type(wrapper$model$weights, "list")
})

test_that("ShallowNN Survival predictions format", {
    set.seed(42)
    n_obs <- 100
    surv_data <- data.frame(
        time = rexp(n_obs, rate = 0.5),
        status = sample(0:1, n_obs, replace = TRUE),
        x1 = rnorm(n_obs),
        x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE))
    )

    task <- ml4t2e_task_surv(surv_data, time = "time", event = "status")
    fit <- ml4t2e_fit(keep_data = TRUE, task, models = "shallownn")

    test_data <- surv_data[1:5, ]

    times <- c(0.5, 1.5)
    preds <- predict(fit, newdata = test_data, times = times)

    expect_true("surv" %in% colnames(preds))
    expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})
