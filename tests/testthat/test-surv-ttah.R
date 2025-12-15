library(testthat)

test_that("ml4t2e_fit fits TTAH Survival model", {
  set.seed(42)
  n_obs <- 100
  surv_data <- data.frame(
    time = rexp(n_obs, rate = 0.5), # Continuous time will be discretized
    status = sample(0:1, n_obs, replace = TRUE),
    x1 = rnorm(n_obs),
    x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE))
  )

  task <- ml4t2e_task_surv(
    surv_data,
    time = "time",
    event = "status"
  )

  # Fit with small grid for speed
  fit <- ml4t2e_fit(keep_data = TRUE, task,
    models = "ttah",
    controls = list(ttah = list(n_time = 10))
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("ttah" %in% fit$model_names)

  wrapper <- fit$models$ttah
  expect_true(inherits(wrapper, "TtahSurvivalModel"))
  expect_true(!is.null(wrapper$model$grid))
})

test_that("TTAH Survival predictions format", {
  set.seed(123)
  n_obs <- 60
  surv_data <- data.frame(
    time = rexp(n_obs),
    status = sample(0:1, n_obs, replace = TRUE),
    x1 = rnorm(n_obs)
  )

  task <- ml4t2e_task_surv(surv_data, time = "time", event = "status")
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "ttah", controls = list(ttah = list(n_time = 10)))

  test_data <- surv_data[1:5, ]
  times <- c(0.5, 1.5)
  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("surv" %in% colnames(preds))
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})
