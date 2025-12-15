library(testthat)
library(mgcv)

test_that("ml4t2e_fit fits GAM Survival model", {
  skip_if_not_installed("mgcv")

  set.seed(42)
  n_obs <- 50
  surv_data <- data.frame(
    time = rexp(n_obs, rate = 0.1),
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
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "gam")

  expect_s3_class(fit, "t2e_fit")
  expect_true("gam" %in% fit$model_names)

  wrapper <- fit$models$gam
  expect_true(inherits(wrapper, "GamSurvivalModel"))
  expect_s3_class(wrapper$model, "gam")
})

test_that("GAM Survival predictions format", {
  skip_if_not_installed("mgcv")

  set.seed(42)
  n_obs <- 40
  surv_data <- data.frame(
    time = rexp(n_obs, rate = 0.1),
    status = sample(0:1, n_obs, replace = TRUE),
    x1 = rnorm(n_obs)
  )

  task <- ml4t2e_task_surv(surv_data, time = "time", event = "status")
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "gam")

  test_data <- surv_data[1:10, ]

  times <- c(1, 5)
  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("surv" %in% colnames(preds))
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})
