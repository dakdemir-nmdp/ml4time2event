library(testthat)
library(survival)
library(xgboost)

test_that("ml4t2e_fit fits XGBoost Survival model", {
  skip_if_not_installed("xgboost")

  # Setup data
  set.seed(789)
  n_obs <- 50
  surv_data <- data.frame(
    time = rexp(n_obs, rate = 0.05),
    status = sample(0:1, n_obs, replace = TRUE, prob = c(0.3, 0.7)),
    x1 = rnorm(n_obs),
    x2 = factor(sample(c("C", "D"), n_obs, replace = TRUE)),
    x3 = rnorm(n_obs, mean = 2),
    stringsAsFactors = FALSE
  )

  task <- ml4t2e_task_surv(
    surv_data,
    time = "time",
    event = "status"
  )

  # Fit via API
  fit <- ml4t2e_fit(keep_data = TRUE, 
    task,
    models = "xgboost",
    ensemble = "none"
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("xgboost" %in% fit$model_names)

  # Check Internal Structure
  xgb_obj <- fit$models$xgboost
  expect_true(inherits(xgb_obj, "XGBoostSurvivalModel"))

  # Check native model storage
  expect_s3_class(xgb_obj$model, "xgb.Booster")
  expect_true(is.data.frame(xgb_obj$baseline_hazard))
  expect_equal(colnames(xgb_obj$baseline_hazard), c("time", "hazard"))
})

test_that("XGBoost Survival predictions format", {
  skip_if_not_installed("xgboost")

  set.seed(999)
  n_obs <- 100
  surv_data <- data.frame(
    time = rexp(n_obs),
    status = rbinom(n_obs, 1, 0.5),
    x1 = rnorm(n_obs)
  )

  train_data <- surv_data[1:80, ]
  test_data <- surv_data[81:100, ]

  task <- ml4t2e_task_surv(train_data, time = "time", event = "status")
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "xgboost")

  # Predict
  times <- c(0.5, 1.0, 1.5)
  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("surv" %in% colnames(preds))
  expect_true("time" %in% colnames(preds))
  expect_true("id" %in% colnames(preds))

  # Values [0, 1]
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))

  # Monotonicity check (approximate due to interpolation)
  # Group by ID and check monotonic decreasing
  first_id_preds <- preds[preds$id == preds$id[1], ]
  first_id_preds <- first_id_preds[order(first_id_preds$time), ]
  # diff(surv) should be <= 0 (allow small epsilon)
  diffs <- diff(first_id_preds$surv)
  expect_true(all(diffs <= 1e-9))
})
