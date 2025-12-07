test_that("conformal prediction works for survival", {
  skip_on_cran()
  library(survival)
  
  # Create simple task
  lung <- survival::lung
  lung <- na.omit(lung)
  lung$status <- lung$status - 1 # 0=censored, 1=dead
  
  task <- ml4t2e_task_surv(lung, time = "time", event = "status")
  
  # Fit with conformal calibration
  set.seed(123)
  fit <- ml4t2e_fit(
    task,
    models = "cox",
    conformal_calibration = 0.5 # Use 50% for calibration to ensure enough data
  )
  
  expect_true(!is.null(fit$conformal))
  expect_true("cox" %in% names(fit$conformal$scores))
  
  # Predict without bands
  times_to_pred <- fit$time_grid[c(10, 20)]
  preds <- predict(fit, newdata = lung[1:5, ], times = times_to_pred)
  expect_false("lower" %in% colnames(preds))
  
  # Predict with bands
  preds_bands <- predict(fit, newdata = lung[1:5, ], times = times_to_pred, conformal_alpha = 0.1)
  expect_true("lower" %in% colnames(preds_bands))
  expect_true("upper" %in% colnames(preds_bands))
  
  # Check bands logic
  expect_true(all(preds_bands$lower <= preds_bands$surv))
  expect_true(all(preds_bands$upper >= preds_bands$surv))
  expect_true(all(preds_bands$lower >= 0))
  expect_true(all(preds_bands$upper <= 1))
})

test_that("conformal prediction works for competing risks", {
  skip_on_cran()
  
  # Create simple CR task
  mgus2 <- survival::mgus2
  mgus2 <- na.omit(mgus2)
  # etime is time, pstat is event (0=cens, 1=PCM, 2=Death)
  
  task <- ml4t2e_task_cr(mgus2, time = "ptime", status = "pstat", cause = "pstat")
  
  # Fit with conformal calibration
  set.seed(123)
  fit <- ml4t2e_fit(
    task,
    models = "cox", # cr_cox
    conformal_calibration = 0.5
  )
  
  expect_true(!is.null(fit$conformal))
  expect_true("cox" %in% names(fit$conformal$scores))
  
  # Predict with bands
  times_to_pred_cr <- fit$time_grid[c(10, 20)]
  preds_bands <- predict(fit, newdata = mgus2[1:5, ], times = times_to_pred_cr, conformal_alpha = 0.1)
  
  expect_true("lower" %in% colnames(preds_bands))
  expect_true("upper" %in% colnames(preds_bands))
  
  # Check bands logic
  expect_true(all(preds_bands$lower <= preds_bands$cif))
  expect_true(all(preds_bands$upper >= preds_bands$cif))
  expect_true(all(preds_bands$lower >= 0))
  expect_true(all(preds_bands$upper <= 1))
})
