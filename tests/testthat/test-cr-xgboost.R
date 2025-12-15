library(testthat)
library(xgboost)

test_that("ml4t2e_fit fits XGBoost Competing Risk model", {
  skip_if_not_installed("xgboost")

  set.seed(42)
  n_train <- 50
  train_data <- data.frame(
    time = rexp(n_train, rate = 0.1),
    event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
    x1 = rnorm(n_train),
    x2 = factor(sample(c("A", "B"), n_train, replace = TRUE)),
    stringsAsFactors = FALSE
  )

  # Ensure we have events of both types
  while (length(unique(train_data$event)) < 3) {
    train_data$event <- sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  }
  train_data$status <- ifelse(train_data$event == 0, 0, 1)

  task <- ml4t2e_task_cr(
    train_data,
    time = "time",
    status = "status",
    cause = "event"
  )

  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cr_xgboost")

  expect_s3_class(fit, "t2e_fit")
  expect_true("cr_xgboost" %in% fit$model_names)

  # Internal Structure
  xgb_wrapper <- fit$models$cr_xgboost
  expect_true(inherits(xgb_wrapper, "XGBoostCompetingRiskModel"))

  expect_true(is.list(xgb_wrapper$model))
  # Should have 1 and 2
  expect_true(all(c("1", "2") %in% names(xgb_wrapper$model)))

  m1 <- xgb_wrapper$model[["1"]]
  expect_s3_class(m1$model, "xgb.Booster")
  expect_true(is.data.frame(m1$baseline_hazard))
})

test_that("XGBoost Competing Risk predictions format", {
  skip_if_not_installed("xgboost")

  set.seed(42)
  n_train <- 60
  train_data <- data.frame(
    time = rexp(n_train, rate = 0.1),
    event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
    x1 = rnorm(n_train),
    x2 = rnorm(n_train)
  )
  while (length(unique(train_data$event)) < 3) {
    train_data$event <- sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  }
  train_data$status <- ifelse(train_data$event == 0, 0, 1)

  test_data <- train_data[1:10, ]

  task <- ml4t2e_task_cr(train_data, time = "time", status = "status", cause = "event")
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cr_xgboost")

  # Predict
  times <- c(1, 5, 10)
  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("cif" %in% colnames(preds))
  expect_true("cause" %in% colnames(preds))

  unique_causes <- unique(preds$cause)
  expect_true(length(unique_causes) >= 2) # Should predict for all causes usually

  # Values [0, 1]
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))

  # Check monotonicity per cause/id
  # Filter for one ID and one cause
  one_curve <- preds[preds$id == preds$id[1] & preds$cause == "1", ]
  one_curve <- one_curve[order(one_curve$time), ]
  if (nrow(one_curve) > 1) {
    expect_true(all(diff(one_curve$cif) >= -1e-8))
  }
})
