library(testthat)
library(survival)

test_that("ml4t2e_fit fits SurvReg Competing Risk model", {
  set.seed(42)
  n_obs <- 40
  train_data <- data.frame(
    time = rexp(n_obs, rate = 0.1),
    event = sample(0:2, n_obs, replace = TRUE),
    x1 = rnorm(n_obs),
    x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE))
  )
  while (length(unique(train_data$event)) < 3) {
    train_data$event <- sample(0:2, n_obs, replace = TRUE)
  }
  train_data$status <- ifelse(train_data$event == 0, 0, 1)

  task <- ml4t2e_task_cr(
    train_data,
    time = "time",
    status = "status",
    cause = "event"
  )

  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cr_survreg")

  expect_s3_class(fit, "t2e_fit")
  expect_true("cr_survreg" %in% fit$model_names)

  wrapper <- fit$models$cr_survreg
  expect_true(inherits(wrapper, "SurvRegCompetingRiskModel"))
  expect_type(wrapper$model, "list")
  expect_true(all(c("1", "2") %in% names(wrapper$model)))
})

test_that("SurvReg Competing Risk predictions format", {
  set.seed(99)
  n_obs <- 30
  train_data <- data.frame(
    time = rexp(n_obs),
    event = sample(0:2, n_obs, replace = TRUE),
    x1 = rnorm(n_obs)
  )
  while (length(unique(train_data$event)) < 3) {
    train_data$event <- sample(0:2, n_obs, replace = TRUE)
  }
  train_data$status <- ifelse(train_data$event == 0, 0, 1)

  task <- ml4t2e_task_cr(train_data, time = "time", status = "status", cause = "event")
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cr_survreg")

  test_data <- train_data[1:5, ]
  times <- c(1, 2)

  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("cif" %in% colnames(preds))
  expect_true("cause" %in% colnames(preds))
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))

  unique_causes <- unique(preds$cause)
  expect_true(length(unique_causes) >= 2)
})
