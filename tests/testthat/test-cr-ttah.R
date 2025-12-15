library(testthat)

test_that("ml4t2e_fit fits TTAH Competing Risk model", {
  set.seed(42)
  n_obs <- 100
  train_data <- data.frame(
    time = rexp(n_obs, rate = 0.5),
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

  # Fit
  fit <- ml4t2e_fit(keep_data = TRUE, task,
    models = "cr_ttah",
    controls = list(cr_ttah = list(n_time = 10))
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("cr_ttah" %in% fit$model_names)

  wrapper <- fit$models$cr_ttah
  expect_true(inherits(wrapper, "TtahCompetingRiskModel"))
  expect_true(!is.null(wrapper$model$multinom))
})

test_that("TTAH Competing Risk predictions format", {
  set.seed(99)
  n_obs <- 60
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
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "cr_ttah", controls = list(cr_ttah = list(n_time = 10)))

  test_data <- train_data[1:5, ]
  times <- c(0.5, 1.5)

  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("cif" %in% colnames(preds))
  expect_true("cause" %in% colnames(preds))
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))
  expect_true(length(unique(preds$cause)) >= 2)
})
