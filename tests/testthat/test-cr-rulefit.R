library(testthat)
library(partykit)
library(glmnet)
library(fastcmprsk)

test_that("ml4t2e_fit fits RuleFit Competing Risk model", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rpart")
  skip_if_not_installed("fastcmprsk")
  skip_if_not_installed("pseudo")

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

  # Use small params for speed
  fit <- ml4t2e_fit(keep_data = TRUE, task,
    models = "cr_rulefit",
    controls = list(cr_rulefit = list(ntree = 5, nsample = 50))
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("cr_rulefit" %in% fit$model_names)

  wrapper <- fit$models$cr_rulefit
  expect_true(inherits(wrapper, "RulefitCompetingRiskModel"))
  expect_type(wrapper$model, "list")
  # Expect sub-models per cause
  expect_true(length(wrapper$model) >= 2)
})

test_that("RuleFit Competing Risk predictions format", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("fastcmprsk")
  skip_if_not_installed("pseudo")

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
  fit <- ml4t2e_fit(keep_data = TRUE, task,
    models = "cr_rulefit",
    controls = list(cr_rulefit = list(ntree = 5, nsample = 50))
  )

  test_data <- train_data[1:5, ]
  times <- c(0.5, 1.5)

  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("cif" %in% colnames(preds))
  expect_true("cause" %in% colnames(preds))
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))

  expect_true(length(unique(preds$cause)) >= 2)
})
