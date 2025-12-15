library(testthat)
library(partykit)
library(glmnet)

test_that("ml4t2e_fit fits RuleFit Survival model", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rpart")

  set.seed(42)
  n_obs <- 100
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

  # Fit with fewer trees for speed
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "rulefit", controls = list(rulefit = list(ntree = 5, nsample = 50)))

  expect_s3_class(fit, "t2e_fit")
  expect_true("rulefit" %in% fit$model_names)

  wrapper <- fit$models$rulefit
  expect_true(inherits(wrapper, "RulefitSurvivalModel"))
  expect_true(!is.null(wrapper$cv_fit))
  expect_type(wrapper$rules_list, "list")
})

test_that("RuleFit Survival predictions format", {
  skip_if_not_installed("partykit")

  set.seed(42)
  n_obs <- 100
  surv_data <- data.frame(
    time = rexp(n_obs, rate = 0.1),
    status = sample(0:1, n_obs, replace = TRUE),
    x1 = rnorm(n_obs),
    x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE))
  )

  task <- ml4t2e_task_surv(surv_data, time = "time", event = "status")
  # Ensure we generate enough valid trees/rules
  fit <- ml4t2e_fit(keep_data = TRUE, task, models = "rulefit", controls = list(rulefit = list(ntree = 5, nsample = 50)))

  test_data <- surv_data[1:5, ]

  times <- c(1, 5)
  preds <- predict(fit, newdata = test_data, times = times)

  expect_true("surv" %in% colnames(preds))
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})
