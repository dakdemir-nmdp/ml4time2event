# ==============================================================================
# Test Suite for Random Forest Survival Model
# Using unified ml4t2e_fit API
# ==============================================================================

library(testthat)
library(survival)

context("Testing Random Forest Survival Model via ml4t2e_fit")

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(789)
n_obs_surv <- 50
surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)), # 0=censored, 1=event
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)),
  x3 = rnorm(n_obs_surv, mean = 2),
  stringsAsFactors = FALSE
)
time_var <- "time"
event_var <- "status"
train_indices_surv <- 1:40
test_indices_surv <- 41:50
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]

train_task <- ml4t2e_task_surv(train_data_surv, time = time_var, event = event_var)
test_task <- ml4t2e_task_surv(test_data_surv, time = time_var, event = event_var)

# ==============================================================================
# Tests via ml4t2e_fit
# ==============================================================================

test_that("ml4t2e_fit fits Random Forest model", {
  skip_if_not_installed("randomForestSRC")

  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "random_forest",
    ensemble = "none",
    controls = list(random_forest = list(ntree = 10)) # Faster tests
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("random_forest" %in% fit$model_names)

  # Check internal model structure
  rf_obj <- fit$models$random_forest
  expect_true(inherits(rf_obj, "RandomForestSurvivalModel"))
  expect_s3_class(rf_obj$model, "rfsrc")
  expect_equal(rf_obj$model$ntree, 10)
})

test_that("Random Forest predictions format", {
  skip_if_not_installed("randomForestSRC")

  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "random_forest",
    ensemble = "none",
    controls = list(random_forest = list(ntree = 10))
  )

  times <- c(1, 5, 10)
  preds <- predict(fit, newdata = test_data_surv, times = times, type = "survival")

  expect_s3_class(preds, "t2e_pred")
  expect_true(all(c("id", "time", "surv", "model") %in% names(preds)))
  expect_equal(nrow(preds), nrow(test_data_surv) * length(times))

  # Check values
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})

test_that("Random Forest handles factor variables without manual intervention", {
  skip_if_not_installed("randomForestSRC")

  # x2 is a factor
  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "random_forest",
    ensemble = "none",
    controls = list(random_forest = list(ntree = 10))
  )

  expect_no_error(predict(fit, newdata = test_data_surv, times = 5))
})
