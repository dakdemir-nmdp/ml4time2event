# ==============================================================================
# Test Suite for Competing Risks Random Forest Model
# Using unified ml4t2e_fit API
# ==============================================================================

library(testthat)

context("Testing Random Forest Competing Risks Model via ml4t2e_fit")

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(123)
n_obs <- 70
cr_data <- data.frame(
  time = pmin(rexp(n_obs, rate = 0.1), rexp(n_obs, rate = 0.15)),
  status = sample(0:2, n_obs, replace = TRUE, prob = c(0.2, 0.4, 0.4)), # 0=censored, 1=event1, 2=event2
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 5),
  stringsAsFactors = FALSE
)

# Explicitly prep columns for CR task if needed, though raw columns should work
cr_data$status_bin <- ifelse(cr_data$status == 0, 0, 1)
cr_data$cause_val <- ifelse(cr_data$status == 0, NA, cr_data$status)

train_indices <- 1:50
test_indices <- 51:70
train_data <- cr_data[train_indices, ]
test_data <- cr_data[test_indices, ]

train_task <- ml4t2e_task_cr(train_data, time = "time", status = "status_bin", cause = "cause_val")
test_task <- ml4t2e_task_cr(test_data, time = "time", status = "status_bin", cause = "cause_val")

# ==============================================================================
# Tests via ml4t2e_fit
# ==============================================================================

test_that("ml4t2e_fit fits CR Random Forest model", {
  skip_if_not_installed("randomForestSRC")

  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "cr_random_forest",
    ensemble = "none",
    controls = list(cr_random_forest = list(ntree = 20))
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("cr_random_forest" %in% fit$model_names)
  expect_equal(attr(fit$task, "task_type"), "competing_risk")

  # Check internal model structure
  # Now using native R6 class
  rf_wrapper <- fit$models$cr_random_forest
  expect_true(inherits(rf_wrapper, "RandomForestCompetingRiskModel"))

  # The internal model is a LIST of RF models (one per cause)
  expect_true(is.list(rf_wrapper$model))
  expect_true(all(c("1", "2") %in% names(rf_wrapper$model)))

  # Check the model for cause "1"
  model_1 <- rf_wrapper$model[["1"]]
  # Direct rfsrc storage (no longer wrapped in list(rf_model=...))
  expect_s3_class(model_1, "rfsrc")
  expect_equal(model_1$ntree, 20)
})

test_that("CR Random Forest predictions format", {
  skip_if_not_installed("randomForestSRC")

  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "cr_random_forest",
    ensemble = "none",
    controls = list(cr_random_forest = list(ntree = 20))
  )

  times <- c(1, 5)
  preds <- predict(fit, newdata = test_data, times = times, type = "cif")

  expect_s3_class(preds, "t2e_pred")
  expect_true(all(c("id", "time", "cif", "model", "cause") %in% names(preds)))
  expect_equal(nrow(preds), nrow(test_data) * length(times) * 2) # 2 causes

  # Check values
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))
})

test_that("CR Random Forest handles factor variables", {
  skip_if_not_installed("randomForestSRC")

  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "cr_random_forest",
    ensemble = "none",
    controls = list(cr_random_forest = list(ntree = 20))
  )
  expect_no_error(predict(fit, newdata = test_data, times = 5))
})
