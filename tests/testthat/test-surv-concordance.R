# ==============================================================================
# Test Suite for Survival Concordance Index
# ==============================================================================

library(testthat)
library(ml4time2event)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated survival data with known concordance
set.seed(42)
n_obs <- 200

# Simulate data where higher x1 leads to higher risk (shorter survival)
beta <- 3.0
x1 <- rnorm(n_obs)
lp <- beta * x1
time <- rexp(n_obs, rate = exp(lp - mean(lp) + log(0.1)))
  censor_time <- rexp(n_obs, rate = 0.01)  # Low censoring
event <- as.integer(time <= censor_time)
observed_time <- pmin(time, censor_time)

surv_data <- data.frame(
  time = observed_time,
  event = event,
  x1 = x1
)

# Split into train/test
train_idx <- sample(seq_len(n_obs), size = floor(0.7 * n_obs))
train_data <- surv_data[train_idx, ]
test_data <- surv_data[-train_idx, ]

# ==============================================================================
# Tests
# ==============================================================================

test_that("Survival concordance index is > 0.5 for models with predictive power", {
  # Create survival task
  surv_task <- ml4t2e_task_surv(
    data = train_data,
    time = "time",
    event = "event",
    features = "x1",
    time_units = "days"
  )

  # Fit a Cox model (should have good concordance since x1 is predictive)
  surv_fit <- ml4t2e_fit(
    task = surv_task,
    models = "cox",
    controls = list(times = seq(0, max(train_data$time), length.out = 10))
  )

  # Evaluate on training data
  train_metrics <- ml4t2e_evaluate(
    surv_fit,
    metrics = "c_index"
  )

  # Concordance should be > 0.5 for a model with predictive power
  c_index <- train_metrics$value[train_metrics$metric == "c_index" & train_metrics$model == "cox"]
  print(paste("C-index:", c_index))
  expect_gt(c_index, 0.5)
  expect_lte(c_index, 1.0)

  # Test on test data
  surv_task_test <- ml4t2e_task_surv(
    data = test_data,
    time = "time",
    event = "event",
    features = "x1",
    time_units = "days"
  )

  test_preds <- predict(
    surv_fit,
    newdata = test_data,
    times = seq(0, max(test_data$time), length.out = 10),
    type = "survival"
  )

  test_metrics <- ml4t2e_evaluate(
    test_preds,
    task = surv_task_test,
    metrics = "c_index"
  )

  c_index_test <- test_metrics$value[test_metrics$metric == "c_index" & test_metrics$model == "cox"]
  expect_gt(c_index_test, 0.5)
  expect_lte(c_index_test, 1.0)
})

test_that("Survival concordance index handles multiple models", {
  surv_task <- ml4t2e_task_surv(
    data = train_data,
    time = "time",
    event = "event",
    features = "x1",
    time_units = "days"
  )

  # Fit multiple models
  surv_fit <- ml4t2e_fit(
    task = surv_task,
    models = c("cox", "random_forest"),
    controls = list(times = seq(0, max(train_data$time), length.out = 10))
  )

  # Evaluate
  metrics <- ml4t2e_evaluate(
    surv_fit,
    metrics = "c_index",
    include = c("cox", "random_forest")
  )

  # All concordance indices should be > 0.5
  c_indices <- metrics$value[metrics$metric == "c_index"]
  expect_true(all(c_indices > 0.5))
  expect_true(all(c_indices <= 1.0))
})