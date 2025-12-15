# ==============================================================================
# Test Suite for Competing Risks Concordance Index
# ==============================================================================

library(testthat)
library(ml4time2event)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated competing risks data with known concordance
set.seed(42)
n_obs <- 200

# Simulate data where higher x1 leads to higher risk for cause 1
beta <- 1.5
x1 <- rnorm(n_obs)
lp <- beta * x1

# Simulate competing risks times
time1 <- rexp(n_obs, rate = exp(lp - mean(lp) + log(0.05))) # Cause 1
time2 <- rexp(n_obs, rate = exp(-lp - mean(-lp) + log(0.05))) # Cause 2 (inverse relationship)
censor_time <- rexp(n_obs, rate = 0.03)

# Determine which event occurs first
min_time <- pmin(time1, time2, censor_time)
event <- ifelse(min_time == time1, 1L,
  ifelse(min_time == time2, 2L, 0L)
)

cr_data <- data.frame(
  time = min_time,
  event = event,
  x1 = x1
)

# Prepare for ml4t2e_task_cr
cr_data$status <- ifelse(cr_data$event == 0, 0, 1)
cr_data$cause <- ifelse(cr_data$event == 0, NA, cr_data$event)

# Split into train/test
train_idx <- sample(seq_len(n_obs), size = floor(0.7 * n_obs))
train_data <- cr_data[train_idx, ]
test_data <- cr_data[-train_idx, ]

# ==============================================================================
# Tests
# ==============================================================================

test_that("Competing risks concordance index is > 0.5 for models with predictive power", {
  # Create competing risks task
  cr_task <- ml4t2e_task_cr(
    data = train_data,
    time = "time",
    status = "status",
    cause = "cause",
    features = "x1",
    time_units = "days"
  )

  # Fit a Fine-Gray model (should have good concordance since x1 is predictive for cause 1)
  cr_fit <- ml4t2e_fit(keep_data = TRUE, 
    task = cr_task,
    models = "cr_fine_gray",
    controls = list(times = seq(0, max(train_data$time), length.out = 10))
  )

  # Evaluate on training data
  train_metrics <- ml4t2e_evaluate(
    cr_fit,
    metrics = "c_index"
  )

  # Concordance should be >= 0.5 for a model with predictive power
  c_indices <- train_metrics$value[train_metrics$metric == "c_index"]
  expect_true(all(c_indices >= 0.5))
  expect_true(all(c_indices <= 1.0))

  # Test on test data
  cr_task_test <- ml4t2e_task_cr(
    data = test_data,
    time = "time",
    status = "status",
    cause = "cause",
    features = "x1",
    time_units = "days"
  )

  test_preds <- predict(
    cr_fit,
    newdata = test_data,
    times = seq(0, max(test_data$time), length.out = 10),
    type = "cif"
  )

  test_metrics <- ml4t2e_evaluate(
    test_preds,
    task = cr_task_test,
    metrics = "c_index"
  )

  c_indices_test <- test_metrics$value[test_metrics$metric == "c_index"]
  expect_true(all(c_indices_test >= 0.5))
  expect_true(all(c_indices_test <= 1.0))
})

test_that("Competing risks concordance index handles multiple causes", {
  cr_task <- ml4t2e_task_cr(
    data = train_data,
    time = "time",
    status = "status",
    cause = "cause",
    features = "x1",
    time_units = "days"
  )

  # Fit cause-specific Cox model
  cr_fit <- ml4t2e_fit(keep_data = TRUE, 
    task = cr_task,
    models = "cox",
    controls = list(times = seq(0, max(train_data$time), length.out = 10))
  )

  # Evaluate
  metrics <- ml4t2e_evaluate(
    cr_fit,
    metrics = "c_index"
  )

  # Should have concordance for each cause
  expect_true("cause" %in% colnames(metrics))
  causes <- unique(metrics$cause)
  expect_length(causes, 2) # Two competing causes

  # All concordance indices should be >= 0.5
  c_indices <- metrics$value[metrics$metric == "c_index"]
  expect_true(all(c_indices >= 0.5))
  expect_true(all(c_indices <= 1.0))
})
