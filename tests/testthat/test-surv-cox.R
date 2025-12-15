# ==============================================================================
# Test Suite for Cox Proportional Hazards Survival Model
# Using unified ml4t2e_fit API
# ==============================================================================

library(testthat)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated survival data for testing
set.seed(42)
n_train <- 200
n_test <- 50

# Simulate data with strong covariate effects for penalized Cox
make_surv_data <- function(n) {
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n, mean = 1)
  x4 <- rnorm(n, mean = -1)
  cat1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE), levels = c("A", "B", "C"))
  cat2 <- factor(sample(c("Low", "High"), n, replace = TRUE), levels = c("Low", "High"))

  beta <- c(x1 = 1.2, x2 = -1.0, x3 = 0.8, x4 = -0.7)

  # Linear predictor
  lp <- beta["x1"] * x1 + beta["x2"] * x2 + beta["x3"] * x3 + beta["x4"] * x4 +
    ifelse(cat1 == "B", 0.8, ifelse(cat1 == "C", -0.8, 0)) +
    ifelse(cat2 == "High", 1.0, 0)

  # Exponential survival times
  time <- rexp(n, rate = exp(lp - mean(lp) + log(0.1)))
  # Censoring
  censor_time <- rexp(n, rate = 0.05)
  event <- as.integer(time <= censor_time)
  observed_time <- pmin(time, censor_time)

  data.frame(
    time = observed_time,
    event = event,
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    cat1 = cat1, cat2 = cat2
  )
}

train_df <- make_surv_data(n_train)
test_df <- make_surv_data(n_test)

# Create Tasks
train_task <- ml4t2e_task_surv(train_df, time = "time", event = "event")
test_task <- ml4t2e_task_surv(test_df, time = "time", event = "event")

# ==============================================================================
# Tests for Cox Model via ml4t2e_fit
# ==============================================================================

test_that("ml4t2e_fit fits basic Cox model", {
  fit <- ml4t2e_fit(
    keep_data = TRUE,
    task = train_task,
    models = "cox",
    ensemble = "none"
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("cox" %in% fit$model_names)

  # Check internal model structure
  cox_obj <- fit$models$cox
  expect_true(inherits(cox_obj, "CoxSurvivalModel"))
  expect_s3_class(cox_obj$model, "coxph")
})

test_that("Cox model handles factor variables", {
  # Train data already has factors cat1, cat2
  fit <- ml4t2e_fit(
    keep_data = TRUE,
    task = train_task,
    models = "cox",
    ensemble = "none"
  )
  expect_no_error(predict(fit, newdata = test_df, times = c(1, 5)))

  # Check coefficient names to verify factors were expanded
  coefs <- coef(fit$models$cox$model)
  expect_true(any(grepl("cat1", names(coefs))))
})

test_that("Cox model predictions match expected structure", {
  fit <- ml4t2e_fit(
    keep_data = TRUE,
    task = train_task,
    models = "cox",
    ensemble = "none"
  )

  times <- c(0.5, 1, 2)
  preds <- predict(fit, newdata = test_df, times = times, type = "survival")

  expect_s3_class(preds, "t2e_pred")
  expect_true(all(c("id", "time", "surv", "model") %in% names(preds)))
  expect_equal(nrow(preds), nrow(test_df) * length(times))

  # Check values
  expect_true(all(preds$surv >= 0 & preds$surv <= 1))
})

test_that("Cox model handles penalized regression (via glmnet) if specified", {
  skip_if_not_installed("glmnet")

  # Note: Currently 'cox' engine maps to survival::coxph.
  # 'glmnet' engine maps to glmnet.
  # If the user wants penalized Cox, they should use model='glmnet'.

  fit_glmnet <- ml4t2e_fit(
    keep_data = TRUE,
    task = train_task,
    models = "glmnet", # Explicitly use glmnet engine
    ensemble = "none",
    controls = list(glmnet = list(alpha = 0.5))
  )

  expect_true("glmnet" %in% fit_glmnet$model_names)
  expect_true(inherits(fit_glmnet$models$glmnet, "GlmnetSurvivalModel"))
})

test_that("Cox model fails gracefully with invalid data", {
  invalid_df <- train_df
  invalid_df$event <- rep(0, nrow(invalid_df)) # No events

  # Task creation itself validates data and should error if no events
  expect_error(
    ml4t2e_task_surv(invalid_df, "time", "event"),
    "All observations are censored.*Cannot fit time-to-event model"
  )

  # Test missing feature columns behavior
  invalid_test <- test_df
  invalid_test$x1 <- NULL

  fit <- ml4t2e_fit(keep_data = TRUE, train_task, models = "cox")
  expect_error(predict(fit, newdata = invalid_test), "`newdata` is missing columns")
})
