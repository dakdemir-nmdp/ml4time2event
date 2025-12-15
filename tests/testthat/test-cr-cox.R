# ==============================================================================
# Test Suite for Competing Risks Cox Model
# Using unified ml4t2e_fit API
# ==============================================================================

library(testthat)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(42)
# Simulate competing risks data
make_cr_data <- function(n) {
  # Create weaker signal to avoid perfect separation
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Very small coefficients
  lp1 <- 0.1 * x1
  lp2 <- -0.1 * x2

  # Simplified simulation
  # We just random assign event type with some dependency on x
  # This avoids perfect separation or lack of convergence
  risk1 <- exp(lp1)
  risk2 <- exp(lp2)
  overall_risk <- risk1 + risk2

  # Ensure rate is reasonable
  rate_val <- 0.1 * overall_risk
  rate_val <- pmax(rate_val, 1e-3)

  time <- rexp(n, rate = rate_val)

  # Event type selection (multinomial)
  # prob(cause1) = risk1 / (risk1 + risk2)
  p1 <- risk1 / (risk1 + risk2)
  event_type <- rbinom(n, 1, p1) + 1 # 1 or 2

  # Censoring
  cens_time <- rexp(n, rate = 0.1)
  status <- ifelse(time < cens_time, event_type, 0)
  obs_time <- pmin(time, cens_time)

  data.frame(
    time = obs_time,
    event = status,
    x1 = x1,
    x2 = x2,
    cat1 = factor(sample(c("A", "B"), n, replace = TRUE)) # Reduced levels
  )
}


n_train <- 150
n_test <- 50

train_df <- make_cr_data(n_train)
test_df <- make_cr_data(n_test)

# Create Tasks
train_df$status_bin <- ifelse(train_df$event == 0, 0, 1)
test_df$status_bin <- ifelse(test_df$event == 0, 0, 1)
train_df$cause_val <- ifelse(train_df$event == 0, NA, train_df$event)
test_df$cause_val <- ifelse(test_df$event == 0, NA, test_df$event)

train_task <- ml4t2e_task_cr(train_df, time = "time", status = "status_bin", cause = "cause_val")
test_task <- ml4t2e_task_cr(test_df, time = "time", status = "status_bin", cause = "cause_val")

# ==============================================================================
# Tests via ml4t2e_fit
# ==============================================================================

test_that("ml4t2e_fit fits Competing Risks Cox model", {
  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "cox",
    ensemble = "none"
  )

  expect_s3_class(fit, "t2e_fit")
  expect_true("cox" %in% fit$model_names)
  expect_equal(attr(fit$task, "task_type"), "competing_risk")

  # Check internal model structure
  # Now using native R6 class
  cox_obj <- fit$models$cox
  expect_true(inherits(cox_obj, "CoxCompetingRiskModel"))

  # The internal model is a list of coxph objects (one per cause)
  expect_true(is.list(cox_obj$model))
  expect_true(all(c("1", "2") %in% names(cox_obj$model)))
  expect_s3_class(cox_obj$model[["1"]], "coxph")
})

test_that("Cox CR model predictions format", {
  fit <- ml4t2e_fit(keep_data = TRUE, 
    task = train_task,
    models = "cox",
    ensemble = "none"
  )

  times <- c(1, 5)
  preds <- predict(fit, newdata = test_df, times = times, type = "cif")

  expect_s3_class(preds, "t2e_pred")
  # CR preds should have 'cause' column
  expect_true(all(c("id", "time", "cif", "model", "cause") %in% names(preds)))
  expect_equal(nrow(preds), nrow(test_df) * length(times) * 2) # 2 causes

  # Check values
  expect_true(all(preds$cif >= 0 & preds$cif <= 1))
})

test_that("Cox CR handles factor variables", {
  fit <- ml4t2e_fit(keep_data = TRUE, train_task, models = "cox")
  expect_no_error(predict(fit, newdata = test_df, times = 1))
})
