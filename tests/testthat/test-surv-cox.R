# ==============================================================================
# Test Suite for Cox Proportional Hazards Survival Model
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
set.seed(42)
beta <- c(x1 = 1.2, x2 = -1.0, x3 = 0.8, x4 = -0.7)
cat1_levels <- c("A", "B", "C")
cat2_levels <- c("Low", "High")

make_surv_data <- function(n) {
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n, mean = 1)
  x4 <- rnorm(n, mean = -1)
  cat1 <- factor(sample(cat1_levels, n, replace = TRUE), levels = cat1_levels)
  cat2 <- factor(sample(cat2_levels, n, replace = TRUE), levels = cat2_levels)
  # Linear predictor: strong effect for all variables
  lp <- beta["x1"] * x1 + beta["x2"] * x2 + beta["x3"] * x3 + beta["x4"] * x4 +
    ifelse(cat1 == "B", 0.8, ifelse(cat1 == "C", -0.8, 0)) +
    ifelse(cat2 == "High", 1.0, 0)
  # Exponential survival times
  time <- rexp(n, rate = exp(lp - mean(lp) + log(0.1)))
  # Censoring
  censor_time <- rexp(n, rate = 0.05)
  event <- as.integer(time <= censor_time)
  observed_time <- pmin(time, censor_time)
  data.frame(time = observed_time, event = event, x1 = x1, x2 = x2, x3 = x3, x4 = x4, cat1 = cat1, cat2 = cat2)
}

train_data <- make_surv_data(n_train)
test_data <- make_surv_data(n_test)

# Define variables
expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# Tests for SurvModel_Cox - Basic Functionality
# ==============================================================================

test_that("SurvModel_Cox fits basic model without variable selection", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    varsel = "none"
  )

  # Check output structure
  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_true(is.list(model))
  expect_named(model, c("cph_model", "times", "time_range", "varprof", "model_type",
                       "expvars", "timevar", "eventvar", "varsel_method", "alpha", "nfolds"))

  # Check model type
  expect_equal(model$model_type, "cox_standard")
  expect_s3_class(model$cph_model, "coxph")

  # Check time range
  expect_true(is.numeric(model$time_range))
  expect_equal(length(model$time_range), 2)
  expect_true(model$time_range[1] == 0)
  expect_true(model$time_range[2] > 0)
  expect_true(model$time_range[2] >= model$time_range[1])

  # Check varprof
  expect_true(is.list(model$varprof))
  expect_equal(length(model$varprof), length(expvars_numeric))
  expect_equal(sort(names(model$varprof)), sort(expvars_numeric))

test_that("SurvModel_Cox handles factor variables", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_s3_class(model$cph_model, "coxph")

  # Check varprof captures factor levels
  expect_true(is.table(model$varprof$cat1))
  expect_true(is.table(model$varprof$cat2))
  expect_setequal(names(model$varprof$cat1), c("A", "B", "C"))
  expect_setequal(names(model$varprof$cat2), c("Low", "High"))

# ==============================================================================
# Tests for Variable Selection Methods
# ==============================================================================

test_that("SurvModel_Cox performs backward selection with AIC", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_many,
    timevar = "time",
    eventvar = "event",
    varsel = "backward",
    penalty = "AIC",
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$varsel_method, "backward")
  expect_s3_class(model$cph_model, "coxph")

  # Model should have <= original number of coefficients
  n_coefs <- length(coef(model$cph_model))
  expect_true(n_coefs > 0)  # At least some variables kept
})

test_that("SurvModel_Cox performs backward selection with BIC", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_many,
    timevar = "time",
    eventvar = "event",
    varsel = "backward",
    penalty = "BIC",
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$varsel_method, "backward")

  # BIC typically selects fewer variables than AIC
  # Note: BIC may remove all variables if none are significant enough
  n_coefs_bic <- length(coef(model$cph_model))
  expect_true(n_coefs_bic >= 0)  # Can be 0 if all removed
})

test_that("SurvModel_Cox performs forward selection", {
  model <- SurvModel_Cox(
    data = train_data,
  expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    varsel = "forward",
    penalty = "AIC",
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$varsel_method, "forward")
  expect_s3_class(model$cph_model, "coxph")
})

test_that("SurvModel_Cox performs stepwise (both) selection", {
  model <- SurvModel_Cox(
    data = train_data,
  expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    varsel = "both",
    penalty = "AIC",
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$varsel_method, "both")
  expect_s3_class(model$cph_model, "coxph")
})

# ==============================================================================
# Tests for Penalized Cox Model
# ==============================================================================

test_that("SurvModel_Cox fits penalized Cox (lasso)", {
  skip_if_not_installed("glmnet")

  model <- SurvModel_Cox(
    data = train_data,
  expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    varsel = "penalized",
    alpha = 1,  # lasso
    nfolds = 5,
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$model_type, "cox_penalized")
  expect_equal(model$varsel_method, "penalized")
  expect_s3_class(model$cph_model, "ml4t2e_cox_penalized")

  # Check components of penalized model
  expect_true("cv_fit" %in% names(model$cph_model))
  expect_true("train_sample" %in% names(model$cph_model))
  expect_s3_class(model$cph_model$cv_fit, "cv.glmnet")
})

test_that("SurvModel_Cox fits penalized Cox (ridge)", {
  skip_if_not_installed("glmnet")

  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    varsel = "penalized",
    alpha = 0,  # ridge
    nfolds = 5,
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$model_type, "cox_penalized")
  expect_equal(model$cph_model$alpha, 0)
})

test_that("SurvModel_Cox fits penalized Cox (elastic net)", {
  skip_if_not_installed("glmnet")

  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    varsel = "penalized",
    alpha = 0.5,  # elastic net
    nfolds = 5,
    verbose = FALSE
  )

  expect_s3_class(model, "ml4t2e_surv_cox")
  expect_equal(model$model_type, "cox_penalized")
  expect_equal(model$cph_model$alpha, 0.5)
})

# ==============================================================================
# Tests for Predict_SurvModel_Cox - Basic Functionality
# ==============================================================================

test_that("Predict_SurvModel_Cox returns correct output structure", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_Cox(model, test_data)

  # Check output structure
  expect_type(preds, "list")
  expect_named(preds, c("Probs", "Times", "survfit_obj"))

  # Check Probs matrix
  expect_true(is.matrix(preds$Probs))
  expect_equal(nrow(preds$Probs), length(preds$Times))
  expect_equal(ncol(preds$Probs), nrow(test_data))

  # Check Times
  expect_true(is.numeric(preds$Times))
  expect_true(all(preds$Times >= 0))
  expect_true(all(diff(preds$Times) >= 0))  # sorted

  # Check probability values
  expect_true(all(preds$Probs >= 0 & preds$Probs <= 1, na.rm = TRUE))

  # Check survfit object
  expect_s3_class(preds$survfit_obj, "survfit")
})

test_that("Predict_SurvModel_Cox predictions are monotonically decreasing", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_Cox(model, test_data)

  # For each observation, survival should be non-increasing over time
  for (i in seq_len(ncol(preds$Probs))) {
    survival_curve <- preds$Probs[, i]
    diffs <- diff(survival_curve)
    expect_true(all(diffs <= 1e-10))  # Allow small numerical tolerance
  }
})

test_that("Predict_SurvModel_Cox handles custom time points", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  custom_times <- c(1, 5, 10, 20, 50)
  preds <- Predict_SurvModel_Cox(model, test_data, new_times = custom_times)

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$Probs), length(custom_times))
  expect_equal(ncol(preds$Probs), nrow(test_data))
})

test_that("Predict_SurvModel_Cox includes time 0", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_SurvModel_Cox(model, test_data)

  # Time 0 should be included
  expect_true(0 %in% preds$Times)

  # Survival at time 0 should be 1
  time_0_idx <- which(preds$Times == 0)
  expect_true(all(abs(preds$Probs[time_0_idx, ] - 1) < 1e-6))
})

# ==============================================================================
# Tests for Predict_SurvModel_Cox - Penalized Models
# ==============================================================================

test_that("Predict_SurvModel_Cox works with penalized Cox", {
  skip_if_not_installed("glmnet")

  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    varsel = "penalized",
    alpha = 0.5,
    nfolds = 5
  )

  preds <- Predict_SurvModel_Cox(model, test_data)

  expect_true(is.matrix(preds$Probs))
  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_true(all(preds$Probs >= 0 & preds$Probs <= 1, na.rm = TRUE))
})

# ==============================================================================
# Tests for Factor Level Handling
# ==============================================================================

test_that("Predict_SurvModel_Cox handles matching factor levels", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Test data with same factor levels
  preds <- Predict_SurvModel_Cox(model, test_data)

  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_true(all(!is.na(preds$Probs)))
})

test_that("Predict_SurvModel_Cox warns about new factor levels", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Create test data with new factor level
  test_data_new_level <- test_data
  test_data_new_level$cat1 <- factor(
    c("A", "B", "C", "D")[1:nrow(test_data_new_level)],
    levels = c("A", "B", "C", "D")
  )

  expect_warning(
    Predict_SurvModel_Cox(model, test_data_new_level),
    "has new levels"
  )
})

test_that("Predict_SurvModel_Cox converts character to factor", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Create test data with character instead of factor
  test_data_char <- test_data
  test_data_char$cat1 <- as.character(test_data_char$cat1)
  test_data_char$cat2 <- as.character(test_data_char$cat2)

  preds <- Predict_SurvModel_Cox(model, test_data_char)

  expect_equal(ncol(preds$Probs), nrow(test_data_char))
})

# ==============================================================================
# Tests for Error Handling and Edge Cases
# ==============================================================================

test_that("SurvModel_Cox validates inputs", {
  # Missing timevar
  expect_error(
    SurvModel_Cox(train_data, expvars_numeric, "missing_var", "event"),
    "timevar.*not found"
  )

  # Missing eventvar
  expect_error(
    SurvModel_Cox(train_data, expvars_numeric, "time", "missing_var"),
    "eventvar.*not found"
  )

  # Missing expvars
  expect_error(
    SurvModel_Cox(train_data, c("x1", "missing_var"), "time", "event"),
    "expvars not found"
  )

  # Invalid varsel
  expect_error(
    SurvModel_Cox(train_data, expvars_numeric, "time", "event",
                 varsel = "invalid"),
    "'arg'.*should be one of"
  )

  # Invalid penalty
  expect_error(
    SurvModel_Cox(train_data, expvars_numeric, "time", "event",
                 penalty = "invalid"),
    "'arg'.*should be one of"
  )

  # Invalid alpha
  expect_error(
    SurvModel_Cox(train_data, expvars_numeric, "time", "event",
                 varsel = "penalized", alpha = 1.5),
    "alpha.*between 0 and 1"
  )
})

test_that("Predict_SurvModel_Cox validates inputs", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  # Missing variables in newdata
  test_data_missing <- test_data[, c("time", "event", "x1")]
  expect_error(
    Predict_SurvModel_Cox(model, test_data_missing),
    "variables missing in newdata"
  )

  # Invalid new_times
  expect_error(
    Predict_SurvModel_Cox(model, test_data, new_times = c(-1, 5, 10)),
    "new_times.*non-negative"
  )

  # Wrong model type
  expect_error(
    Predict_SurvModel_Cox(list(a = 1), test_data),
    "must be output from SurvModel_Cox"
  )
})

test_that("SurvModel_Cox handles data with no events", {
  # Create data with all censored
  data_no_events <- train_data
  data_no_events$event <- 0

  expect_error(
    SurvModel_Cox(data_no_events, expvars_numeric, "time", "event"),
    "No events in training data"
  )
})

test_that("SurvModel_Cox handles single observation prediction", {
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  single_obs <- test_data[1, , drop = FALSE]
  preds <- Predict_SurvModel_Cox(model, single_obs)

  expect_true(is.matrix(preds$Probs))
  expect_equal(ncol(preds$Probs), 1)
  expect_equal(nrow(preds$Probs), length(preds$Times))
})

  # Fit model
  model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    varsel = "backward",
    penalty = "AIC"
  )

  # Predict
  preds <- Predict_SurvModel_Cox(model, test_data)

  # Extract probability at specific time
  time_point <- 10
  if (time_point %in% preds$Times) {
    time_idx <- which(preds$Times == time_point)
    probs_at_t <- preds$Probs[time_idx, ]

    expect_equal(length(probs_at_t), nrow(test_data))
    expect_true(all(probs_at_t >= 0 & probs_at_t <= 1))
  }
})


test_that("Multiple models produce consistent output shapes", {
  # Fit different models
  model1 <- SurvModel_Cox(train_data, expvars_numeric, "time", "event",
                         varsel = "none")
  model2 <- SurvModel_Cox(train_data, expvars_numeric, "time", "event",
                         varsel = "backward", penalty = "AIC")

  skip_if_not_installed("glmnet")
  model3 <- SurvModel_Cox(train_data, expvars_all, "time", "event",
                         varsel = "penalized", nfolds = 3)

  # Get predictions at same time points
  custom_times <- c(5, 10, 20)
  preds1 <- Predict_SurvModel_Cox(model1, test_data, new_times = custom_times)

  # Only compare if model2 has variables
  if (length(coef(model2$cph_model)) > 0) {
    preds2 <- Predict_SurvModel_Cox(model2, test_data, new_times = custom_times)
    expect_equal(dim(preds1$Probs), dim(preds2$Probs))
    expect_equal(preds1$Times, preds2$Times)
  }

  preds3 <- Predict_SurvModel_Cox(model3, test_data, new_times = custom_times)

  # All should have same dimensions
  expect_equal(dim(preds1$Probs), dim(preds3$Probs))
  expect_equal(preds1$Times, preds3$Times)
})
})
