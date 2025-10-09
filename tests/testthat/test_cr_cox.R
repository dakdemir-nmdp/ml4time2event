# ==============================================================================
# Test Suite for Competing Risks Cox Model
# ==============================================================================

library(testthat)
library(survival)

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create simulated competing risks data for testing
set.seed(42)
n_train <- 50  # Reduced for faster testing
n_test <- 20   # Reduced for faster testing

# Training data - competing risks format
train_data <- data.frame(
  time = rexp(n_train, rate = 0.1),
  event = sample(0:2, n_train, replace = TRUE, prob = c(0.3, 0.4, 0.3)), # 0=censored, 1=cause1, 2=cause2
  x1 = rnorm(n_train),
  x2 = rnorm(n_train),
  x3 = rnorm(n_train, mean = 1),
  x4 = rnorm(n_train, mean = -1),
  cat1 = factor(sample(c("A", "B", "C"), n_train, replace = TRUE)),
  cat2 = factor(sample(c("Low", "High"), n_train, replace = TRUE))
)

# Test data
test_data <- data.frame(
  time = rexp(n_test, rate = 0.1),
  event = sample(0:2, n_test, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
  x1 = rnorm(n_test),
  x2 = rnorm(n_test),
  x3 = rnorm(n_test, mean = 1),
  x4 = rnorm(n_test, mean = -1),
  cat1 = factor(sample(c("A", "B", "C"), n_test, replace = TRUE),
                levels = c("A", "B", "C")),
  cat2 = factor(sample(c("Low", "High"), n_test, replace = TRUE),
                levels = c("Low", "High"))
)

# Define variables
expvars_numeric <- c("x1", "x2", "x3")
expvars_all <- c("x1", "x2", "x3", "cat1", "cat2")
expvars_many <- c("x1", "x2", "x3", "x4", "cat1", "cat2")

# ==============================================================================
# Tests for CRModel_Cox - Basic Functionality
# ==============================================================================

test_that("CRModel_Cox fits basic model", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  # Check output structure
  expect_type(model, "list")
  expect_s3_class(model, "ml4t2e_cr_cox")
  expect_setequal(names(model),
                  c("cph_models_all_causes", "times", "varprof", "model_type",
                    "expvars", "timevar", "eventvar", "event_codes",
                    "varsel_method", "time_range"))

  # Check model type
  expect_equal(model$model_type, "cr_cox")

  # Ensure we have a model per event code
  non_null_models <- Filter(Negate(is.null), model$cph_models_all_causes)
  expect_true(length(non_null_models) >= 1)
  expect_true(all(names(non_null_models) %in% model$event_codes))
  expect_true(all(vapply(non_null_models, inherits, logical(1), what = "coxph")))

  # Check time range
  expect_true(is.numeric(model$time_range))
  expect_equal(length(model$time_range), 2)
  expect_true(all(is.finite(model$time_range)))

  # Check times
  expect_true(is.numeric(model$times))
  expect_true(length(model$times) > 0)
  expect_true(all(model$times >= 0))

  # Check varprof
  expect_true(is.list(model$varprof))
  expect_equal(length(model$varprof), length(expvars_numeric))
  expect_true(all(model$event_codes %in% c("1", "2")))
})

test_that("CRModel_Cox handles factor variables", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  expect_s3_class(model, "ml4t2e_cr_cox")
  expect_true(is.list(model$varprof))
  expect_true(all(names(Filter(Negate(is.null), model$cph_models_all_causes)) %in% model$event_codes))
})

# ==============================================================================
# Tests for Predict_CRModel_Cox - Basic Functionality
# ==============================================================================

test_that("Predict_CRModel_Cox returns correct output structure", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_CRModel_Cox(model, test_data, event_of_interest = "1")

  # Check output structure (no longer returns cph_modelTestPredict)
  expect_type(preds, "list")
  expect_named(preds, c("CIFs", "Times"))

  # Check CIFs matrix
  expect_true(is.matrix(preds$CIFs))
  expect_equal(ncol(preds$CIFs), nrow(test_data))

  # Check Times
  expect_true(is.numeric(preds$Times))
  expect_true(all(preds$Times >= 0))

  # Check CIF values (should be between 0 and 1)
  expect_true(all(preds$CIFs >= 0 & preds$CIFs <= 1, na.rm = TRUE))
})

test_that("Predict_CRModel_Cox CIFs are monotonically non-decreasing", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_CRModel_Cox(model, test_data, event_of_interest = "1")

  # For each observation, CIF should be non-decreasing over time
  for (i in seq_len(ncol(preds$CIFs))) {
    cif_curve <- preds$CIFs[, i]
    diffs <- diff(cif_curve)
    expect_true(all(diffs >= -1e-10))  # Allow small numerical tolerance
  }
})

test_that("Predict_CRModel_Cox handles custom time points", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  custom_times <- c(1, 5, 10, 20, 50)
  preds <- Predict_CRModel_Cox(model, test_data, newtimes = custom_times, event_of_interest = "1")

  expect_equal(preds$Times, custom_times)
  expect_equal(nrow(preds$CIFs), length(custom_times))
  expect_equal(ncol(preds$CIFs), nrow(test_data))
})

test_that("Predict_CRModel_Cox validates requested event", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  expect_error(
    Predict_CRModel_Cox(model, test_data, event_of_interest = "999"),
    "Requested event_of_interest"
  )
})

test_that("Predict_CRModel_Cox includes time 0", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  preds <- Predict_CRModel_Cox(model, test_data, event_of_interest = "1")

  # Time 0 should be included
  expect_true(0 %in% preds$Times)

  # CIF at time 0 should be 0
  time_0_idx <- which(preds$Times == 0)
  expect_true(all(abs(preds$CIFs[time_0_idx, ]) < 1e-6))
})

# ==============================================================================
# Tests for Factor Level Handling
# ==============================================================================

test_that("Predict_CRModel_Cox handles matching factor levels", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Test data with same factor levels
  preds <- Predict_CRModel_Cox(model, test_data, event_of_interest = "1")

  expect_equal(ncol(preds$CIFs), nrow(test_data))
  expect_true(all(!is.na(preds$CIFs)))
})

test_that("Predict_CRModel_Cox warns about new factor levels", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Create test data with new factor level
  test_data_new_level <- test_data
  test_data_new_level$cat1 <- factor(
    c("A", "B", "C", "D")[seq_len(nrow(test_data_new_level))],
    levels = c("A", "B", "C", "D")
  )

  expect_warning(
    preds_new_level <- Predict_CRModel_Cox(model, test_data_new_level, event_of_interest = "1"),
    "has new levels"
  )
  expect_true(any(colSums(is.na(preds_new_level$CIFs)) > 0))
})

test_that("Predict_CRModel_Cox converts character to factor", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event"
  )

  # Create test data with character instead of factor
  test_data_char <- test_data
  test_data_char$cat1 <- as.character(test_data_char$cat1)
  test_data_char$cat2 <- as.character(test_data_char$cat2)

  preds <- Predict_CRModel_Cox(model, test_data_char, event_of_interest = "1")

  expect_equal(ncol(preds$CIFs), nrow(test_data_char))
  expect_true(all(!is.na(preds$CIFs)))
})

# ==============================================================================
# Tests for Error Handling
# ==============================================================================

test_that("CRModel_Cox handles missing data appropriately", {
  # Create data with missing values
  train_data_na <- train_data
  train_data_na$x1[1:10] <- NA

  # Should handle missing data by removing rows
  model <- CRModel_Cox(
    data = train_data_na,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  expect_s3_class(model, "ml4t2e_cr_cox")
  fitted_models <- Filter(Negate(is.null), model$cph_models_all_causes)
  expect_true(length(fitted_models) >= 1)
  expect_true(all(vapply(fitted_models, inherits, logical(1), what = "coxph")))
})

test_that("CRModel_Cox requires valid inputs", {
  expect_error(CRModel_Cox(data = train_data, expvars = "nonexistent_var",
                          timevar = "time", eventvar = "event"))
  expect_error(CRModel_Cox(data = train_data, expvars = expvars_numeric,
                          timevar = "nonexistent_time", eventvar = "event"))
  expect_error(CRModel_Cox(data = train_data, expvars = expvars_numeric,
                          timevar = "time", eventvar = "nonexistent_event"))
})

test_that("Predict_CRModel_Cox requires valid inputs", {
  model <- CRModel_Cox(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event"
  )

  expect_error(Predict_CRModel_Cox(model, test_data[, !names(test_data) %in% "x1"], event_of_interest = "1"))  # Missing required variable (x1)
  expect_error(Predict_CRModel_Cox(list(), test_data, event_of_interest = "1"))  # Invalid model
})

test_that("CRModel_Cox requires data, expvars, timevar, eventvar", {
  skip_if_not_installed("survival")
  expect_error(CRModel_Cox(expvars = expvars_numeric, timevar = "time", eventvar = "event"), "argument \"data\" is missing")
  expect_error(CRModel_Cox(data = train_data, timevar = "time", eventvar = "event"), "argument \"expvars\" is missing")
  expect_error(CRModel_Cox(data = train_data, expvars = expvars_numeric, eventvar = "event"), "argument \"timevar\" is missing")
  expect_error(CRModel_Cox(data = train_data, expvars = expvars_numeric, timevar = "time"), "argument \"eventvar\" is missing")
})


# --- Tests for Predict_CRModel_Cox ---

test_that("Predict_CRModel_Cox returns predictions in correct format", {
  skip_if_not_installed("survival")

  # Fit the model first
  model_cox <- CRModel_Cox(data = train_data,
                          expvars = expvars_numeric,
                          timevar = "time",
                          eventvar = "event")

  # Get predictions
  predictions <- Predict_CRModel_Cox(modelout = model_cox, newdata = test_data, event_of_interest = "1")

  # Check output structure (no longer returns cph_modelTestPredict)
  expect_type(predictions, "list")
  expect_named(predictions, c("CIFs", "Times"))

  # Check CIFs matrix
  expect_true(is.matrix(predictions$CIFs))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))  # columns = observations
  expect_equal(nrow(predictions$CIFs), length(predictions$Times))  # rows = times

  # Check times
  expect_true(is.numeric(predictions$Times))
  expect_true(all(predictions$Times >= 0))
  expect_true(predictions$Times[1] == 0)  # Should start at 0

  # Check CIF values are probabilities
  expect_true(all(predictions$CIFs >= 0 & predictions$CIFs <= 1, na.rm = TRUE))

  # Check monotonicity (CIF should be non-decreasing)
  expect_true(all(diff(predictions$CIFs[, 1]) >= 0))  # Check first observation
})

test_that("Predict_CRModel_Cox handles custom time points", {
  skip_if_not_installed("survival")

  model_cox <- CRModel_Cox(data = train_data,
                          expvars = expvars_numeric,
                          timevar = "time",
                          eventvar = "event")

  custom_times <- c(1, 5, 10)
  predictions <- Predict_CRModel_Cox(modelout = model_cox, newdata = test_data,
                                     newtimes = custom_times, event_of_interest = "1")

  expect_type(predictions, "list")
  expect_true(is.matrix(predictions$CIFs))
  expect_equal(length(predictions$Times), length(custom_times))
  expect_equal(nrow(predictions$CIFs), length(custom_times))
  expect_equal(ncol(predictions$CIFs), nrow(test_data))
})

test_that("Predict_CRModel_Cox requires model and data", {
  skip_if_not_installed("survival")

  model_cox <- CRModel_Cox(data = train_data,
                          expvars = expvars_numeric,
                          timevar = "time",
                          eventvar = "event")

  expect_error(Predict_CRModel_Cox(modelout = model_cox), "argument \"newdata\" is missing")
  expect_error(Predict_CRModel_Cox(newdata = test_data), "argument \"modelout\" is missing")
})
