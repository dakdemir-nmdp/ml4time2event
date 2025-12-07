# ==============================================================================
# Test Suite for Time-Varying Additive Hazard (TTAH) Survival Model
# ==============================================================================

library(testthat)

# Source required files explicitly for standalone execution
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/general_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/ttah_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_ttah.R")

# ==============================================================================
# Test Data Setup
# ==============================================================================

set.seed(123)
n_train <- 180
n_test <- 40

sim_survival_data <- function(n) {
  baseline_time <- rexp(n, rate = 0.04)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n, mean = 1)
  cat_a <- factor(sample(c("A", "B", "C"), n, TRUE))
  cat_b <- factor(sample(c("Low", "Mid", "High"), n, TRUE),
                  levels = c("Low", "Mid", "High"))

  # Inject a time-varying effect: hazard increases with z1 for large times
  event_prob <- plogis(-2 + 0.5 * z1 + 0.2 * z2 + 0.3 * z1 * baseline_time)
  event <- rbinom(n, 1, pmin(0.95, pmax(0.05, event_prob)))

  censor_time <- rexp(n, rate = 0.02)
  observed_time <- pmin(baseline_time, censor_time)
  status <- ifelse(baseline_time <= censor_time, event, 0)

  data.frame(
    time = observed_time,
    event = status,
    z1 = z1,
    z2 = z2,
    z3 = z3,
    cat_a = cat_a,
    cat_b = cat_b
  )
}

train_data <- sim_survival_data(n_train)
test_data <- sim_survival_data(n_test)

expvars_numeric <- c("z1", "z2", "z3")
expvars_all <- c("z1", "z2", "z3", "cat_a", "cat_b")

# ==============================================================================
# Tests for SurvModel_TTAH - Basic Functionality
# ==============================================================================

test_that("SurvModel_TTAH fits and returns expected structure", {
  model <- SurvModel_TTAH(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    n_time = 18,
    latent_dim = 3
  )

  expect_s3_class(model, "ml4t2e_surv_ttah")
  expect_named(
    model,
    c("model", "times", "varprof", "expvars", "factor_levels",
      "time_grid", "basis_specs", "latent_projection", "engine")
  )

  expect_true(is.list(model$model))
  expect_true(is.numeric(model$times))
  expect_true(is.list(model$varprof))
  expect_equal(model$expvars, expvars_numeric)
  expect_true(length(model$time_grid) >= 5)
  expect_equal(model$engine, "base")
})

test_that("SurvModel_TTAH handles factor variables", {
  model <- SurvModel_TTAH(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    n_time = 15
  )

  expect_s3_class(model, "ml4t2e_surv_ttah")
  expect_true(is.table(model$varprof$cat_a))
  expect_true(is.table(model$varprof$cat_b))
  expect_setequal(levels(train_data$cat_a), names(model$varprof$cat_a))
  expect_setequal(levels(train_data$cat_b), names(model$varprof$cat_b))
})

# ==============================================================================
# Tests for Predict_SurvModel_TTAH - Basic Functionality
# ==============================================================================

test_that("Predict_SurvModel_TTAH returns monotone survival curves", {
  model <- SurvModel_TTAH(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    n_time = 20,
    latent_dim = 2
  )

  preds <- Predict_SurvModel_TTAH(model, test_data)

  expect_type(preds, "list")
  expect_true(is.matrix(preds$Probs))
  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_equal(preds$Probs[1, ], rep(1, nrow(test_data)))

  # Survival should be non-increasing in time
  decreasing <- apply(preds$Probs, 2, function(col) all(diff(col) <= 1e-6))
  expect_true(all(decreasing))

  # Hazards derived from survival should be between 0 and 1
  hazards <- 1 - preds$Probs[-1, , drop = FALSE] / preds$Probs[-nrow(preds$Probs), , drop = FALSE]
  expect_true(all(hazards >= -1e-6 & hazards <= 1 + 1e-6))
})

test_that("Predict_SurvModel_TTAH supports custom time grid interpolation", {
  model <- SurvModel_TTAH(
    data = train_data,
    expvars = expvars_numeric,
    timevar = "time",
    eventvar = "event",
    n_time = 12
  )

  new_times <- seq(0, max(train_data$time), length.out = 8)
  preds <- Predict_SurvModel_TTAH(model, test_data, new_times = new_times)

  expect_true(all(abs(preds$Times - new_times) < 1e-8))
  expect_equal(ncol(preds$Probs), nrow(test_data))
  expect_equal(length(preds$Times), nrow(preds$Probs))
})
