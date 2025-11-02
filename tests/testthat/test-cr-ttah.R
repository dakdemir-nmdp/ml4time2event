# ==============================================================================
# Test Suite for Time-Varying Additive Hazard (TTAH) Competing Risks Model
# ==============================================================================

library(testthat)

source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/general_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/ttah_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_ttah.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_ttah.R")

set.seed(222)
n_train <- 240
n_test <- 60

sim_competing_risks <- function(n) {
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- runif(n, -1, 1)
  cat_a <- factor(sample(c("A", "B"), n, TRUE))
  cat_b <- factor(sample(c("Low", "Medium", "High"), n, TRUE),
                  levels = c("Low", "Medium", "High"))

  # Cause-specific latent times with covariate effects
  rate1 <- exp(-0.3 + 0.4 * x1 - 0.2 * x2)
  rate2 <- exp(-0.4 + 0.3 * x2 + 0.2 * x3)
  t1 <- rexp(n, rate = rate1 / 20)
  t2 <- rexp(n, rate = rate2 / 15)
  censor <- rexp(n, rate = 1 / 30)

  time <- pmin(t1, t2, censor)
  event <- ifelse(time == t1 & time < censor, 1,
                  ifelse(time == t2 & time < censor, 2, 0))

  data.frame(
    time = time,
    event = event,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    cat_a = cat_a,
    cat_b = cat_b
  )
}

train_data <- sim_competing_risks(n_train)
test_data <- sim_competing_risks(n_test)

expvars_all <- c("x1", "x2", "x3", "cat_a", "cat_b")

test_that("CRModel_TTAH fits and persists metadata", {
  model <- CRModel_TTAH(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    n_time = 16,
    latent_dim = 2,
    maxit = 200
  )

  expect_s3_class(model, "ml4t2e_cr_ttah")
  expect_named(
    model,
    c("model", "times", "varprof", "expvars", "factor_levels", "time_grid",
      "basis_specs", "latent_projection", "cause_codes", "engine")
  )
  expect_true(is.list(model$model))
  expect_true(is.numeric(model$times))
  expect_equal(model$engine, "base")
  expect_true(length(model$cause_codes) >= 2)
})

test_that("Predict_CRModel_TTAH returns coherent hazards and CIFs", {
  model <- CRModel_TTAH(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    n_time = 18,
    latent_dim = 3
  )

  preds <- Predict_CRModel_TTAH(model, test_data)

  expect_type(preds, "list")
  expect_true(all(c("CauseSpecificHazard", "CauseSpecificCIF",
                    "TotalSurvival", "Times", "AllCauseHazard") %in% names(preds)))

  haz_df <- preds$CauseSpecificHazard
  cif_df <- preds$CauseSpecificCIF
  surv_df <- preds$TotalSurvival

  expect_true(is.matrix(haz_df))
  expect_equal(ncol(haz_df), nrow(test_data))
  expect_true(all(haz_df >= -1e-6 & haz_df <= 1 + 1e-6))

  expect_true(is.matrix(cif_df))
  expect_equal(ncol(cif_df), nrow(test_data))
  expect_true(all(diff(preds$Times) >= 0))

  # CIFs should be monotone increasing and bounded
  increasing <- apply(cif_df, 2, function(col) all(diff(col) >= -1e-6))
  expect_true(all(increasing))
  expect_true(all(cif_df >= -1e-6 & cif_df <= 1 + 1e-6))

  # Total survival should start at 1 and be non-increasing
  expect_equal(surv_df[1, ], rep(1, nrow(test_data)))
  decreasing_surv <- apply(surv_df, 2, function(col) all(diff(col) <= 1e-6))
  expect_true(all(decreasing_surv))

  # Stay probabilities plus total hazard mass should sum to 1
  stay <- preds$StayProbability
  hazard_total <- apply(preds$AllCauseHazard, c(1, 2), sum)
  combined <- stay + hazard_total
  expect_true(all(abs(combined - 1) < 1e-5))
})

test_that("Predict_CRModel_TTAH handles custom event code and time grid", {
  model <- CRModel_TTAH(
    data = train_data,
    expvars = expvars_all,
    timevar = "time",
    eventvar = "event",
    n_time = 12
  )

  new_times <- seq(0, max(train_data$time), length.out = 7)
  preds_cause2 <- Predict_CRModel_TTAH(
    model, test_data,
    new_times = new_times,
    event_of_interest = "2"
  )

  expect_equal(length(preds_cause2$Times), length(new_times))
  expect_equal(nrow(preds_cause2$CauseSpecificHazard), length(new_times) - 1)
  expect_true(all(preds_cause2$CauseSpecificCIF >= -1e-6))
  expect_true(all(preds_cause2$CauseSpecificCIF <= 1 + 1e-6))
})
