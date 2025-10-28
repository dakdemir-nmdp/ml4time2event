# Comprehensive ETL Validation Tests for Competing Risks Models
# Purpose: Diagnose potential issues with reversed age effects in SHAP analysis

library(testthat)
library(ml4time2event)

# ==============================================================================
# Test 1: Matrix Orientation Verification
# ==============================================================================

test_that("CR predictions have correct matrix orientation (rows=times, cols=obs)", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE,
    prediction_times = seq(0, 100, length.out = 20)
  )

  # Get predictions directly
  test_data <- bmt_small[1:5, ]
  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = test_data,
    new_times = seq(0, 100, length.out = 20),
    ensemble_method = "average"
  )

  # Check matrix dimensions
  expect_true(!is.null(preds$NewProbs))
  expect_equal(nrow(preds$NewProbs), 20)  # 20 time points
  expect_equal(ncol(preds$NewProbs), 5)   # 5 observations

  message("Matrix orientation: ", nrow(preds$NewProbs), " times x ",
          ncol(preds$NewProbs), " observations")
})

# ==============================================================================
# Test 2: CIF Values Logical Checks
# ==============================================================================

test_that("CR CIF values are non-negative", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = bmt_small[1:10, ],
    new_times = seq(0, 100, length.out = 20),
    ensemble_method = "average"
  )

  # Check all CIF values are >= 0
  expect_true(all(preds$NewProbs >= 0, na.rm = TRUE))

  message("Min CIF value: ", min(preds$NewProbs, na.rm = TRUE))
  message("Max CIF value: ", max(preds$NewProbs, na.rm = TRUE))
})

test_that("CR CIF values are bounded by [0, 1]", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = bmt_small[1:10, ],
    new_times = seq(0, 100, length.out = 20),
    ensemble_method = "average"
  )

  # Check all CIF values are <= 1
  expect_true(all(preds$NewProbs <= 1.0, na.rm = TRUE))

  # Check CIF at t=0 is close to 0
  expect_true(preds$NewProbs[1, 1] < 0.01)
})

test_that("CR CIF values are monotonically non-decreasing over time", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = bmt_small[1:10, ],
    new_times = seq(0, 100, length.out = 20),
    ensemble_method = "average"
  )

  # Check each observation's CIF is non-decreasing
  for (j in 1:ncol(preds$NewProbs)) {
    cif_obs <- preds$NewProbs[, j]
    diffs <- diff(cif_obs)
    # Allow small numerical errors (1e-6)
    expect_true(all(diffs >= -1e-6, na.rm = TRUE),
                info = paste("Observation", j, "has decreasing CIF"))
  }
})

# ==============================================================================
# Test 3: ETL Calculation Verification
# ==============================================================================

test_that("ETL calculation produces sensible values for CR models", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Create prediction function
  time_horizon <- 100
  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = time_horizon)

  # Get ETL predictions
  test_data <- bmt_small[1:10, ]
  etl_values <- pred_fn(test_data)

  # ETL should be between 0 and time_horizon
  expect_true(all(etl_values >= 0))
  expect_true(all(etl_values <= time_horizon))
  expect_true(all(is.finite(etl_values)))

  message("ETL range: [", min(etl_values), ", ", max(etl_values), "]")
  message("Mean ETL: ", mean(etl_values))
})

test_that("Manual ETL calculation matches prediction function", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("pracma")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Get one observation
  test_obs <- bmt_small[1, , drop = FALSE]
  time_horizon <- 100

  # Method 1: Use prediction function
  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = time_horizon)
  etl_from_fn <- pred_fn(test_obs)

  # Method 2: Manual calculation
  time_points <- seq(0, time_horizon, length.out = 100)
  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = test_obs,
    new_times = time_points,
    ensemble_method = "average"
  )

  cif_values <- preds$NewProbs[, 1]
  etl_manual <- pracma::trapz(time_points, cif_values)

  # Should match closely (within 1% relative error)
  rel_error <- abs(etl_from_fn - etl_manual) / (etl_manual + 1e-10)
  expect_true(rel_error < 0.01,
              info = paste("ETL mismatch: fn=", etl_from_fn,
                          "manual=", etl_manual,
                          "rel_error=", rel_error))

  message("ETL from function: ", etl_from_fn)
  message("ETL manual: ", etl_manual)
  message("Relative error: ", rel_error)
})

# ==============================================================================
# Test 4: Age Effect Direction Verification
# ==============================================================================

test_that("Higher age leads to higher ETL (worse outcome) for CR models", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()

  # Create two synthetic patients: identical except for age
  base_patient <- bmt_df[1, ]

  young_patient <- base_patient
  young_patient$age <- 20

  old_patient <- base_patient
  old_patient$age <- 60

  test_data <- rbind(young_patient, old_patient)

  # Fit pipeline on full data
  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_df[1:50, ],
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Get ETL predictions
  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = 100)
  etl_values <- pred_fn(test_data)

  etl_young <- etl_values[1]
  etl_old <- etl_values[2]

  message("ETL for age=20: ", etl_young)
  message("ETL for age=60: ", etl_old)
  message("Difference (old - young): ", etl_old - etl_young)

  # Clinical expectation: older age → worse outcome → higher ETL
  # This test will FAIL if the relationship is reversed
  expect_true(etl_old > etl_young,
              info = paste("Expected higher ETL for older age, but got: young=",
                          etl_young, "old=", etl_old))
})

test_that("Age SHAP values have expected sign for CR models", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:40, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Calculate SHAP for subset with age variation
  explain_data <- bmt_small[1:20, ]
  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = explain_data,
    time_horizon = 100,
    nsim = 30
  )

  # Get age SHAP values and age values
  age_shap <- shap_result$shap_values[, "age"]
  age_values <- explain_data$age

  # Check correlation between age and age SHAP
  # Positive SHAP = increases ETL = worse outcome
  # If higher age → higher SHAP, that makes clinical sense
  cor_age_shap <- cor(age_values, age_shap, use = "complete.obs")

  message("Correlation between age and age SHAP: ", cor_age_shap)
  message("Mean age SHAP: ", mean(age_shap))
  message("Age range: [", min(age_values), ", ", max(age_values), "]")

  # Clinical expectation: positive correlation
  # (higher age → higher SHAP → higher ETL → worse outcome)
  expect_true(cor_age_shap > 0,
              info = paste("Expected positive correlation, got:", cor_age_shap))
})

# ==============================================================================
# Test 5: Comparison with Survival Models
# ==============================================================================

test_that("Survival model ETL calculation works as expected", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:40, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  # Create two patients with different ages
  base_patient <- lung_small[1, ]

  young_patient <- base_patient
  young_patient$age <- 40

  old_patient <- base_patient
  old_patient$age <- 75

  test_data <- rbind(young_patient, old_patient)

  # Get ETL predictions
  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = 365)
  etl_values <- pred_fn(test_data)

  etl_young <- etl_values[1]
  etl_old <- etl_values[2]

  message("Survival ETL for age=40: ", etl_young)
  message("Survival ETL for age=75: ", etl_old)
  message("Difference (old - young): ", etl_old - etl_young)

  # For lung cancer: older age → worse survival → higher ETL
  expect_true(etl_old > etl_young,
              info = "Survival model age effect check")
})

# ==============================================================================
# Test 6: Event of Interest Verification
# ==============================================================================

test_that("CR models predict the correct event of interest", {
  bmt_df <- get_bmt_competing_risks_data()

  # Check what events exist
  event_codes <- unique(bmt_df$status)
  message("Event codes in BMT data: ", paste(event_codes, collapse = ", "))
  message("Event counts: ", paste(table(bmt_df$status), collapse = ", "))

  # Fit pipeline and check what event it's predicting
  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_df[1:50, ],
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Check pipeline metadata
  expect_true(!is.null(pipeline$model))

  # The event codes should be stored
  model_output <- pipeline$model
  if ("Cox_Model" %in% names(model_output)) {
    cox_model <- model_output$Cox_Model
    message("Event codes in model: ", paste(cox_model$event_codes, collapse = ", "))
  }
})

# ==============================================================================
# Test 7: SHAP Additivity for CR Models
# ==============================================================================

test_that("SHAP additivity holds for CR models (prediction = baseline + sum(SHAP))", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = bmt_small[1:10, ],
    time_horizon = 100,
    nsim = 30
  )

  # Check additivity: prediction = baseline + sum(SHAP)
  shap_sums <- rowSums(shap_result$shap_values)
  expected <- shap_result$predictions - shap_result$baseline

  rel_errors <- abs(shap_sums - expected) / (abs(expected) + 1e-10)

  message("Max relative error in SHAP additivity: ", max(rel_errors))

  # Additivity should hold within numerical precision
  expect_true(all(rel_errors < 1e-5),
              info = paste("SHAP additivity violated, max error:", max(rel_errors)))
})

# ==============================================================================
# Test 8: CIF vs Time Relationship
# ==============================================================================

test_that("CIF increases with time for all patients", {
  skip_if_not_installed("fastshap")

  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:20, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE
  )

  # Get predictions at multiple time points
  time_points <- c(0, 25, 50, 75, 100)
  test_obs <- bmt_small[1, , drop = FALSE]

  preds <- PredictCRModels(
    models = pipeline$model,
    newdata = test_obs,
    new_times = time_points,
    ensemble_method = "average"
  )

  cif_values <- preds$NewProbs[, 1]

  message("CIF at times ", paste(time_points, collapse = ", "), ":")
  message(paste(round(cif_values, 4), collapse = ", "))

  # Each subsequent CIF should be >= previous (monotonic)
  for (i in 2:length(cif_values)) {
    expect_true(cif_values[i] >= cif_values[i-1] - 1e-6,
                info = paste("CIF not monotonic at time", time_points[i]))
  }
})
