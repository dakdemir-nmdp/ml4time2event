# test-shap-analysis.R
# Test suite for SHAP-based explainability features

# Setup: Load required libraries and create test fixtures
# We'll use the lung survival and BMT competing risks datasets

library(testthat)
library(ml4time2event)

# ============================================================================
# Phase 2: TDD Cycle 1 - Prediction Wrapper Tests (RED)
# ============================================================================

test_that("ml4t2e_shap_predict_fn: function exists", {
  expect_true(exists("ml4t2e_shap_predict_fn"))
  expect_true(is.function(ml4t2e_shap_predict_fn))
})

test_that("ml4t2e_shap_predict_fn: works with survival pipeline", {
  skip_if_not_installed("fastshap")

  # Create a minimal survival pipeline
  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:50, ]  # Small dataset for speed

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE,
    prediction_times = seq(0, 1000, length.out = 30)
  )

  # Create prediction function
  pred_fn <- ml4t2e_shap_predict_fn(
    pipeline = pipeline,
    time_horizon = 365  # 1-year expected time lost
  )

  # Test that it returns a function
  expect_true(is.function(pred_fn))

  # Test that the function works on raw data
  test_data <- lung_small[1:5, ]
  predictions <- pred_fn(test_data)

  # Should return numeric vector
  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), nrow(test_data))

  # All predictions should be non-negative (time lost >= 0)
  expect_true(all(predictions >= 0))

  # Predictions should be finite
  expect_true(all(is.finite(predictions)))
})

test_that("ml4t2e_shap_predict_fn: works with competing risks pipeline", {
  skip_if_not_installed("fastshap")

  # Create a minimal CR pipeline
  bmt_df <- get_bmt_competing_risks_data()
  bmt_small <- bmt_df[1:50, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = bmt_small,
    analysis_type = "competing_risks",
    timevar = "ftime",
    eventvar = "status",
    models = c("cox"),
    include_rf = FALSE,
    prediction_times = seq(0, 150, length.out = 30)
  )

  pred_fn <- ml4t2e_shap_predict_fn(
    pipeline = pipeline,
    time_horizon = 100
  )

  expect_true(is.function(pred_fn))

  test_data <- bmt_small[1:5, ]
  predictions <- pred_fn(test_data)

  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), nrow(test_data))
  expect_true(all(predictions >= 0))
  expect_true(all(is.finite(predictions)))
})

test_that("ml4t2e_shap_predict_fn: handles preprocessing correctly", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:50, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = 365)

  # The prediction function should work on raw data (before preprocessing)
  # Get the original variable names that went into the pipeline
  raw_data <- lung_small[1:3, ]

  predictions <- pred_fn(raw_data)

  # Should successfully process raw data through the pipeline's recipe
  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), 3)
})

test_that("ml4t2e_shap_predict_fn: validates inputs correctly", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  # Should error with invalid pipeline
  expect_error(
    ml4t2e_shap_predict_fn(pipeline = "not_a_pipeline", time_horizon = 365),
    "pipeline.*ml4t2e_pipeline"
  )

  # Should error with invalid time_horizon
  expect_error(
    ml4t2e_shap_predict_fn(pipeline = pipeline, time_horizon = -1),
    "time_horizon.*positive"
  )

  expect_error(
    ml4t2e_shap_predict_fn(pipeline = pipeline, time_horizon = NA),
    "time_horizon"
  )
})

test_that("ml4t2e_shap_predict_fn: returns consistent predictions", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = 365)

  test_data <- lung_small[1:5, ]

  # Calling the same function twice should give same results
  preds1 <- pred_fn(test_data)
  preds2 <- pred_fn(test_data)

  expect_equal(preds1, preds2)
})

# ============================================================================
# Phase 3: TDD Cycle 2 - SHAP Calculation Tests (RED)
# ============================================================================

test_that("ml4t2e_calculate_shap: function exists", {
  expect_true(exists("ml4t2e_calculate_shap"))
  expect_true(is.function(ml4t2e_calculate_shap))
})

test_that("ml4t2e_calculate_shap: works with survival pipeline", {
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

  # Calculate SHAP values for a subset
  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10  # Small nsim for speed in tests
  )

  # Should return a list/object with SHAP values
  expect_true(!is.null(shap_result))

  # Should contain SHAP values matrix
  expect_true("shap_values" %in% names(shap_result))

  # SHAP values should be a matrix or data.frame
  expect_true(is.matrix(shap_result$shap_values) || is.data.frame(shap_result$shap_values))

  # Should have correct dimensions (rows = observations, cols = features)
  expect_equal(nrow(shap_result$shap_values), 10)

  # Should have feature names
  expect_true(length(colnames(shap_result$shap_values)) > 0)
})

test_that("ml4t2e_calculate_shap: works with competing risks pipeline", {
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

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = bmt_small[1:10, ],
    time_horizon = 100,
    nsim = 10
  )

  expect_true(!is.null(shap_result))
  expect_true("shap_values" %in% names(shap_result))
  expect_true(is.matrix(shap_result$shap_values) || is.data.frame(shap_result$shap_values))
  expect_equal(nrow(shap_result$shap_values), 10)
})

test_that("ml4t2e_calculate_shap: includes baseline prediction", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:8, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should include baseline (expected value)
  expect_true("baseline" %in% names(shap_result))
  expect_true(is.numeric(shap_result$baseline))
  expect_equal(length(shap_result$baseline), 1)
})

test_that("ml4t2e_calculate_shap: includes predictions", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:8, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should include actual predictions
  expect_true("predictions" %in% names(shap_result))
  expect_true(is.numeric(shap_result$predictions))
  expect_equal(length(shap_result$predictions), 8)
})

test_that("ml4t2e_calculate_shap: SHAP values sum to prediction - baseline", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:5, ],
    time_horizon = 365,
    nsim = 10
  )

  # SHAP values should approximately sum to (prediction - baseline) for each observation
  shap_sums <- rowSums(shap_result$shap_values)
  expected <- shap_result$predictions - shap_result$baseline

  # Allow for small numerical differences
  expect_true(all(abs(shap_sums - expected) < 1e-6))
})

test_that("ml4t2e_calculate_shap: validates inputs", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  # Invalid pipeline
  expect_error(
    ml4t2e_calculate_shap(
      pipeline = "not_a_pipeline",
      data = lung_small[1:5, ],
      time_horizon = 365
    ),
    "pipeline.*ml4t2e_pipeline"
  )

  # Invalid data
  expect_error(
    ml4t2e_calculate_shap(
      pipeline = pipeline,
      data = "not_a_dataframe",
      time_horizon = 365
    ),
    "data.*data.frame"
  )

  # Invalid time_horizon
  expect_error(
    ml4t2e_calculate_shap(
      pipeline = pipeline,
      data = lung_small[1:5, ],
      time_horizon = -100
    ),
    "time_horizon.*positive"
  )
})

test_that("ml4t2e_calculate_shap: respects nsim parameter", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  # Should accept nsim parameter
  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:5, ],
    time_horizon = 365,
    nsim = 5
  )

  expect_true(!is.null(shap_result))
  expect_equal(nrow(shap_result$shap_values), 5)
})

# ============================================================================
# Phase 4: TDD Cycle 3 - Variable Importance Plot Tests (RED)
# ============================================================================

test_that("ml4t2e_shap_importance: function exists", {
  expect_true(exists("ml4t2e_shap_importance"))
  expect_true(is.function(ml4t2e_shap_importance))
})

test_that("ml4t2e_shap_importance: returns a ggplot object", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("ggplot2")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  plot <- ml4t2e_shap_importance(shap_result)

  expect_s3_class(plot, "ggplot")
})

test_that("ml4t2e_shap_importance: handles max_features parameter", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("ggplot2")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should work with max_features specified
  plot <- ml4t2e_shap_importance(shap_result, max_features = 5)

  expect_s3_class(plot, "ggplot")
})

test_that("ml4t2e_shap_importance: validates inputs", {
  expect_error(
    ml4t2e_shap_importance("not_shap_result"),
    "shap_result"
  )
})

# ============================================================================
# Phase 5: TDD Cycle 4 - Dependence Plot Tests (RED)
# ============================================================================

test_that("ml4t2e_shap_dependence: function exists", {
  expect_true(exists("ml4t2e_shap_dependence"))
  expect_true(is.function(ml4t2e_shap_dependence))
})

test_that("ml4t2e_shap_dependence: returns a ggplot object", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("ggplot2")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  # Get a valid feature name
  feature_name <- colnames(shap_result$shap_values)[1]

  plot <- ml4t2e_shap_dependence(shap_result, feature = feature_name)

  expect_s3_class(plot, "ggplot")
})

test_that("ml4t2e_shap_dependence: validates feature parameter", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should error with invalid feature name
  expect_error(
    ml4t2e_shap_dependence(shap_result, feature = "nonexistent_feature"),
    "feature.*not found"
  )
})

# ============================================================================
# Phase 6: TDD Cycle 5 - Waterfall Plot Tests (RED)
# ============================================================================

test_that("ml4t2e_shap_waterfall: function exists", {
  expect_true(exists("ml4t2e_shap_waterfall"))
  expect_true(is.function(ml4t2e_shap_waterfall))
})

test_that("ml4t2e_shap_waterfall: returns a ggplot object", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("ggplot2")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  plot <- ml4t2e_shap_waterfall(shap_result, obs_id = 1)

  expect_s3_class(plot, "ggplot")
})

test_that("ml4t2e_shap_waterfall: validates obs_id parameter", {
  skip_if_not_installed("fastshap")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:5, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should error with invalid obs_id
  expect_error(
    ml4t2e_shap_waterfall(shap_result, obs_id = 100),
    "obs_id.*out of range"
  )

  expect_error(
    ml4t2e_shap_waterfall(shap_result, obs_id = 0),
    "obs_id.*positive"
  )
})

test_that("ml4t2e_shap_waterfall: handles max_features parameter", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("ggplot2")

  lung_df <- get_lung_survival_data()
  lung_small <- lung_df[1:30, ]

  pipeline <- ml4t2e_fit_pipeline(
    data = lung_small,
    analysis_type = "survival",
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    include_rf = FALSE
  )

  shap_result <- ml4t2e_calculate_shap(
    pipeline = pipeline,
    data = lung_small[1:10, ],
    time_horizon = 365,
    nsim = 10
  )

  # Should work with max_features
  plot <- ml4t2e_shap_waterfall(shap_result, obs_id = 1, max_features = 5)

  expect_s3_class(plot, "ggplot")
})
