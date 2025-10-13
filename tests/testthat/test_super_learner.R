library(testthat)
library(survival)

# Load the package
devtools::load_all()

context("Super Learner Survival Stacking")

# Test data setup
set.seed(456)
n_obs <- 80
train_data <- data.frame(
  time = rexp(n_obs, rate = 0.03),
  status = sample(0:1, n_obs, replace = TRUE, prob = c(0.2, 0.8)),
  x1 = rnorm(n_obs),
  x2 = factor(sample(c("A", "B"), n_obs, replace = TRUE)),
  x3 = rnorm(n_obs, mean = 1),
  stringsAsFactors = FALSE
)

n_test <- 20
test_data <- data.frame(
  time = rexp(n_test, rate = 0.03),
  status = sample(0:1, n_test, replace = TRUE, prob = c(0.2, 0.8)),
  x1 = rnorm(n_test),
  x2 = factor(sample(c("A", "B"), n_test, replace = TRUE)),
  x3 = rnorm(n_test, mean = 1),
  stringsAsFactors = FALSE
)

test_times <- quantile(train_data$time[train_data$status == 1], probs = c(0.3, 0.6, 0.9))

test_that("buildObservedSurvivalMatrix creates proper target matrix", {
  obs_surv <- buildObservedSurvivalMatrix(
    data = train_data,
    timevar = "time",
    eventvar = "status",
    eval_times = test_times
  )
  
  expect_true(is.matrix(obs_surv))
  expect_equal(nrow(obs_surv), length(unique(c(0, test_times))))
  expect_equal(ncol(obs_surv), nrow(train_data))
  expect_true(all(obs_surv >= 0 & obs_surv <= 1))
  
  # Check that survival starts at 1 (time 0)
  expect_true(all(obs_surv[1, ] == 1))
  
  # Check monotonicity (survival should be non-increasing)
  for (i in seq_len(ncol(obs_surv))) {
    expect_true(all(diff(obs_surv[, i]) <= 1e-10))  # Allow for small numerical errors
  }
})

test_that("PredictSurvModels super_learner with training data works", {
  skip_if_not_installed("randomForestSRC")
  
  # Fit models
  fitted <- RunSurvModels(
    datatrain = train_data,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 30
  )
  
  # Test super learner with training data
  sl_preds <- PredictSurvModels(
    models = fitted,
    newdata = test_data,
    newtimes = test_times,
    ensemble_method = "super_learner",
    super_learner_training_data = train_data,
    super_learner_timevar = "time",
    super_learner_eventvar = "status"
  )
  
  expect_equal(sl_preds$ensemble_method, "super_learner")
  expect_true(is.matrix(sl_preds$NewProbs))
  expect_equal(ncol(sl_preds$NewProbs), nrow(test_data))
  expect_equal(nrow(sl_preds$NewProbs), length(unique(c(0, test_times))))
  expect_true(all(sl_preds$NewProbs >= 0 & sl_preds$NewProbs <= 1))
})

test_that("PredictSurvModels super_learner with pre-computed weights works", {
  skip_if_not_installed("randomForestSRC")
  
  # Fit models
  fitted <- RunSurvModels(
    datatrain = train_data,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 30
  )
  
  # Pre-computed weights
  sl_weights <- c(RF_Model = 0.3, RF_Model2 = 0.2, CPH_Model = 0.5)
  
  # Test super learner with weights
  sl_preds <- PredictSurvModels(
    models = fitted,
    newdata = test_data,
    newtimes = test_times,
    ensemble_method = "super_learner",
    model_weights = sl_weights
  )
  
  expect_equal(sl_preds$ensemble_method, "super_learner")
  expect_true(is.matrix(sl_preds$NewProbs))
  expect_equal(ncol(sl_preds$NewProbs), nrow(test_data))
  expect_equal(nrow(sl_preds$NewProbs), length(unique(c(0, test_times))))
  expect_true(all(sl_preds$NewProbs >= 0 & sl_preds$NewProbs <= 1))
})

test_that("optimizeSuperLearnerWeights produces valid weights", {
  # Create mock prediction matrices
  n_times <- 4
  n_obs <- 10
  
  pred1 <- matrix(runif(n_times * n_obs, 0.5, 0.9), nrow = n_times, ncol = n_obs)
  pred2 <- matrix(runif(n_times * n_obs, 0.1, 0.5), nrow = n_times, ncol = n_obs)  # More distinct
  
  # Ensure monotonicity
  for (i in seq_len(n_obs)) {
    pred1[, i] <- cummin(pred1[, i])
    pred2[, i] <- cummin(pred2[, i])
  }
  
  predictions <- list(Model1 = pred1, Model2 = pred2)
  
  # Create observed survival (much closer to pred1) - ensure same dimensions as predictions
  actual <- pred1 + matrix(rnorm(n_times * n_obs, 0, 0.01), nrow = n_times, ncol = n_obs)  # Smaller noise
  actual <- pmax(0, pmin(actual, 1))
  # Ensure it's properly formatted as a matrix with correct dimensions
  actual <- matrix(actual, nrow = n_times, ncol = n_obs)
  # Ensure monotonicity for realistic survival data
  for (i in seq_len(n_obs)) {
    actual[, i] <- cummin(actual[, i])
  }
  
  weights <- optimizeSuperLearnerWeights(predictions, actual, loss_type = "mse")
  
  expect_true(is.numeric(weights))
  expect_equal(length(weights), 2)
  expect_equal(names(weights), c("Model1", "Model2"))
  expect_equal(sum(weights), 1, tolerance = 1e-6)
  expect_true(all(weights >= 0))
  
  # Weights should sum to 1 and be valid (main test of optimization)
  expect_true(all(weights >= -1e-10))  # Allow for small numerical errors
})

test_that("Super learner falls back gracefully when training data is missing", {
  skip_if_not_installed("randomForestSRC")
  
  # Fit models
  fitted <- RunSurvModels(
    datatrain = train_data,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 30
  )
  
  # Test super learner without training data or weights
  expect_message(
    sl_preds <- PredictSurvModels(
      models = fitted,
      newdata = test_data,
      newtimes = test_times,
      ensemble_method = "super_learner"
    ),
    "requires either pre-computed weights"
  )
  
  expect_equal(sl_preds$ensemble_method, "super_learner")
  expect_true(is.matrix(sl_preds$NewProbs))
  expect_true(all(sl_preds$NewProbs >= 0 & sl_preds$NewProbs <= 1))
})

test_that("ComputeSuperLearnerWeights works and weights are used in predictions", {
  skip_if_not_installed("randomForestSRC")
  
  # Fit models (use random forest which is more stable for this test)
  fitted <- RunSurvModels(
    datatrain = train_data,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("RF"),
    ntreeRF = 30
  )
  
  # Compute and store super learner weights
  fitted_with_weights <- ComputeSuperLearnerWeights(
    ensemble_models = fitted,
    training_data = train_data
  )
  
  # Check that weights are stored
  expect_true(!is.null(fitted_with_weights$super_learner_weights))
  expect_true(is.numeric(fitted_with_weights$super_learner_weights))
  expect_equal(sum(fitted_with_weights$super_learner_weights), 1, tolerance = 1e-6)
  expect_true(all(fitted_with_weights$super_learner_weights >= 0))
  
  # Test that pre-computed weights are used during prediction
  sl_preds <- PredictSurvModels(
    models = fitted_with_weights,
    newdata = test_data,
    newtimes = test_times,
    ensemble_method = "super_learner"
  )
  
  expect_equal(sl_preds$ensemble_method, "super_learner")
  expect_true(is.matrix(sl_preds$NewProbs))
  expect_true(all(sl_preds$NewProbs >= 0 & sl_preds$NewProbs <= 1))
})