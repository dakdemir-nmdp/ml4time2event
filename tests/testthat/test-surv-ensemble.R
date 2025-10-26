library(testthat)
library(survival)

set.seed(789)

n_obs_surv <- 60

surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)),
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)),
  x3 = rnorm(n_obs_surv, mean = 2),
  stringsAsFactors = FALSE
)

train_indices_surv <- seq_len(45)
test_indices_surv <- setdiff(seq_len(n_obs_surv), train_indices_surv)
train_data_surv <- surv_data[train_indices_surv, ]
test_data_surv <- surv_data[test_indices_surv, ]
event_times <- train_data_surv$time[train_data_surv$status == 1]
if (length(event_times) < 3L) {
  event_times <- quantile(train_data_surv$time, probs = c(0.25, 0.5, 0.75))
} else {
  event_times <- quantile(event_times, probs = c(0.25, 0.5, 0.75))
}

test_that("RunSurvModels returns SurvEnsemble with expected components", {
  skip_if_not_installed("randomForestSRC")

  fitted <- RunSurvModels(
    datatrain = train_data_surv,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 50
  )

  expect_s3_class(fitted, "SurvEnsemble")
  expect_true(all(c("input", "model_status", "RF_Model", "RF_Model2") %in% names(fitted)))
  expect_named(fitted$model_status)
  expect_true("RF_Model" %in% names(fitted$model_status))
  expect_true("RF_Model2" %in% names(fitted$model_status))
  expect_true("CPH_Model" %in% names(fitted$model_status))
})

test_that("PredictSurvModels excludes requested models that failed to train", {
  skip_if_not_installed("randomForestSRC")

  fitted <- RunSurvModels(
    datatrain = train_data_surv,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 40
  )

  fitted$CPH_Model <- NULL
  fitted$model_status["CPH_Model"] <- FALSE

  expect_warning(
    preds <- PredictSurvModels(
      models = fitted,
      newdata = test_data_surv,
      new_times = event_times,
      models_to_use = c("RF_Model", "CPH_Model")
    ),
    "failed during training"
  )

  expect_type(preds, "list")
  expect_true("models_used" %in% names(preds))
  expect_true("RF_Model" %in% preds$models_used)
  expect_false("CPH_Model" %in% preds$models_used)
})

test_that("PredictSurvModels supports weighted ensembles", {
  skip_if_not_installed("randomForestSRC")

  fitted <- RunSurvModels(
    datatrain = train_data_surv,
    ExpVars = c("x1", "x2", "x3"),
    timevar = "time",
    eventvar = "status",
    models = c("coxph"),
    ntreeRF = 50
  )

  weights <- c(RF_Model = 0.25, RF_Model2 = 0.25, CPH_Model = 0.5)

  preds <- PredictSurvModels(
    models = fitted,
    newdata = test_data_surv,
    new_times = event_times,
    models_to_use = names(weights),
    ensemble_method = "weighted",
    model_weights = weights
  )

  expect_equal(preds$ensemble_method, "weighted")
  expect_true(all(names(preds$ModelPredictions) %in% names(weights)))
  expect_true(is.matrix(preds$NewProbs))
  expect_equal(ncol(preds$NewProbs), nrow(test_data_surv))
  expect_equal(nrow(preds$NewProbs), length(unique(c(0, event_times))))
  expect_true(all(preds$NewProbs >= 0 & preds$NewProbs <= 1))
})
