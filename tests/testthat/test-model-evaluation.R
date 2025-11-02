library(testthat)
library(survival)

# Ensure package functions are available for testing
# This handles both installed package and development modes
if (!exists("SurvModel_Cox", mode = "function")) {
  # Try loading via devtools (development mode)
  if (requireNamespace("devtools", quietly = TRUE)) {
    # Find package root (tests/testthat -> ../.. -> package root)
    test_path <- getwd()
    if (basename(test_path) == "testthat") {
      pkg_root <- normalizePath(file.path(test_path, "..", ".."))
    } else if (file.exists(file.path(test_path, "tests", "testthat"))) {
      pkg_root <- test_path
    } else {
      pkg_root <- normalizePath(file.path(test_path, ".."))
    }
    devtools::load_all(pkg_root, quiet = TRUE, export_all = FALSE)
  } else {
    # Fall back to library() if devtools not available
    library(ml4time2event)
  }
} else {
  # Ensure ml4time2event is loaded if functions already exist
  if (!"package:ml4time2event" %in% search()) {
    library(ml4time2event)
  }
}

test_that("EvaluateSurvModels returns leaderboard for ensembles", {
  data("lung_survival_prepared", package = "ml4time2event")
  set.seed(42)
  n <- nrow(lung_survival_prepared)
  idx <- sample(seq_len(n), size = floor(0.7 * n))
  train_data <- lung_survival_prepared[idx, , drop = FALSE]
  test_data <- lung_survival_prepared[-idx, , drop = FALSE]

  expvars <- c("age", "sex", "ph.ecog", "ph.karno", "meal.cal", "wt.loss")
  timevar <- "time"
  eventvar <- "status"

  cox_model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    varsel = "none",
    ntimes = 20
  )
  rf_model <- SurvModel_RF(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    ntree = 50,
    samplesize = min(100, ceiling(0.4 * nrow(train_data)))
  )

  surv_ensemble <- list(
    input = list(timevar = timevar, eventvar = eventvar),
    model_status = c(CPH_Model = TRUE, RF_Model = TRUE),
    CPH_Model = cox_model,
    RF_Model = rf_model
  )
  class(surv_ensemble) <- c("SurvEnsemble", "list")

  eval_times <- seq(0, max(test_data[[timevar]]), length.out = 15)
  eval_times <- unique(sort(c(0, eval_times)))

  leaderboard <- EvaluateSurvModels(
    models = surv_ensemble,
    data = test_data,
    timevar = timevar,
    eventvar = eventvar,
    eval_times = eval_times,
    prediction_times = eval_times
  )

  expect_s3_class(leaderboard, "data.frame")
  expect_true(all(c("model", "integrated_brier", "integrated_c") %in% names(leaderboard)))
  expect_true("Ensemble" %in% leaderboard$model)
  expect_equal(nrow(leaderboard), 3)

  surv_preds <- PredictSurvModels(
    models = surv_ensemble,
    newdata = test_data,
    new_times = eval_times
  )
  direct_ibrier <- integratedBrier(
    predsurv = surv_preds$ModelPredictions$CPH_Model,
    pred_times = surv_preds$NewTimes,
    obstimes = test_data[[timevar]],
    obsevents = test_data[[eventvar]]
  )
  cox_row <- leaderboard[leaderboard$model == "CPH_Model", , drop = FALSE]
  expect_equal(cox_row$integrated_brier, direct_ibrier, tolerance = 1e-6)
})

test_that("EvaluateSurvModels supports single model objects", {
  data("lung_survival_prepared", package = "ml4time2event")
  set.seed(123)
  n <- nrow(lung_survival_prepared)
  idx <- sample(seq_len(n), size = floor(0.7 * n))
  train_data <- lung_survival_prepared[idx, , drop = FALSE]
  test_data <- lung_survival_prepared[-idx, , drop = FALSE]

  expvars <- c("age", "sex", "ph.ecog", "ph.karno", "meal.cal", "wt.loss")
  timevar <- "time"
  eventvar <- "status"

  cox_model <- SurvModel_Cox(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    varsel = "none",
    ntimes = 15
  )

  eval_times <- seq(0, max(test_data[[timevar]]), length.out = 12)
  eval_times <- unique(sort(c(0, eval_times)))

  leaderboard <- EvaluateSurvModels(
    models = cox_model,
    data = test_data,
    timevar = timevar,
    eventvar = eventvar,
    eval_times = eval_times,
    prediction_times = eval_times
  )

  expect_equal(nrow(leaderboard), 1)
  expect_true(all(c("model", "integrated_brier", "integrated_c") %in% names(leaderboard)))

  direct_preds <- Predict_SurvModel_Cox(cox_model, test_data, new_times = eval_times)
  direct_ibrier <- integratedBrier(
    predsurv = direct_preds$Probs,
    pred_times = direct_preds$Times,
    obstimes = test_data[[timevar]],
    obsevents = test_data[[eventvar]]
  )
  expect_equal(leaderboard$integrated_brier[[1]], direct_ibrier, tolerance = 1e-6)
})

test_that("EvaluateCRModels returns leaderboard for ensembles", {
  data("bmt_competing_risks_prepared", package = "ml4time2event")
  set.seed(99)
  n <- nrow(bmt_competing_risks_prepared)
  idx <- sample(seq_len(n), size = floor(0.7 * n))
  train_data <- bmt_competing_risks_prepared[idx, , drop = FALSE]
  test_data <- bmt_competing_risks_prepared[-idx, , drop = FALSE]

  expvars <- c("sex", "d", "phase", "age", "source")
  timevar <- "ftime"
  eventvar <- "status"

  fg_model <- CRModel_FineGray(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = "1"
  )
  cox_model <- CRModel_Cox(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar
  )

  cr_ensemble <- list(
    input = list(timevar = timevar, eventvar = eventvar),
    model_status = c(FG_Model = TRUE, Cox_Model = TRUE),
    FG_Model = fg_model,
    Cox_Model = cox_model
  )
  class(cr_ensemble) <- c("CREnsemble", "list")

  eval_times <- seq(0, max(test_data[[timevar]]), length.out = 12)
  eval_times <- unique(sort(c(0, eval_times)))

  leaderboard <- EvaluateCRModels(
    models = cr_ensemble,
    data = test_data,
    timevar = timevar,
    eventvar = eventvar,
    eval_times = eval_times,
    prediction_times = eval_times,
    cause = 1
  )

  expect_s3_class(leaderboard, "data.frame")
  expect_true(all(c("model", "integrated_brier", "integrated_c") %in% names(leaderboard)))
  expect_true("Ensemble" %in% leaderboard$model)
  expect_equal(nrow(leaderboard), 3)

  cr_preds <- PredictCRModels(
    models = cr_ensemble,
    newdata = test_data,
    new_times = eval_times
  )
  surv_obj <- Surv(test_data[[timevar]], test_data[[eventvar]], type = "mstate")
  direct_ibrier <- integratedBrierCR(
    SurvObj = surv_obj,
    Predictions = cr_preds$ModelPredictions$Cox_Model,
    eval_times = eval_times,
    cause = 1,
    pred_times = eval_times
  )
  cox_row <- leaderboard[leaderboard$model == "Cox_Model", , drop = FALSE]
  expect_equal(cox_row$integrated_brier, direct_ibrier, tolerance = 1e-6)
})

test_that("EvaluateCRModels supports single model objects", {
  data("bmt_competing_risks_prepared", package = "ml4time2event")
  set.seed(321)
  n <- nrow(bmt_competing_risks_prepared)
  idx <- sample(seq_len(n), size = floor(0.7 * n))
  train_data <- bmt_competing_risks_prepared[idx, , drop = FALSE]
  test_data <- bmt_competing_risks_prepared[-idx, , drop = FALSE]

  expvars <- c("sex", "d", "phase", "age", "source")
  timevar <- "ftime"
  eventvar <- "status"

  fg_model <- CRModel_FineGray(
    data = train_data,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = "1"
  )

  eval_times <- seq(0, max(test_data[[timevar]]), length.out = 10)
  eval_times <- unique(sort(c(0, eval_times)))

  leaderboard <- EvaluateCRModels(
    models = fg_model,
    data = test_data,
    timevar = timevar,
    eventvar = eventvar,
    eval_times = eval_times,
    prediction_times = eval_times,
    cause = 1
  )

  expect_equal(nrow(leaderboard), 1)
  expect_true(all(c("model", "integrated_brier", "integrated_c") %in% names(leaderboard)))

  fg_preds <- Predict_CRModel_FineGray(fg_model, newdata = test_data, new_times = eval_times)
  surv_obj <- Surv(test_data[[timevar]], test_data[[eventvar]], type = "mstate")
  direct_ibrier <- integratedBrierCR(
    SurvObj = surv_obj,
    Predictions = fg_preds$CIFs,
    eval_times = eval_times,
    cause = 1,
    pred_times = eval_times
  )
  expect_equal(leaderboard$integrated_brier[[1]], direct_ibrier, tolerance = 1e-6)
})
