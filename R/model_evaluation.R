#' Evaluate survival models on held-out data
#'
#' @description
#' Computes time-to-event accuracy diagnostics (Brier scores and concordance
#' indices) for single survival models, `SurvEnsemble` objects, or survival
#' pipelines created with [ml4t2e_fit_pipeline()]. Predictions are generated
#' internally and summarised in a leaderboard ordered by the integrated Brier
#' score.
#'
#' @param models A fitted survival model, a `SurvEnsemble`, or an
#'   `ml4t2e_pipeline` with `analysis_type = "survival"`.
#' @param data Evaluation data frame containing the time and event columns used
#'   during training.
#' @param timevar Name of the time-to-event column. Optional when `models`
#'   already stores this information (e.g., `SurvEnsemble`, pipeline).
#' @param eventvar Name of the event indicator column (0 = censored, 1 = event).
#'   Optional when stored within `models`.
#' @param eval_times Numeric vector of times at which to summarise the Brier
#'   score. Defaults to the prediction grid used for scoring.
#' @param prediction_times Numeric vector of times passed to the prediction
#'   routine. Defaults to `eval_times` when supplied, otherwise a grid spanning
#'   the evaluation data.
#' @param ensemble_method Ensemble strategy forwarded to [PredictSurvModels()]
#'   (ignored for single models). Defaults to `"average"`.
#'
#' @return A data frame sorted by `integrated_brier` containing:
#'   \itemize{
#'     \item `model` – identifier for the model.
#'     \item `integrated_brier` – integrated Brier score over `eval_times`.
#'     \item `integrated_c` – time-averaged concordance index.
#'     \item `brier_scores` – named numeric vector of Brier scores per time.
#'     \item `c_index` – named numeric vector of concordance values per time.
#'     \item `eval_times` – numeric vector of evaluation times.
#'     \item `n_observations` – number of subjects used in scoring.
#'     \item `type` – `"individual"` or `"ensemble"`.
#'     \item `rank` – rank by integrated Brier score (lower is better).
#'   }
#' @export
EvaluateSurvModels <- function(models,
                               data,
                               timevar = NULL,
                               eventvar = NULL,
                               eval_times = NULL,
                               prediction_times = NULL,
                               ensemble_method = "average") {

  if (missing(models) || is.null(models)) {
    stop("'models' must be supplied for evaluation.")
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("'data' must be a data.frame for evaluation.")
  }

  prepared <- .ml4t2e_prepare_surv_evaluation(
    models = models,
    data = data,
    timevar = timevar,
    eventvar = eventvar,
    prediction_times = prediction_times,
    ensemble_method = ensemble_method
  )

  predictions <- prepared$predictions
  pred_times <- prepared$prediction_times
  timevar <- prepared$timevar
  eventvar <- prepared$eventvar
  eval_data <- prepared$evaluation_data

  if (length(predictions) == 0) {
    stop("No valid survival predictions available for evaluation.")
  }

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }
  pred_times <- .ml4t2e_compute_time_grid(
    observed_times = if (timevar %in% colnames(eval_data)) eval_data[[timevar]] else NULL,
    override = pred_times,
    fallback = pred_times
  )
  eval_times <- .ml4t2e_compute_time_grid(
    observed_times = pred_times,
    override = eval_times,
    fallback = pred_times
  )

  obstimes <- as.numeric(eval_data[[timevar]])
  obsevents <- eval_data[[eventvar]]

  if (any(!is.finite(obstimes))) {
    valid_idx <- which(is.finite(obstimes))
    if (length(valid_idx) == 0) {
      stop("No finite observed times available for evaluation.")
    }
    obstimes <- obstimes[valid_idx]
    obsevents <- obsevents[valid_idx]
    predictions <- lapply(predictions, function(mat) mat[, valid_idx, drop = FALSE])
  }

  obsevents <- ml4t2e_validate_events(obsevents)

  result_rows <- lapply(names(predictions), function(name) {
    pred_mat <- predictions[[name]]

    brier_obj <- BrierScore(
      predsurv = pred_mat,
      pred_times = pred_times,
      obstimes = obstimes,
      obsevents = obsevents,
      eval_times = eval_times
    )
    brier_values <- as.numeric(brier_obj$AppErr$model)
    names(brier_values) <- brier_obj$time

    integrated_brier <- integratedBrier(
      predsurv = pred_mat,
      pred_times = pred_times,
      obstimes = obstimes,
      obsevents = obsevents,
      eval_times = eval_times
    )

    cindex_obj <- timedepConcordance(
      predsurv = pred_mat,
      pred_times = pred_times,
      obstimes = obstimes,
      obsevents = obsevents
    )
    cindex_values <- as.numeric(cindex_obj$AppCindex$matrix)
    names(cindex_values) <- cindex_obj$AppCindex$time

    integrated_c <- integratedC(
      predsurv = pred_mat,
      pred_times = pred_times,
      obstimes = obstimes,
      obsevents = obsevents
    )

    data.frame(
      model = name,
      integrated_brier = integrated_brier,
      integrated_c = integrated_c,
      brier_scores = I(list(brier_values)),
      c_index = I(list(cindex_values)),
      eval_times = I(list(eval_times)),
      n_observations = length(obstimes),
      type = ifelse(name == "Ensemble", "ensemble", "individual"),
      stringsAsFactors = FALSE
    )
  })

  leaderboard <- do.call(rbind, result_rows)
  leaderboard <- leaderboard[order(leaderboard$integrated_brier), , drop = FALSE]
  leaderboard$rank <- seq_len(nrow(leaderboard))
  rownames(leaderboard) <- NULL
  leaderboard
}

#' Evaluate competing risks models on held-out data
#'
#' @description
#' Computes cumulative-incidence accuracy diagnostics (Brier scores and
#' concordance indices) for single competing-risks models, `CREnsemble`
#' objects, or competing-risks pipelines.
#'
#' @inheritParams EvaluateSurvModels
#' @param cause Integer code of the event of interest (default `1`).
#'
#' @return A leaderboard data frame with the same columns described for
#'   [EvaluateSurvModels()], ranked by integrated Brier score.
#' @export
EvaluateCRModels <- function(models,
                             data,
                             timevar = NULL,
                             eventvar = NULL,
                             eval_times = NULL,
                             prediction_times = NULL,
                             ensemble_method = "average",
                             cause = 1) {

  if (missing(models) || is.null(models)) {
    stop("'models' must be supplied for evaluation.")
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("'data' must be a data.frame for evaluation.")
  }

  prepared <- .ml4t2e_prepare_cr_evaluation(
    models = models,
    data = data,
    timevar = timevar,
    eventvar = eventvar,
    prediction_times = prediction_times,
    ensemble_method = ensemble_method,
    cause = cause
  )

  predictions <- prepared$predictions
  pred_times <- prepared$prediction_times
  timevar <- prepared$timevar
  eventvar <- prepared$eventvar
  eval_data <- prepared$evaluation_data

  if (length(predictions) == 0) {
    stop("No valid competing risks predictions available for evaluation.")
  }

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }

  obstimes <- as.numeric(eval_data[[timevar]])

  pred_times <- .ml4t2e_compute_time_grid(
    observed_times = obstimes,
    override = pred_times,
    fallback = pred_times
  )
  eval_times <- .ml4t2e_compute_time_grid(
    observed_times = pred_times,
    override = eval_times,
    fallback = pred_times
  )

  obsevents <- eval_data[[eventvar]]

  valid_idx <- which(is.finite(obstimes) & !is.na(obsevents))
  if (length(valid_idx) == 0) {
    stop("No valid observations available for competing risks evaluation.")
  }
  if (length(valid_idx) < length(obstimes)) {
    obstimes <- obstimes[valid_idx]
    obsevents <- obsevents[valid_idx]
    predictions <- lapply(predictions, function(mat) mat[, valid_idx, drop = FALSE])
  }

  surv_obj <- survival::Surv(obstimes, obsevents, type = "mstate")

  result_rows <- lapply(names(predictions), function(name) {
    pred_mat <- predictions[[name]]

    brier_values <- BrierScoreCR(
      SurvObj = surv_obj,
      Predictions = pred_mat,
      eval_times = eval_times,
      cause = cause,
      pred_times = pred_times
    )
    names(brier_values) <- eval_times

    integrated_brier <- integratedBrierCR(
      SurvObj = surv_obj,
      Predictions = pred_mat,
      eval_times = eval_times,
      cause = cause,
      pred_times = pred_times
    )

    concordance_values <- sapply(eval_times, function(t) {
      timedepConcordanceCR(
        SurvObj = surv_obj,
        Predictions = pred_mat,
        time = t,
        cause = cause,
        pred_times = pred_times
      )
    })
    names(concordance_values) <- eval_times

    integrated_c <- integratedConcordanceCR(
      SurvObj = surv_obj,
      Predictions = pred_mat,
      eval_times = eval_times,
      cause = cause,
      pred_times = pred_times
    )

    data.frame(
      model = name,
      integrated_brier = integrated_brier,
      integrated_c = integrated_c,
      brier_scores = I(list(brier_values)),
      c_index = I(list(concordance_values)),
      eval_times = I(list(eval_times)),
      n_observations = length(obstimes),
      type = ifelse(name == "Ensemble", "ensemble", "individual"),
      stringsAsFactors = FALSE
    )
  })

  leaderboard <- do.call(rbind, result_rows)
  leaderboard <- leaderboard[order(leaderboard$integrated_brier), , drop = FALSE]
  leaderboard$rank <- seq_len(nrow(leaderboard))
  rownames(leaderboard) <- NULL
  leaderboard
}

.ml4t2e_prepare_surv_evaluation <- function(models,
                                            data,
                                            timevar,
                                            eventvar,
                                            prediction_times,
                                            ensemble_method) {
  if (inherits(models, "ml4t2e_pipeline")) {
    if (!identical(models$analysis_type, "survival")) {
      stop("Pipeline analysis_type must be 'survival' for EvaluateSurvModels().")
    }
    if (is.null(timevar)) timevar <- models$timevar
    if (is.null(eventvar)) eventvar <- models$eventvar
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = models$prediction_grid
    )
    pred_obj <- predict(
      models,
      newdata = data,
      new_times = pred_grid,
      ensemble_method = ensemble_method
    )
    pred_list <- pred_obj$predictions$ModelPredictions
    if (!is.null(pred_obj$predictions$NewProbs)) {
      pred_list[["Ensemble"]] <- pred_obj$predictions$NewProbs
    }
    pred_list <- Filter(Negate(is.null), pred_list)
    pred_times <- pred_obj$predictions$NewTimes
    if (is.null(pred_times)) {
      pred_times <- pred_grid
    }
    evaluation_data <- data
  } else if (inherits(models, "SurvEnsemble")) {
    if (is.null(timevar)) timevar <- models$input$timevar
    if (is.null(eventvar)) eventvar <- models$input$eventvar
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = NULL
    )
    pred_obj <- PredictSurvModels(
      models = models,
      newdata = data,
      new_times = pred_grid,
      ensemble_method = ensemble_method
    )
    pred_list <- pred_obj$ModelPredictions
    if (!is.null(pred_obj$NewProbs)) {
      pred_list[["Ensemble"]] <- pred_obj$NewProbs
    }
    pred_list <- Filter(Negate(is.null), pred_list)
    pred_times <- pred_obj$NewTimes
    if (is.null(pred_times)) {
      pred_times <- pred_grid
    }
    evaluation_data <- data
  } else {
    predictor <- .ml4t2e_lookup_surv_predictor(models)
    if (is.null(timevar)) timevar <- predictor$timevar(models)
    if (is.null(eventvar)) eventvar <- predictor$eventvar(models)
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = predictor$times(models)
    )
    prediction <- predictor$predict(models, data, pred_grid)
    pred_list <- setNames(list(prediction$Probs), prediction$label)
    pred_times <- prediction$Times
    evaluation_data <- data
  }

  if (is.null(timevar) || is.null(eventvar)) {
    stop("Unable to determine 'timevar' and 'eventvar' for survival evaluation.")
  }

  list(
    predictions = pred_list,
    prediction_times = pred_times,
    timevar = timevar,
    eventvar = eventvar,
    evaluation_data = evaluation_data
  )
}

.ml4t2e_prepare_cr_evaluation <- function(models,
                                          data,
                                          timevar,
                                          eventvar,
                                          prediction_times,
                                          ensemble_method,
                                          cause) {
  if (inherits(models, "ml4t2e_pipeline")) {
    if (!identical(models$analysis_type, "competing_risks")) {
      stop("Pipeline analysis_type must be 'competing_risks' for EvaluateCRModels().")
    }
    if (is.null(timevar)) timevar <- models$timevar
    if (is.null(eventvar)) eventvar <- models$eventvar
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = models$prediction_grid
    )
    pred_obj <- predict(
      models,
      newdata = data,
      new_times = pred_grid,
      ensemble_method = ensemble_method
    )
    pred_list <- pred_obj$predictions$ModelPredictions
    if (!is.null(pred_obj$predictions$NewProbs)) {
      pred_list[["Ensemble"]] <- pred_obj$predictions$NewProbs
    }
    pred_list <- Filter(Negate(is.null), pred_list)
    pred_times <- pred_obj$predictions$NewTimes
    if (is.null(pred_times)) {
      pred_times <- pred_grid
    }
    evaluation_data <- data
  } else if (inherits(models, "CREnsemble")) {
    if (is.null(timevar)) timevar <- models$input$timevar
    if (is.null(eventvar)) eventvar <- models$input$eventvar
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = NULL
    )
    pred_obj <- PredictCRModels(
      models = models,
      newdata = data,
      new_times = pred_grid,
      ensemble_method = ensemble_method
    )
    pred_list <- pred_obj$ModelPredictions
    if (!is.null(pred_obj$NewProbs)) {
      pred_list[["Ensemble"]] <- pred_obj$NewProbs
    }
    pred_list <- Filter(Negate(is.null), pred_list)
    pred_times <- pred_obj$NewTimes
    if (is.null(pred_times)) {
      pred_times <- pred_grid
    }
    evaluation_data <- data
  } else {
    predictor <- .ml4t2e_lookup_cr_predictor(models)
    if (is.null(timevar)) timevar <- predictor$timevar(models)
    if (is.null(eventvar)) eventvar <- predictor$eventvar(models)
    observed_times <- if (!is.null(timevar) && timevar %in% colnames(data)) data[[timevar]] else NULL
    pred_grid <- .ml4t2e_compute_time_grid(
      observed_times = observed_times,
      override = prediction_times,
      fallback = predictor$times(models)
    )
    prediction <- predictor$predict(models, data, pred_grid, cause)
    pred_list <- setNames(list(prediction$CIFs), prediction$label)
    pred_times <- prediction$Times
    evaluation_data <- data
  }

  if (is.null(timevar) || is.null(eventvar)) {
    stop("Unable to determine 'timevar' and 'eventvar' for competing risks evaluation.")
  }

  list(
    predictions = pred_list,
    prediction_times = pred_times,
    timevar = timevar,
    eventvar = eventvar,
    evaluation_data = evaluation_data
  )
}

.ml4t2e_lookup_surv_predictor <- function(model) {
  predictor_map <- list(
    ml4t2e_surv_cox = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_Cox(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_Cox"
    ),
    ml4t2e_surv_rf = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_RF(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_RF"
    ),
    ml4t2e_surv_bart = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_BART(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_BART"
    ),
    ml4t2e_surv_shallownn = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_ShallowNN(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "SurvModel_ShallowNN"
    ),
    ml4t2e_surv_gam = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_GAM(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_GAM"
    ),
    ml4t2e_surv_gbm = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_gbm(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "SurvModel_gbm"
    ),
    ml4t2e_surv_rulefit = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_rulefit(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_rulefit"
    ),
    ml4t2e_surv_xgboost = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_xgboost(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "SurvModel_xgboost"
    ),
    ml4t2e_surv_ttah = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_TTAH(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "SurvModel_TTAH"
    ),
    ml4t2e_surv_survreg = list(
      predict = function(object, newdata, new_times) {
        Predict_SurvModel_SurvReg(object, newdata, new_times)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$times)) object$times else NULL,
      label = "SurvModel_SurvReg"
    )
  )

  classes <- class(model)
  for (cls in classes) {
    if (!is.null(predictor_map[[cls]])) {
      entry <- predictor_map[[cls]]
      return(list(
        predict = function(object, newdata, new_times) {
          res <- tryCatch(
            entry$predict(object, newdata, new_times),
            error = function(e) {
              if (!is.null(new_times)) {
                entry$predict(object, newdata, NULL)
              } else {
                stop(e)
              }
            }
          )
          res$label <- entry$label
          res
        },
        timevar = function(object) entry$timevar(object),
        eventvar = function(object) entry$eventvar(object),
        times = function(object) entry$times(object)
      ))
    }
  }

  stop("Unsupported survival model class: ", paste(classes, collapse = ", "))
}

.ml4t2e_lookup_cr_predictor <- function(model) {
  predictor_map <- list(
    ml4t2e_cr_finegray = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_FineGray(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_FineGray"
    ),
    ml4t2e_cr_cox = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_Cox(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_Cox"
    ),
    ml4t2e_cr_rf = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_RF(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_RF"
    ),
    ml4t2e_cr_bart = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_BART(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_BART"
    ),
    ml4t2e_cr_rulefit = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_rulefit(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_rulefit"
    ),
    ml4t2e_cr_gam = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_GAM(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_GAM"
    ),
    ml4t2e_cr_survreg = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_SurvReg(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_SurvReg"
    ),
    ml4t2e_cr_xgboost = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_xgboost(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_xgboost"
    ),
    ml4t2e_cr_shallownn = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_ShallowNN(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_ShallowNN"
    ),
    ml4t2e_cr_ttah = list(
      predict = function(object, newdata, new_times, cause) {
        Predict_CRModel_TTAH(object, newdata, new_times, event_of_interest = cause)
      },
      timevar = function(object) if (!is.null(object$timevar)) object$timevar else NULL,
      eventvar = function(object) if (!is.null(object$eventvar)) object$eventvar else NULL,
      times = function(object) if (!is.null(object$Times)) object$Times else NULL,
      label = "CRModel_TTAH"
    )
  )

  classes <- class(model)
  for (cls in classes) {
    if (!is.null(predictor_map[[cls]])) {
      entry <- predictor_map[[cls]]
      return(list(
        predict = function(object, newdata, new_times, cause) {
          res <- tryCatch(
            entry$predict(object, newdata, new_times, cause),
            error = function(e) {
              if (!is.null(new_times)) {
                entry$predict(object, newdata, NULL, cause)
              } else {
                stop(e)
              }
            }
          )
          cif_matrix <- res$CIFs
          times <- res$Times
          if (!is.null(new_times)) {
            cif_matrix <- cifMatInterpolaltor(
              probsMat = cif_matrix,
              times = times,
              new_times = new_times
            )
            times <- new_times
          }
          list(
            CIFs = cif_matrix,
            Times = times,
            label = entry$label
          )
        },
        timevar = function(object) entry$timevar(object),
        eventvar = function(object) entry$eventvar(object),
        times = function(object) entry$times(object)
      ))
    }
  }

  stop("Unsupported competing risks model class: ", paste(classes, collapse = ", "))
}
