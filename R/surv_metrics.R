ml4t2e_align_surv_predictions <- function(predsurv, pred_times, obstimes, context = "predsurv") {
  if (is.null(predsurv)) {
    stop(sprintf("'%s' cannot be NULL.", context))
  }

  pred_matrix <- as.matrix(predsurv)
  if (!is.numeric(pred_matrix)) {
    stop(sprintf("'%s' must be numeric.", context))
  }
  if (any(!is.finite(pred_matrix))) {
    stop(sprintf("'%s' contains non-finite values.", context))
  }

  if (missing(pred_times) || is.null(pred_times)) {
    stop("'pred_times' must be supplied and cannot be NULL.")
  }
  if (!is.numeric(pred_times)) {
    stop("'pred_times' must be numeric.")
  }
  if (length(pred_times) == 0) {
    stop("'pred_times' cannot be empty.")
  }

  obstimes <- as.numeric(obstimes)
  if (!all(is.finite(obstimes))) {
    stop("'obstimes' must be numeric and finite.")
  }
  n_obs <- length(obstimes)

  n_times <- length(pred_times)
  dims <- dim(pred_matrix)

  orientation_fixed <- FALSE
  if (nrow(pred_matrix) == n_times && ncol(pred_matrix) == n_obs) {
    orientation_fixed <- TRUE
  } else if (ncol(pred_matrix) == n_times && nrow(pred_matrix) == n_obs) {
    pred_matrix <- t(pred_matrix)
    orientation_fixed <- TRUE
  }

  if (!orientation_fixed) {
    stop(sprintf(
      "Unable to align predictions. Expected matrix with %d time points and %d observations, got %dx%d.",
      n_times, n_obs, dims[1], dims[2]
    ))
  }

  # Ensure times are strictly increasing for interpolation
  if (any(is.na(pred_times))) {
    stop("'pred_times' contains NA values.")
  }
  order_idx <- order(pred_times)
  pred_times_sorted <- pred_times[order_idx]
  if (any(diff(pred_times_sorted) == 0)) {
    keep_idx <- !duplicated(pred_times_sorted)
    pred_times_sorted <- pred_times_sorted[keep_idx]
    pred_matrix <- pred_matrix[order_idx, , drop = FALSE][keep_idx, , drop = FALSE]
  } else {
    pred_matrix <- pred_matrix[order_idx, , drop = FALSE]
  }

  rownames(pred_matrix) <- NULL

  list(
    matrix = pred_matrix,
    times = pred_times_sorted,
    obstimes = obstimes,
    n_obs = n_obs
  )
}

ml4t2e_interpolate_survival <- function(pred_matrix, pred_times, eval_time) {
  if (length(eval_time) != 1L || !is.finite(eval_time)) {
    stop("'eval_time' must be a single finite numeric value.")
  }

  if (eval_time <= pred_times[1]) {
    return(pred_matrix[1, ])
  }
  if (eval_time >= pred_times[length(pred_times)]) {
    return(pred_matrix[nrow(pred_matrix), ])
  }

  apply(pred_matrix, 2, function(curve) {
    if (all(is.na(curve))) {
      return(NA_real_)
    }
    stats::approx(
      x = pred_times,
      y = curve,
      xout = eval_time,
      method = "linear",
      rule = 2
    )$y
  })
}

ml4t2e_validate_events <- function(obsevents) {
  if (is.null(obsevents)) {
    stop("'obsevents' cannot be NULL.")
  }
  if (!is.numeric(obsevents) && !is.logical(obsevents)) {
    stop("'obsevents' must be numeric or logical.")
  }
  ev <- as.numeric(obsevents)
  ev[is.na(ev)] <- NA_real_
  unique_vals <- unique(ev[!is.na(ev)])
  if (!all(unique_vals %in% c(0, 1))) {
    warning("'obsevents' contains values other than 0/1. Treating non-1 values as censored.")
    ev[!is.na(ev) & ev != 1] <- 0
  }
  ev
}

ml4t2e_cindex_at_time <- function(pred_surv, eval_time, obstimes, obsevents) {
  events_idx <- which(obsevents == 1 & obstimes <= eval_time & !is.na(obstimes))
  if (length(events_idx) == 0) {
    return(NA_real_)
  }

  concordant <- 0
  comparable <- 0

  risk_scores <- 1 - pred_surv

  for (i in events_idx) {
    risk_i <- risk_scores[i]
    if (!is.finite(risk_i)) next

    for (j in seq_along(obstimes)) {
      if (i == j) next
      if (!is.finite(risk_scores[j])) next
      if (obstimes[j] <= obstimes[i]) next

      comparable <- comparable + 1
      if (risk_i > risk_scores[j]) {
        concordant <- concordant + 1
      } else if (risk_i == risk_scores[j]) {
        concordant <- concordant + 0.5
      }
    }
  }

  if (comparable == 0) {
    return(NA_real_)
  }
  concordant / comparable
}

#' @title timedepConcordance
#'
#' @description Calculate time-dependent concordance for survival predictions without relying on external pec methods.
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations). If the matrix is observations x times it will be transposed automatically.
#' @param pred_times Numeric vector of prediction times matching the rows or columns of `predsurv`.
#' @param obstimes Observed follow-up times.
#' @param obsevents Observed event indicator (0=censored, 1=event).
#' @param TestMat Optional test dataset (ignored, preserved for backward compatibility).
#'
#' @return A list with element `AppCindex$matrix` containing the time-specific concordance values.
#' @export
timedepConcordance <- function(predsurv, pred_times, obstimes, obsevents, TestMat = NULL) {
  alignment <- ml4t2e_align_surv_predictions(predsurv, pred_times, obstimes, context = "predsurv")
  obsevents_numeric <- ml4t2e_validate_events(obsevents)

  c_values <- vapply(
    seq_along(alignment$times),
    function(idx) {
      surv_at_time <- alignment$matrix[idx, ]
      ml4t2e_cindex_at_time(
        pred_surv = surv_at_time,
        eval_time = alignment$times[idx],
        obstimes = alignment$obstimes,
        obsevents = obsevents_numeric
      )
    },
    numeric(1)
  )

  result <- list(
    AppCindex = list(
      matrix = c_values,
      time = alignment$times
    ),
    time = alignment$times,
    call = match.call()
  )
  class(result) <- c("ml4time2event_cindex", "pecCindex")
  result
}


#' @title BrierScore
#'
#' @description Calculate Brier score for survival predictions at specific times.
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations). If needed, orientation is fixed automatically.
#' @param pred_times Numeric vector of prediction times matching the rows or columns of `predsurv`.
#' @param obstimes Observed follow-up times.
#' @param obsevents Observed event indicator (0=censored, 1=event).
#' @param eval_times Optional numeric vector of evaluation times (defaults to `pred_times`).
#' @param TestMat Optional test dataset (ignored, preserved for backward compatibility).
#'
#' @return A `pec`-like list with element `AppErr$model` storing the Brier scores.
BrierScore <- function(predsurv, pred_times, obstimes, obsevents,
                       eval_times = NULL, TestMat = NULL) {

  alignment <- ml4t2e_align_surv_predictions(predsurv, pred_times, obstimes, context = "predsurv")
  if (is.null(eval_times)) {
    eval_times <- alignment$times
  } else {
    if (!is.numeric(eval_times)) {
      stop("'eval_times' must be numeric.")
    }
    if (any(!is.finite(eval_times))) {
      stop("'eval_times' must contain only finite values.")
    }
  }

  obsevents_numeric <- ml4t2e_validate_events(obsevents)

  brier_values <- vapply(
    eval_times,
    function(t_eval) {
      pred_surv <- ml4t2e_interpolate_survival(alignment$matrix, alignment$times, t_eval)
      valid <- is.finite(pred_surv) & is.finite(obstimes) & !is.na(obsevents_numeric)
      if (!any(valid)) {
        return(NA_real_)
      }

      event_by_t <- obsevents_numeric[valid] == 1 & obstimes[valid] <= t_eval
      surv_pred <- pred_surv[valid]

      scores <- ifelse(event_by_t, surv_pred^2, (1 - surv_pred)^2)
      if (all(is.na(scores))) {
        return(NA_real_)
      }
      mean(scores, na.rm = TRUE)
    },
    numeric(1)
  )

  result <- list(
    AppErr = list(
      model = brier_values,
      time = eval_times
    ),
    time = eval_times,
    call = match.call()
  )
  class(result) <- c("ml4time2event_brier", "pec")
  result
}


#' @title integratedBrier
#'
#' @description Calculate integrated Brier score over time range
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations)
#' @param pred_times The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param eval_times Optional vector of evaluation times for integration
#' @param TestMat Optional test dataset
#'
#' @return Integrated Brier score (scalar)
#' @export
integratedBrier <- function(predsurv, pred_times, obstimes, obsevents,
                            eval_times = NULL, TestMat = NULL) {

  # Get Brier scores at all time points
  brier_obj <- BrierScore(
    predsurv = predsurv,
    pred_times = pred_times,
    obstimes = obstimes,
    obsevents = obsevents,
    eval_times = eval_times,
    TestMat = TestMat
  )

  # Extract the Brier score values
  # pec objects have AppErr slot with the prediction error
  if (is.null(brier_obj$AppErr)) {
    stop("BrierScore output is missing the 'AppErr' component.")
  }
  brier_values <- brier_obj$AppErr$model
  eval_times <- brier_obj$time

  # Integrate using trapezoidal rule
  if (length(eval_times) < 2) {
    return(brier_values[1])
  }

  # Calculate integrated Brier score (area under curve)
  time_diffs <- diff(eval_times)
  avg_brier <- (brier_values[-1] + brier_values[-length(brier_values)]) / 2
  integrated_bs <- sum(time_diffs * avg_brier) / (max(eval_times) - min(eval_times))

  return(integrated_bs)
}


#' @title integratedC
#'
#' @description Calculate integrated concordance index over time range
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations)
#' @param pred_times The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param TestMat Optional test dataset
#'
#' @return Integrated concordance index (scalar)
#' @export
integratedC <- function(predsurv, pred_times, obstimes, obsevents, TestMat = NULL) {

  # Get time-dependent concordance
  cindex_obj <- timedepConcordance(
    predsurv = predsurv,
    pred_times = pred_times,
    obstimes = obstimes,
    obsevents = obsevents,
    TestMat = TestMat
  )

  # Extract concordance values at each time point
  if (is.null(cindex_obj$AppCindex)) {
    stop("timedepConcordance output is missing the 'AppCindex' component.")
  }
  cindex_values <- cindex_obj$AppCindex$matrix

  # Calculate mean concordance over time (integrated C-index)
  integrated_c <- mean(cindex_values, na.rm = TRUE)

  return(integrated_c)
}
