ml4t2e_align_surv_predictions <- function(predsurv, pred_times, obstimes, context = "predsurv") {
  if (is.null(predsurv)) {
    rlang::abort(sprintf("'%s' cannot be NULL.", context))
  }

  pred_matrix <- as.matrix(predsurv)
  if (!is.numeric(pred_matrix)) {
    rlang::abort(sprintf("'%s' must be numeric.", context))
  }
  if (any(!is.finite(pred_matrix))) {
    rlang::abort(sprintf("'%s' contains non-finite values.", context))
  }

  if (missing(pred_times) || is.null(pred_times)) {
    rlang::abort("'pred_times' must be supplied and cannot be NULL.")
  }
  if (!is.numeric(pred_times)) {
    rlang::abort("'pred_times' must be numeric.")
  }
  if (length(pred_times) == 0) {
    rlang::abort("'pred_times' cannot be empty.")
  }

  obstimes <- as.numeric(obstimes)
  if (!all(is.finite(obstimes))) {
    rlang::abort("'obstimes' must be numeric and finite.")
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
    rlang::abort(sprintf(
      "Unable to align predictions. Expected matrix with %d time points and %d observations, got %dx%d.",
      n_times, n_obs, dims[1], dims[2]
    ))
  }

  # Ensure times are strictly increasing for interpolation
  if (any(is.na(pred_times))) {
    rlang::abort("'pred_times' contains NA values.")
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
    rlang::abort("'eval_time' must be a single finite numeric value.")
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
    rlang::abort("'obsevents' cannot be NULL.")
  }
  if (!is.numeric(obsevents) && !is.logical(obsevents)) {
    rlang::abort("'obsevents' must be numeric or logical.")
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

#' Time-Dependent Concordance Index for Survival Analysis
#'
#' Computes Harrell's C-index at each time point in the prediction grid.
#' This assesses how well predicted survival probabilities discriminate between
#' pairs of observations where one had an earlier event.
#'
#' **Formula**:
#' At time \eqn{t}, \eqn{C(t) = \frac{\sum_{i,j} I(t_i < t_j) \cdot I(\hat{S}_i(t) > \hat{S}_j(t)) \cdot \delta_j}{\sum_{i,j} I(t_i < t_j) \cdot \delta_j}}
#'
#' **Interpretation**:
#' - Range: 0.5 (random predictions) to 1.0 (perfect discrimination)
#' - 0.5 = worse than random; 1.0 = perfect ranking
#' - Values > 0.7 typically considered good discrimination
#'
#' @param predsurv Matrix of predicted survival probabilities (rows = times, cols = observations)
#' @param pred_times Numeric vector of time points at which predictions were made
#' @param obstimes Observed event/censoring times
#' @param obsevents Binary event indicator (1 = event, 0 = censored)
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return List of class `ml4time2event_cindex` containing:
#'   - `AppCindex$matrix`: Vector of C-index values at each time point
#'   - `AppCindex$time`: Corresponding time points
#'   - `time`: Time grid
#'   - `call`: Original function call
#'
#' @references
#' Harrell, F. E., Lee, K. L., & Mark, D. B. (1996).
#' "Multivariable prognostic models." *Statistics in Medicine*, 15(4), 361–387.
#'
#' @keywords internal
#'
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


#' Brier Score for Survival Predictions
#'
#' Computes mean squared error between predicted survival probabilities and
#' observed binary event status at specified time points. Lower values indicate
#' better calibration.
#'
#' **Formula**:
#' \eqn{BS(t) = E[(S(t|X) - Y(t))^2]}
#'
#' where \eqn{S(t|X)} is predicted survival probability and \eqn{Y(t)} is binary
#' event status at time \eqn{t}.
#'
#' **Interpretation**:
#' - Range: 0 (perfect calibration) to 1 (worst)
#' - Combines discrimination and calibration
#' - 0.25 at all times = baseline for completely random predictions
#'
#' @param predsurv Matrix of predicted survival probabilities (rows = times, cols = observations)
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param eval_times Optional evaluation time points. If NULL, computed at all `pred_times`.
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return List of class `ml4time2event_brier` containing Brier scores and time grid
#'
#' @references
#' Brier, G. W. (1950). "Verification of forecasts expressed in terms of
#' probability." *Monthly Weather Review*, 78(1), 1–3.
#'
#' @keywords internal
#'
BrierScore <- function(predsurv, pred_times, obstimes, obsevents,
                       eval_times = NULL, TestMat = NULL) {

  alignment <- ml4t2e_align_surv_predictions(predsurv, pred_times, obstimes, context = "predsurv")
  if (is.null(eval_times)) {
    eval_times <- alignment$times
  } else {
    if (!is.numeric(eval_times)) {
      rlang::abort("'eval_times' must be numeric.")
    }
    if (any(!is.finite(eval_times))) {
      rlang::abort("'eval_times' must contain only finite values.")
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

      # Only include observations at risk at time t_eval
      # At risk if: obstimes > t_eval (still at risk) OR (obstimes <= t_eval AND obsevents == 1) (event occurred)
      at_risk <- (obstimes[valid] > t_eval) | (obstimes[valid] <= t_eval & obsevents_numeric[valid] == 1)
      event_by_t <- obsevents_numeric[valid] == 1 & obstimes[valid] <= t_eval & at_risk
      surv_pred <- pred_surv[valid][at_risk]
      
      if (length(surv_pred) == 0) {
        return(NA_real_)
      }

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


#' Integrated Brier Score (IBS)
#'
#' Integrates Brier scores over the time range using trapezoidal rule.
#' Provides a single summary measure of prediction error across all time points.
#'
#' **Formula**:
#' \eqn{IBS = \frac{1}{t_{max} - t_{min}} \int_0^{t_{max}} BS(t) dt}
#'
#' Computed numerically via trapezoidal rule over the time grid.
#'
#' **Interpretation**:
#' - Range: 0 (perfect) to 1 (worst)
#' - Combines calibration across entire follow-up
#' - Weights later times equally (unlike summary C-index)
#'
#' @param predsurv Matrix of predicted survival probabilities
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param eval_times Optional evaluation times. If NULL, uses `pred_times`.
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return Scalar integrated Brier score value
#'
#' @references
#' Mogensen, U. B., Ishwaran, H., & Gerds, T. A. (2012).
#' "Evaluating random forests for survival analysis using prediction error curves."
#' *Journal of Statistical Software*, 50(11), 1–23.
#'
#' @keywords internal
#'
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
    rlang::abort("BrierScore output is missing the 'AppErr' component.")
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


#' Integrated Concordance Index
#'
#' Averages time-dependent C-index values over the time grid.
#' Provides a summary measure of discrimination across follow-up.
#'
#' **Interpretation**:
#' - Average of all C-index values at supplied time points
#' - Single summary value (unlike time-dependent C-index)
#' - Range: 0.5 to 1.0
#'
#' @param predsurv Matrix of predicted survival probabilities
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return Scalar integrated C-index value
#'
#' @keywords internal
#'
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
    rlang::abort("timedepConcordance output is missing the 'AppCindex' component.")
  }
  cindex_values <- cindex_obj$AppCindex$matrix

  # Calculate mean concordance over time (integrated C-index)
  integrated_c <- mean(cindex_values, na.rm = TRUE)

  return(integrated_c)
}
