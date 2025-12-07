.cr_prepare_prediction_matrix <- function(Predictions,
                                          n_obs,
                                          pred_times = NULL,
                                          default_time = NULL,
                                          context = "Predictions") {

  original_obj <- Predictions

  if (is.data.frame(Predictions)) {
    Predictions <- as.matrix(Predictions)
  } else if (is.vector(Predictions)) {
    if (length(Predictions) != n_obs) {
      rlang::abort(sprintf(
        "%s vector must have length equal to the number of observations (%d).",
        context, n_obs
      ))
    }
    Predictions <- matrix(Predictions, nrow = 1)
  } else {
    Predictions <- as.matrix(Predictions)
  }

  if (!is.matrix(Predictions)) {
    rlang::abort(sprintf("%s must be coercible to a matrix.", context))
  }

  mat <- Predictions

  if (ncol(mat) == n_obs) {
    # Already time-by-observation
  } else {
    rlang::abort(sprintf(
      "%s must have rows=times and columns=observations. Got %dx%d, expected nrow=%d (times), ncol=%d (observations).",
      context, nrow(mat), ncol(mat), length(pred_times) %||% nrow(mat), n_obs
    ))
  }

  times <- pred_times
  if (is.null(times)) {
    times <- attr(original_obj, "Times")
  }
  if (is.null(times)) {
    times <- attr(original_obj, "times")
  }
  if (is.null(times)) {
    rownames_numeric <- suppressWarnings(as.numeric(rownames(mat)))
    if (!any(is.na(rownames_numeric))) {
      times <- rownames_numeric
    }
  }
  if (is.null(times)) {
    if (nrow(mat) == 1 && !is.null(default_time)) {
      times <- default_time
    }
  }
  if (is.null(times)) {
    rlang::abort(sprintf(
      "Unable to determine the time grid for %s. Provide 'pred_times' or supply numeric rownames/'Times' attribute.",
      tolower(context)
    ))
  }

  times <- suppressWarnings(as.numeric(times))
  if (any(is.na(times))) {
    rlang::abort(sprintf("'pred_times' must be numeric for %s.", tolower(context)))
  }

  if (length(times) != nrow(mat)) {
    rlang::abort(sprintf(
      "'pred_times' must have length %d to match the number of rows in %s (rows represent time points).",
      nrow(mat), tolower(context)
    ))
  }

  list(matrix = mat, times = times)
}

.cr_select_prediction_row <- function(pred_matrix,
                                      pred_times,
                                      eval_time,
                                      context = "Predictions") {
  if (is.null(pred_times) || length(pred_times) != nrow(pred_matrix)) {
    rlang::abort(sprintf(
      "Prediction times for %s must be provided and match the number of matrix rows.",
      tolower(context)
    ))
  }

  idx <- which.min(abs(pred_times - eval_time))
  list(
    pred = as.numeric(pred_matrix[idx, ]),
    matched_time = pred_times[idx],
    times = pred_times
  )
}

#' Time-Dependent Concordance Index for Competing Risks
#'
#' Computes cause-specific concordance index for competing risks data.
#' Measures how well predicted cumulative incidence functions rank pairs
#' of observations where one had the competing event of interest.
#'
#' **Formula**:
#' For cause \eqn{j} at time \eqn{t}:
#' \eqn{C_j(t) = \frac{\sum_{i,j} I(t_i < t_j) \cdot I(\hat{CIF}_{j,i}(t) > \hat{CIF}_{j,j}(t)) \cdot I(\delta_j = j)}{\sum_{i,j} I(t_i < t_j) \cdot I(\delta_j = j)}}
#'
#' **Interpretation**:
#' - Range: 0.5 (random) to 1.0 (perfect)
#' - Computed only among pairs where the event of interest occurred
#' - Competing events exclude observations from risk set
#'
#' @param SurvObj Surv object with structure: `Surv(time, status=cause)`
#' @param Predictions Matrix of predicted CIFs (rows = times, cols = observations)
#' @param time Scalar evaluation time
#' @param cause Integer cause of interest (default: 1)
#' @param TestMat Optional test matrix (currently unused)
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar C-index value at specified time, or NA if no pairs comparable
#'
#' @references
#' Fine, J. P., & Gray, R. J. (1999). "A proportional hazards model for the
#' subdistribution of a competing risk." *Journal of the American Statistical Association*,
#' 94(446), 496–509.
#'
#' @keywords internal
#'
timedepConcordanceCR<-function(SurvObj, Predictions, time, cause=1, TestMat=NULL, pred_times = NULL){
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(time) || length(time) != 1) {
    rlang::abort("'time' must be a single numeric value")
  }
  if (missing(cause)) {
    rlang::abort("argument 'cause' is missing, with no default")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  obsevents <- SurvObj[, "status"]

  n_obs <- length(obstimes)
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    default_time = time,
    context = "Predictions"
  )
  pred_row <- .cr_select_prediction_row(
    pred_matrix = pred_prep$matrix,
    pred_times = pred_prep$times,
    eval_time = time,
    context = "Predictions"
  )
  pred_at_time <- pred_row$pred

  event_idx <- which(obsevents == cause & obstimes <= time & !is.na(obsevents))
  if (length(event_idx) == 0) {
    return(NA_real_)
  }

  concordant <- 0
  total_pairs <- 0

  for (i in event_idx) {
    if (is.na(pred_at_time[i])) next
    for (j in seq_along(obstimes)) {
      if (i == j) next
      if (is.na(pred_at_time[j])) next
      if (obstimes[j] < obstimes[i]) next
      if (obsevents[j] == cause && obstimes[j] <= obstimes[i]) next
      total_pairs <- total_pairs + 1
      # Standard concordance: higher risk = earlier event
      if (pred_at_time[i] > pred_at_time[j]) {
        concordant <- concordant + 1
      } else if (pred_at_time[i] == pred_at_time[j]) {
        concordant <- concordant + 0.5
      }
    }
  }

  if (total_pairs == 0) {
    return(NA_real_)
  }

  concordant / total_pairs
}


#' Brier Score for Competing Risks Predictions
#'
#' Computes mean squared error between predicted cumulative incidence functions
#' and observed event status for a specified cause. Accounts for competing events.
#'
#' **Formula**:
#' \eqn{BS_j(t) = E[(CIF_j(t|X) - Y_j(t))^2]}
#'
#' where \eqn{CIF_j(t|X)} is predicted CIF for cause \eqn{j} and \eqn{Y_j(t)}
#' is binary indicator of cause \eqn{j} occurring by time \eqn{t}.
#'
#' **Interpretation**:
#' - Range: 0 (perfect) to 1 (worst)
#' - Combines calibration and discrimination
#' - Computed only for observations at risk (event of interest or censored/competing)
#'
#' @param SurvObj Surv object with competing risks structure
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Scalar or vector of evaluation time(s)
#' @param cause Integer cause of interest (default: 1)
#' @param TestMat Optional test matrix (currently unused)
#' @param pred_times Numeric vector of time grid from predictions
#' @param time Deprecated parameter for backward compatibility; use `eval_times`
#'
#' @return Numeric Brier score at specified time(s)
#'
#' @keywords internal
#'
BrierScoreCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL, time = NULL) {
  
  # Handle backward compatibility: support both 'time' and 'eval_times'
  if (!is.null(time) && is.null(eval_times)) {
    # Warn user about deprecation
    message("Parameter 'time' is deprecated. Please use 'eval_times' instead.")
    eval_times <- time
  }
  
  if (is.null(eval_times)) {
    rlang::abort("'eval_times' must be specified (or 'time' for backward compatibility)")
  }
  
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(eval_times)) {
    rlang::abort("'eval_times' must be numeric")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  obsevents <- SurvObj[, "status"]

  n_obs <- length(obstimes)
  
  # Handle vector of evaluation times - compute Brier score for each
  if (length(eval_times) > 1) {
    brier_scores <- sapply(eval_times, function(t) {
      pred_prep <- .cr_prepare_prediction_matrix(
        Predictions = Predictions,
        n_obs = n_obs,
        pred_times = pred_times,
        default_time = t,
        context = "Predictions"
      )
      pred_row <- .cr_select_prediction_row(
        pred_matrix = pred_prep$matrix,
        pred_times = pred_prep$times,
        eval_time = t,
        context = "Predictions"
      )
      pred_at_time <- pred_row$pred

      # Only include observations at risk at time t
      # At risk if: time > t (still at risk) OR (time <= t AND event == cause) (event of interest occurred)
      # Note: Those with competing events before t are not in the risk set
      at_risk <- (obstimes > t) | (obstimes <= t & obsevents == cause & !is.na(obsevents))
      event_indicator <- as.numeric(obsevents == cause & obstimes <= t & !is.na(obsevents))
      valid_idx <- at_risk & !is.na(pred_at_time)
      
      if (!any(valid_idx)) {
        return(NA_real_)
      }
      mean((event_indicator[valid_idx] - pred_at_time[valid_idx])^2)
    })
    return(brier_scores)
  }
  
  # Single evaluation time - return scalar
  time <- eval_times[1]
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    default_time = time,
    context = "Predictions"
  )
  pred_row <- .cr_select_prediction_row(
    pred_matrix = pred_prep$matrix,
    pred_times = pred_prep$times,
    eval_time = time,
    context = "Predictions"
  )
  pred_at_time <- pred_row$pred

  # Only include observations at risk at time t
  # At risk if: time > t (still at risk) OR (time <= t AND event == cause) (event of interest occurred)
  # Note: Those with competing events before t are not in the risk set
  at_risk <- (obstimes > time) | (obstimes <= time & obsevents == cause & !is.na(obsevents))
  event_indicator <- as.numeric(obsevents == cause & obstimes <= time & !is.na(obsevents))
  valid_idx <- at_risk & !is.na(pred_at_time)
  
  if (!any(valid_idx)) {
    return(NA_real_)
  }

  mean((event_indicator[valid_idx] - pred_at_time[valid_idx])^2)
}


#' Integrated Concordance Index for Competing Risks
#'
#' Averages time-dependent concordance across evaluation times.
#' Provides a single summary measure of discrimination over follow-up.
#'
#' **Interpretation**:
#' - Mean of C-index values across time points
#' - Range: 0.5 (random) to 1.0 (perfect)
#' - Accounts for competing events in all comparisons
#'
#' @param SurvObj Surv object with competing risks
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Vector of time points at which to evaluate
#' @param cause Integer cause of interest
#' @param TestMat Optional test matrix
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar integrated C-index value
#'
#' @keywords internal
#'
integratedConcordanceCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]

  n_obs <- length(obstimes)
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    context = "Predictions"
  )
  pred_matrix <- pred_prep$matrix
  pred_times <- pred_prep$times

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }

  # Calculate concordance at each time point
  concordance_values <- sapply(eval_times, function(t) {
    timedepConcordanceCR(
      SurvObj = SurvObj,
      Predictions = pred_matrix,
      time = t,
      cause = cause,
      TestMat = TestMat,
      pred_times = pred_times
    )
  })

  # Remove NA values
  valid_idx <- !is.na(concordance_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  concordance_values <- concordance_values[valid_idx]
  # eval_times_valid <- eval_times[valid_idx]  # Not used

  # Calculate integrated concordance (mean over time)
  integrated_c <- mean(concordance_values, na.rm = TRUE)

  return(integrated_c)
}


#' Integrated Brier Score for Competing Risks
#'
#' Integrates cause-specific Brier scores over time using trapezoidal rule.
#' Provides a single summary measure of calibration and discrimination.
#'
#' **Formula**:
#' \eqn{IBS_j = \frac{1}{t_{max} - t_{min}} \int_0^{t_{max}} BS_j(t) dt}
#'
#' @param SurvObj Surv object with competing risks
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Vector of evaluation times. If NULL, uses time grid from predictions.
#' @param cause Integer cause of interest
#' @param TestMat Optional test matrix
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar integrated Brier score
#'
#' @keywords internal
#'
integratedBrierCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]

  n_obs <- length(obstimes)
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    context = "Predictions"
  )
  pred_matrix <- pred_prep$matrix
  pred_times <- pred_prep$times

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }

  # Calculate Brier score at each time point
  # Use the non-deprecated 'eval_times' argument to avoid warnings
  brier_values <- sapply(eval_times, function(t) {
    BrierScoreCR(
      SurvObj = SurvObj,
      Predictions = pred_matrix,
      eval_times = t,
      cause = cause,
      TestMat = TestMat,
      pred_times = pred_times
    )
  })

  # Remove NA values
  valid_idx <- !is.na(brier_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  brier_values <- brier_values[valid_idx]
  eval_times_valid <- eval_times[valid_idx]

  # Calculate integrated Brier score (area under curve using trapezoidal rule)
  if (length(eval_times_valid) < 2) {
    return(brier_values[1])
  }

  time_diffs <- diff(eval_times_valid)
  avg_brier <- (brier_values[-1] + brier_values[-length(brier_values)]) / 2
  integrated_bs <- sum(time_diffs * avg_brier) / (max(eval_times_valid) - min(eval_times_valid))

  return(integrated_bs)
}


#'
#'
restrictedMeanTimeLostCR <- function(Predictions, times, UL, LL = 0) {
  # Input validation
  if (!is.matrix(Predictions)) {
    rlang::abort("'Predictions' must be a matrix")
  }
  if (!is.numeric(times)) {
    rlang::abort("'times' must be numeric")
  }
  if (length(times) != nrow(Predictions)) {
    rlang::abort("'times' must match the number of rows in Predictions (rows represent time points)")
  }
  if (!is.numeric(UL) || length(UL) != 1 || UL <= LL) {
    rlang::abort("'UL' must be a single numeric value greater than 'LL'")
  }
  if (!is.numeric(LL) || length(LL) != 1 || LL < 0) {
    rlang::abort("'LL' must be a single non-negative numeric value")
  }

  # Ensure times are sorted
  if (!all(diff(times) >= 0)) {
    rlang::abort("'times' must be sorted in ascending order")
  }

  # Calculate ETL using unified function
  time_lost_vector <- calculate_expected_time_lost(
    times = times,
    event_probs = Predictions,
    upper_limit = UL,
    lower_limit = LL
  )

  return(time_lost_vector)
}
