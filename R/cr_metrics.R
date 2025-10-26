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
      stop(sprintf(
        "%s vector must have length equal to the number of observations (%d).",
        context, n_obs
      ))
    }
    Predictions <- matrix(Predictions, nrow = 1)
  } else {
    Predictions <- as.matrix(Predictions)
  }

  if (!is.matrix(Predictions)) {
    stop(sprintf("%s must be coercible to a matrix.", context))
  }

  mat <- Predictions

  if (ncol(mat) == n_obs) {
    # Already time-by-observation
  } else {
    stop(sprintf(
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
    stop(sprintf(
      "Unable to determine the time grid for %s. Provide 'pred_times' or supply numeric rownames/'Times' attribute.",
      tolower(context)
    ))
  }

  times <- suppressWarnings(as.numeric(times))
  if (any(is.na(times))) {
    stop(sprintf("'pred_times' must be numeric for %s.", tolower(context)))
  }

  if (length(times) != nrow(mat)) {
    stop(sprintf(
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
    stop(sprintf(
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

#' @title timedepConcordanceCR
#'
#' @description Get time dependent concordance for competing risk outcomes (only event 1)
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=times, cols=observations) for the cause of interest
#' @param time evaluation time point
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#' @param pred_times optional numeric vector of times corresponding to columns in `Predictions`
#'
#' @return concordance index value
#' @export
timedepConcordanceCR<-function(SurvObj, Predictions, time, cause=1, TestMat=NULL, pred_times = NULL){
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(time) || length(time) != 1) {
    stop("'time' must be a single numeric value")
  }
  if (missing(cause)) {
    stop("argument 'cause' is missing, with no default")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
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


#' @title BrierScoreCR
#'
#' @description Calculate Brier score for competing risk predictions at evaluation time(s)
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=times, cols=observations) for the cause of interest
#' @param eval_times evaluation time point(s). If a vector, returns a vector of scores.
#'   For backward compatibility, also accepts 'time' as a deprecated alias.
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#' @param pred_times optional numeric vector of times corresponding to rows in `Predictions`
#' @param time deprecated, use eval_times instead
#'
#' @return Brier score value(s). If eval_times is a vector, returns a vector; if scalar, returns a scalar.
#' @export
BrierScoreCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL, time = NULL) {
  
  # Handle backward compatibility: support both 'time' and 'eval_times'
  if (!is.null(time) && is.null(eval_times)) {
    # Warn user about deprecation
    message("Parameter 'time' is deprecated. Please use 'eval_times' instead.")
    eval_times <- time
  }
  
  if (is.null(eval_times)) {
    stop("'eval_times' must be specified (or 'time' for backward compatibility)")
  }
  
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(eval_times)) {
    stop("'eval_times' must be numeric")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
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

      event_indicator <- as.numeric(obsevents == cause & obstimes <= t & !is.na(obsevents))
      valid_idx <- !is.na(pred_at_time)
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

  event_indicator <- as.numeric(obsevents == cause & obstimes <= time & !is.na(obsevents))
  valid_idx <- !is.na(pred_at_time)
  if (!any(valid_idx)) {
    return(NA_real_)
  }

  mean((event_indicator[valid_idx] - pred_at_time[valid_idx])^2)
}


#' @title integratedConcordanceCR
#'
#' @description Calculate integrated concordance index for competing risk predictions over time range
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=times, cols=observations) for the cause of interest
#' @param eval_times vector of evaluation time points
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return integrated concordance index (scalar)
#' @export
integratedConcordanceCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
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


#' @title integratedBrierCR
#'
#' @description Calculate integrated Brier score for competing risk predictions over time range
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=times, cols=observations) for the cause of interest
#' @param eval_times vector of evaluation time points
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return integrated Brier score (scalar)
#' @export
integratedBrierCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
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
  brier_values <- sapply(eval_times, function(t) {
    BrierScoreCR(
      SurvObj = SurvObj,
      Predictions = pred_matrix,
      time = t,
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


#' @title restrictedMeanTimeLostCR
#'
#' @description Calculate Restricted Mean Time Lost (RMTL) for competing risk predictions
#' @param Predictions matrix of predicted CIFs (rows=times, cols=observations) for the cause of interest
#' @param times vector of time points corresponding to rows of Predictions
#' @param UL upper limit of integration
#' @param LL lower limit of integration (default: 0)
#'
#' @return vector of RMTL values for each observation
#' @export
restrictedMeanTimeLostCR <- function(Predictions, times, UL, LL = 0) {
  # Input validation
  if (!is.matrix(Predictions)) {
    stop("'Predictions' must be a matrix")
  }
  if (!is.numeric(times)) {
    stop("'times' must be numeric")
  }
  if (length(times) != nrow(Predictions)) {
    stop("'times' must match the number of rows in Predictions (rows represent time points)")
  }
  if (!is.numeric(UL) || length(UL) != 1 || UL <= LL) {
    stop("'UL' must be a single numeric value greater than 'LL'")
  }
  if (!is.numeric(LL) || length(LL) != 1 || LL < 0) {
    stop("'LL' must be a single non-negative numeric value")
  }

  # Ensure times are sorted
  if (!all(diff(times) >= 0)) {
    stop("'times' must be sorted in ascending order")
  }

  # Integrate CIF for each observation using trapezoidal rule
  time_lost_vector <- apply(Predictions, 2, function(cif_curve) {
    # Find indices within [LL, UL]
    valid_idx <- times >= LL & times <= UL
    if (sum(valid_idx) < 2) {
      return(NA)  # Need at least 2 points for integration
    }

    times_valid <- times[valid_idx]
    cif_valid <- cif_curve[valid_idx]

    # Trapezoidal integration
    time_diffs <- diff(times_valid)
    avg_cif <- (cif_valid[-1] + cif_valid[-length(cif_valid)]) / 2
    integral <- sum(time_diffs * avg_cif)

    return(integral)
  })

  return(time_lost_vector)
}
