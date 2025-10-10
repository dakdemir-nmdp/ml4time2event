.select_prediction_vector <- function(Predictions, eval_time, pred_times = NULL) {
  mat <- as.matrix(Predictions)
  if (is.null(pred_times)) {
    attr_times <- attr(Predictions, "Times")
    if (!is.null(attr_times)) {
      pred_times <- attr_times
    } else if (!is.null(colnames(mat))) {
      numeric_colnames <- suppressWarnings(as.numeric(colnames(mat)))
      if (!all(is.na(numeric_colnames))) {
        pred_times <- numeric_colnames
      }
    }
  }

  if (ncol(mat) == 1 && is.null(pred_times)) {
    return(list(pred = as.numeric(mat[, 1]), times = eval_time))
  }

  if (is.null(pred_times)) {
    stop("Unable to determine prediction time grid. Provide 'pred_times' or numeric column names/Times attribute.")
  }

  if (length(pred_times) != ncol(mat)) {
    stop("'pred_times' must have the same length as the number of columns in 'Predictions'.")
  }

  idx <- which.min(abs(pred_times - eval_time))
  list(pred = as.numeric(mat[, idx]), times = pred_times)
}

#' @title timedepConcordanceCR
#'
#' @description Get time dependent concordance for competing risk outcomes (only event 1)
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
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
  if (!is.matrix(Predictions) && !is.vector(Predictions) && !is.data.frame(Predictions)) {
    stop("'Predictions' must be a matrix, vector, or data.frame")
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

  if (is.data.frame(Predictions)) {
    Predictions <- as.matrix(Predictions)
  } else if (is.vector(Predictions)) {
    Predictions <- matrix(Predictions, ncol = 1)
  }

  pred_info <- .select_prediction_vector(Predictions, time, pred_times)
  pred_at_time <- pred_info$pred

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
#' @description Calculate Brier score for competing risk predictions at specific times
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
#' @param time evaluation time point
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return Brier score value
#' @param pred_times optional numeric vector of times corresponding to columns in `Predictions`
#' @export
BrierScoreCR <- function(SurvObj, Predictions, time, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.matrix(Predictions) && !is.vector(Predictions) && !is.data.frame(Predictions)) {
    stop("'Predictions' must be a matrix, vector, or data.frame")
  }
  if (!is.numeric(time) || length(time) != 1) {
    stop("'time' must be a single numeric value")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  obsevents <- SurvObj[, "status"]

  if (is.vector(Predictions)) {
    Predictions <- matrix(Predictions, ncol = 1)
  }

  pred_info <- .select_prediction_vector(Predictions, time, pred_times)
  pred_at_time <- pred_info$pred

  event_indicator <- as.numeric(obsevents == cause & obstimes <= time)
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
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
#' @param eval.times vector of evaluation time points
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return integrated concordance index (scalar)
#' @export
integratedConcordanceCR <- function(SurvObj, Predictions, eval.times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  Predictions <- as.matrix(Predictions)
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  # obsevents <- SurvObj[, "status"]  # Not used in this function

  if (is.null(pred_times)) {
    attr_times <- attr(Predictions, "Times")
    if (!is.null(attr_times)) {
      pred_times <- attr_times
    } else if (!is.null(colnames(Predictions))) {
      numeric_colnames <- suppressWarnings(as.numeric(colnames(Predictions)))
      if (!all(is.na(numeric_colnames))) {
        pred_times <- numeric_colnames
      }
    }
  }

  if (is.null(eval.times)) {
    if (!is.null(pred_times)) {
      eval.times <- pred_times
    } else {
      eval.times <- seq(min(obstimes), max(obstimes), length.out = ncol(Predictions))
    }
  }

  # Calculate concordance at each time point
  concordance_values <- sapply(eval.times, function(t) {
    timedepConcordanceCR(SurvObj, Predictions, t, cause, TestMat, pred_times = pred_times)
  })

  # Remove NA values
  valid_idx <- !is.na(concordance_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  concordance_values <- concordance_values[valid_idx]
  # eval_times_valid <- eval.times[valid_idx]  # Not used

  # Calculate integrated concordance (mean over time)
  integrated_c <- mean(concordance_values, na.rm = TRUE)

  return(integrated_c)
}


#' @title integratedBrierCR
#'
#' @description Calculate integrated Brier score for competing risk predictions over time range
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
#' @param eval.times vector of evaluation time points
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return integrated Brier score (scalar)
#' @export
integratedBrierCR <- function(SurvObj, Predictions, eval.times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  Predictions <- as.matrix(Predictions)
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  # obsevents <- SurvObj[, "status"]  # Not used in this function

  if (is.null(pred_times)) {
    attr_times <- attr(Predictions, "Times")
    if (!is.null(attr_times)) {
      pred_times <- attr_times
    } else if (!is.null(colnames(Predictions))) {
      numeric_colnames <- suppressWarnings(as.numeric(colnames(Predictions)))
      if (!all(is.na(numeric_colnames))) {
        pred_times <- numeric_colnames
      }
    }
  }

  if (is.null(eval.times)) {
    if (!is.null(pred_times)) {
      eval.times <- pred_times
    } else {
      eval.times <- seq(min(obstimes), max(obstimes), length.out = ncol(Predictions))
    }
  }

  # Calculate Brier score at each time point
  brier_values <- sapply(eval.times, function(t) {
    BrierScoreCR(SurvObj, Predictions, t, cause, TestMat, pred_times = pred_times)
  })

  # Remove NA values
  valid_idx <- !is.na(brier_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  brier_values <- brier_values[valid_idx]
  eval_times_valid <- eval.times[valid_idx]

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
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
#' @param times vector of time points corresponding to columns of Predictions
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
  if (!is.numeric(times) || length(times) != ncol(Predictions)) {
    stop("'times' must be numeric and match number of columns in Predictions")
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
  time_lost_vector <- apply(Predictions, 1, function(cif_curve) {
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
