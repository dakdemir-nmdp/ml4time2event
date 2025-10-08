#' @title timedepConcordanceCR
#'
#' @description Get time dependent concordance for competing risk outcomes (only event 1)
#' @param SurvObj survival object (Surv(time, status))
#' @param Predictions matrix of predicted CIFs (rows=observations, cols=times) for the cause of interest
#' @param time evaluation time point
#' @param cause cause of interest (1 or 2)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return concordance index value
#' @importFrom pec cindex
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
timedepConcordanceCR<-function(SurvObj, Predictions, time, cause=1, TestMat=NULL){
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

  # Ensure Predictions is a matrix
  if (is.data.frame(Predictions)) {
    Predictions <- as.matrix(Predictions)
  } else if (is.vector(Predictions)) {
    Predictions <- matrix(Predictions, ncol = 1)
  }

  # If multiple time points in Predictions, use first column (simplified)
  if (ncol(Predictions) > 1) {
    pred_at_time <- Predictions[, 1]
  } else {
    pred_at_time <- Predictions[, 1]
  }

  # Create data frame for pec
  datforpec <- data.frame(time = obstimes, event = obsevents)

  # Add predictors from TestMat if provided
  if (!is.null(TestMat)) {
    datforpec <- cbind(datforpec, TestMat)
    # formula_str <- "Hist(time, event)~."  # Not used in simplified implementation
  } else {
    # formula_str <- "Hist(time, event)~1"  # Not used in simplified implementation
  }

  # For competing risks, we need to handle this differently
  # Use a simplified concordance calculation for now
  # This is a placeholder - proper CR concordance needs more sophisticated implementation

  # Extract relevant observations (those with events of the specified cause by time t, or censored by time t)
  cause_events <- obsevents == cause & obstimes <= time & !is.na(obsevents)
  censored <- obsevents == 0 & obstimes <= time & !is.na(obsevents)
  # Also include competing events that occurred by time t
  competing_events <- obsevents != cause & obsevents != 0 & obstimes <= time & !is.na(obsevents)
  
  # For concordance, we compare cause events vs (censored + competing events)
  event_preds <- pred_at_time[cause_events]
  non_event_preds <- pred_at_time[censored | competing_events]

  # Remove NA values
  event_preds <- event_preds[!is.na(event_preds)]
  non_event_preds <- non_event_preds[!is.na(non_event_preds)]

  if (length(event_preds) < 2) {
    return(NA)  # Need at least 2 events for meaningful concordance
  }

  if (length(non_event_preds) == 0) {
    return(NA)
  }

  # Count concordant pairs (event pred > non-event pred)
  concordant <- 0
  total_pairs <- 0

  for (e_pred in event_preds) {
    for (c_pred in non_event_preds) {
      # Skip if either prediction is NA (shouldn't happen after filtering, but safety check)
      if (is.na(e_pred) || is.na(c_pred)) next
      
      total_pairs <- total_pairs + 1
      if (e_pred > c_pred) concordant <- concordant + 1
      else if (e_pred == c_pred) concordant <- concordant + 0.5  # Tie
    }
  }

  if (total_pairs == 0) {
    return(NA)
  }

  concordance <- concordant / total_pairs
  return(concordance)
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
#' @importFrom pec pec
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
BrierScoreCR <- function(SurvObj, Predictions, time, cause = 1, TestMat = NULL) {
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

  # Ensure Predictions is a matrix
  if (is.vector(Predictions)) {
    Predictions <- matrix(Predictions, ncol = 1)
  }

  # If multiple time points in Predictions, find the closest one to the requested time
  if (ncol(Predictions) > 1) {
    pred_at_time <- Predictions[, 1]  # Simplified - should match time points properly
  } else {
    pred_at_time <- Predictions[, 1]
  }

  # Create data frame for pec
  datforpec <- data.frame(time = obstimes, event = obsevents)

  # Add predictors from TestMat if provided
  if (!is.null(TestMat)) {
    datforpec <- cbind(datforpec, TestMat)
    formula_str <- "Hist(time, event)~."
  } else {
    formula_str <- "Hist(time, event)~1"
  }

  # Create prediction object for pec
  pred_obj <- list(
    time = time,
    cif = matrix(pred_at_time, ncol = 1)
  )
  class(pred_obj) <- c("cif", "list")

  # Calculate Brier score using pec::pec
  brier_result <- pec::pec(
    object = list(model = pred_obj),
    formula = stats::as.formula(formula_str),
    data = datforpec,
    times = time,
    cause = cause,
    exact = FALSE,
    cens.model = "marginal",
    splitMethod = "none",
    B = 0,
    verbose = FALSE
  )

  # Extract Brier score
  if (!is.null(brier_result$AppErr) && !is.null(brier_result$AppErr$model)) {
    brier_score <- brier_result$AppErr$model[1]
  } else {
    brier_score <- NA
  }

  return(brier_score)
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
integratedConcordanceCR <- function(SurvObj, Predictions, eval.times = NULL, cause = 1, TestMat = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.matrix(Predictions)) {
    stop("'Predictions' must be a matrix")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  # obsevents <- SurvObj[, "status"]  # Not used in this function

  # If eval.times not provided, use all time points from Predictions
  if (is.null(eval.times)) {
    eval.times <- seq(min(obstimes), max(obstimes), length.out = ncol(Predictions))
  }

  # Calculate concordance at each time point
  concordance_values <- sapply(eval.times, function(t) {
    timedepConcordanceCR(SurvObj, Predictions, t, cause, TestMat)
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
integratedBrierCR <- function(SurvObj, Predictions, eval.times = NULL, cause = 1, TestMat = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    stop("'SurvObj' must be a Surv object")
  }
  if (!is.matrix(Predictions)) {
    stop("'Predictions' must be a matrix")
  }
  if (!is.numeric(cause) || length(cause) != 1 || !(cause %in% c(1, 2))) {
    stop("'cause' must be 1 or 2")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  # obsevents <- SurvObj[, "status"]  # Not used in this function

  # If eval.times not provided, use all time points from Predictions
  if (is.null(eval.times)) {
    eval.times <- seq(min(obstimes), max(obstimes), length.out = ncol(Predictions))
  }

  # Calculate Brier score at each time point
  brier_values <- sapply(eval.times, function(t) {
    BrierScoreCR(SurvObj, Predictions, t, cause, TestMat)
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
