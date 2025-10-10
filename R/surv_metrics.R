#' @title timedepConcordance
#'
#' @description Calculate time-dependent concordance for survival predictions
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations)
#' @param predsurvtimes The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param TestMat Optional test dataset (used if formula needs predictors)
#'
#' @return Output from the pec::cindex function
#' @importFrom pec cindex
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
timedepConcordance <- function(predsurv, predsurvtimes, obstimes, obsevents, TestMat = NULL) {

  # Input validation
  if (!is.matrix(predsurv)) {
    stop("'predsurv' must be a matrix")
  }
  if (!is.numeric(predsurvtimes)) {
    stop("'predsurvtimes' must be numeric")
  }
  if (!is.numeric(obstimes)) {
    stop("'obstimes' must be numeric")
  }
  if (!is.numeric(obsevents) && !is.logical(obsevents)) {
    stop("'obsevents' must be numeric or logical")
  }

  # Ensure obsevents is numeric (0=censored, 1=event)
  obsevents_numeric <- as.numeric(obsevents)

  # Create data frame for pec
  datforpec <- data.frame(time = obstimes, event = obsevents_numeric)

  # Add predictors from TestMat if provided
  if (!is.null(TestMat)) {
    datforpec <- cbind(datforpec, TestMat)
    formula_str <- "Hist(time, event) ~ ."
  } else {
    formula_str <- "Hist(time, event) ~ 1"
  }

  # Calculate c-index using pec::cindex
  # Note: pec expects a matrix with rows=observations, cols=times
  # Input predsurv is rows=times, cols=observations, so we need to transpose
  pred_for_pec <- t(predsurv)

  cindexTest <- pec::cindex(
    object = pred_for_pec,
    formula = stats::as.formula(formula_str),
    data = datforpec,
    eval.times = predsurvtimes,
    cens.model = "marginal"
  )

  return(cindexTest)
}


#' @title BrierScore
#'
#' @description Calculate Brier score for survival predictions at specific times
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations)
#' @param predsurvtimes The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param eval.times Optional vector of evaluation times (defaults to predsurvtimes)
#' @param TestMat Optional test dataset (used if formula needs predictors)
#'
#' @return Output from the pec::pec function containing Brier scores
#' @importFrom pec pec
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
BrierScore <- function(predsurv, predsurvtimes, obstimes, obsevents,
                       eval.times = NULL, TestMat = NULL) {

  # Input validation
  if (!is.matrix(predsurv)) {
    stop("'predsurv' must be a matrix")
  }
  if (!is.numeric(predsurvtimes)) {
    stop("'predsurvtimes' must be numeric")
  }
  if (!is.numeric(obstimes)) {
    stop("'obstimes' must be numeric")
  }
  if (!is.numeric(obsevents) && !is.logical(obsevents)) {
    stop("'obsevents' must be numeric or logical")
  }

  # Use predsurvtimes if eval.times not specified
  if (is.null(eval.times)) {
    eval.times <- predsurvtimes
  }

  # Ensure obsevents is numeric (0=censored, 1=event)
  obsevents_numeric <- as.numeric(obsevents)

  # Create data frame for pec
  datforpec <- data.frame(time = obstimes, event = obsevents_numeric)

  # Add predictors from TestMat if provided
  if (!is.null(TestMat)) {
    datforpec <- cbind(datforpec, TestMat)
    formula_str <- "Hist(time, event) ~ ."
  } else {
    formula_str <- "Hist(time, event) ~ 1"
  }

  # Create a list object that pec can use
  # pec expects a list with survfit-like structure
  # Need to add required attributes for survfit objects
  pred_obj <- list(
    time = predsurvtimes,
    surv = predsurv,
    type = "right",
    n = ncol(predsurv),
    n.event = rep(0, length(predsurvtimes)),
    n.censor = rep(0, length(predsurvtimes)),
    n.risk = rep(ncol(predsurv), length(predsurvtimes))
  )
  class(pred_obj) <- c("survfit", "list")

  # Calculate Brier score using pec::pec
  brier_result <- pec::pec(
    object = list(model = pred_obj),
    formula = stats::as.formula(formula_str),
    data = datforpec,
    times = eval.times,
    exact = FALSE,
    cens.model = "marginal",
    splitMethod = "none",
    B = 0,
    verbose = FALSE
  )

  # Always extract Brier scores for requested eval.times
  if (!is.null(brier_result$AppErr)) {
    # brier_result$time gives the times for which Brier scores are available
    # brier_result$AppErr$model gives the Brier scores at those times
    brier_times <- as.numeric(brier_result$time)
    brier_values <- as.numeric(brier_result$AppErr$model)
    # Match requested eval.times to available times (allowing for floating point imprecision)
    req_times <- as.numeric(eval.times)
    idx <- vapply(req_times, function(t) {
      which.min(abs(brier_times - t))
    }, integer(1))
    # If the closest time is not close enough, set NA
    tol <- 1e-8
    matched <- abs(brier_times[idx] - req_times) < tol | brier_times[idx] == req_times
    out <- rep(NA_real_, length(req_times))
    out[matched] <- brier_values[idx[matched]]
    names(out) <- as.character(req_times)
    # Return scalar if only one time requested
    if (length(out) == 1) return(out[[1]])
    return(out)
  }
  stop("Could not extract Brier scores from pec object")
}


#' @title integratedBrier
#'
#' @description Calculate integrated Brier score over time range
#' @param predsurv Predicted survival probability matrix (rows=times, cols=observations)
#' @param predsurvtimes The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param eval.times Optional vector of evaluation times for integration
#' @param TestMat Optional test dataset
#'
#' @return Integrated Brier score (scalar)
#' @export
integratedBrier <- function(predsurv, predsurvtimes, obstimes, obsevents,
                            eval.times = NULL, TestMat = NULL) {

  # Get Brier scores at all time points
  brier_obj <- BrierScore(
    predsurv = predsurv,
    predsurvtimes = predsurvtimes,
    obstimes = obstimes,
    obsevents = obsevents,
    eval.times = eval.times,
    TestMat = TestMat
  )

  # Extract the Brier score values
  # pec objects have AppErr slot with the prediction error
  if (!is.null(brier_obj$AppErr)) {
    brier_values <- brier_obj$AppErr$model
    eval_times <- brier_obj$time
  } else {
    stop("Could not extract Brier scores from pec object")
  }

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
#' @param predsurvtimes The times for which the survival probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event)
#' @param TestMat Optional test dataset
#'
#' @return Integrated concordance index (scalar)
#' @export
integratedC <- function(predsurv, predsurvtimes, obstimes, obsevents, TestMat = NULL) {

  # Get time-dependent concordance
  cindex_obj <- timedepConcordance(
    predsurv = predsurv,
    predsurvtimes = predsurvtimes,
    obstimes = obstimes,
    obsevents = obsevents,
    TestMat = TestMat
  )

  # Extract concordance values at each time point
  if (!is.null(cindex_obj$AppCindex)) {
    cindex_values <- cindex_obj$AppCindex$matrix
    eval_times <- cindex_obj$time
  } else {
    stop("Could not extract concordance values from pec object")
  }

  # Calculate mean concordance over time (integrated C-index)
  integrated_c <- mean(cindex_values, na.rm = TRUE)

  return(integrated_c)
}
