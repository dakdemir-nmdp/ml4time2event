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

  # Number of observations and evaluation times
  n_obs <- length(obstimes)
  n_times <- length(eval.times)

  # Initialize Brier score vector
  brier_scores <- numeric(n_times)

  # For each evaluation time
  for (t_idx in seq_len(n_times)) {
    t <- eval.times[t_idx]

    tryCatch({
      # Get predicted survival probabilities at time t for all observations
      # predsurv is times x observations, so find the closest time
      time_diffs <- abs(predsurvtimes - t)
      closest_time_idx <- which.min(time_diffs)
      pred_surv_at_t <- predsurv[closest_time_idx, ]

      # Calculate Brier score components for each observation
      brier_components <- numeric(n_obs)

      for (i in seq_len(n_obs)) {
        pred_surv_i <- pred_surv_at_t[i]
        obs_time_i <- obstimes[i]
        obs_event_i <- obsevents_numeric[i]

        # For right-censored data, the Brier score at time t is:
        # If individual had event by time t: (Ŝ(t))^2
        # If individual censored or event after t: (1 - Ŝ(t))^2
        if (obs_event_i == 1 && obs_time_i <= t) {
          # Event occurred by time t
          brier_components[i] <- pred_surv_i^2
        } else {
          # No event by time t (censored or event later)
          brier_components[i] <- (1 - pred_surv_i)^2
        }
      }

      # Average Brier score at time t
      brier_scores[t_idx] <- mean(brier_components)
    }, error = function(e) {
      warning("Error calculating Brier score at time ", t, ": ", e$message)
      brier_scores[t_idx] <- NA
    })
  }

  # Create a pec-like object structure for compatibility
  result <- list(
    AppErr = list(
      model = brier_scores,
      time = eval.times
    ),
    time = eval.times,
    call = match.call()
  )

  class(result) <- "pec"
  return(result)
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
  } else {
    stop("Could not extract concordance values from pec object")
  }

  # Calculate mean concordance over time (integrated C-index)
  integrated_c <- mean(cindex_values, na.rm = TRUE)

  return(integrated_c)
}
