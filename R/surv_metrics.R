#' @title timedepConcordance
#'
#' @description Get time dependent concordance index (C-index) for survival models.
#' @param predsurv Predicted matrix of survival probabilities (rows=times, cols=observations).
#' @param predsurvtimes The times for which the probabilities are predicted.
#' @param obstimes Observed times vector.
#' @param obsevents Observed event indicator vector (0=censored, 1=event).
#' @param ctimes Optional numeric vector of times at which to calculate the C-index. If NULL, uses predsurvtimes.
#'
#' @return Output from the 'pec::cindex' function, containing C-index values over time.
#' @importFrom pec cindex
#' @importFrom survival Surv
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
timedepConcordance<-function(predsurv, predsurvtimes, obstimes, obsevents, ctimes=NULL){
  # Ensure event is numeric 0/1
  obsevents_numeric <- as.numeric(obsevents == 1)
  # Create data frame for pec
  datforpec<-data.frame(time=obstimes, event=obsevents_numeric)

  # Determine evaluation times
  if (is.null(ctimes)){
    eval_times <- predsurvtimes
  } else {
    eval_times <- ctimes
  }
  # Ensure times are sorted and unique, and >= 0
  eval_times <- sort(unique(eval_times[eval_times >= 0]))
  if (length(eval_times) == 0) {
      stop("No valid evaluation times (ctimes >= 0) provided or derived.")
  }

  # Calculate C-index using pec::cindex
  # pec::cindex expects predictions matrix with rows=observations, cols=times
  # Input predsurv is rows=times, cols=observations, so transpose it.
  cindexTest<-pec::cindex(object=t(predsurv),
                          formula=Surv(time, event)~1, # Formula for baseline C-index
                          data=datforpec,
                          eval.times=eval_times)
  cindexTest
}


#' @title integratedC
#'
#' @description Calculate the integrated C-index over a specified time range.
#'
#' @param times The times at which the C-index scores were calculated.
#' @param scores Time-dependent C-index scores corresponding to 'times'.
#' @param minmax Numeric vector of length 2 specifying the integration limits [min, max].
#'
#' @return A scalar value representing the integrated C-index, scaled by the interval length.
#' @export
integratedC<-function(times, scores, minmax=c(1,35)){
  # Assuming Integrator function is loaded/available
  AUCMean <- Integrator(times = times, scores = scores, minmax = minmax, scale = TRUE)
  AUCMean
}




#' @title BrierScore
#'
#' @description Calculate the time-dependent Brier score using pec::pec.
#'
#' @param predsurv Predicted matrix of survival probabilities (rows=times, cols=observations).
#' @param predsurvtimes The times for which the probabilities are predicted.
#' @param obstimes Observed times vector.
#' @param obsevents Observed event indicator vector (0=censored, 1=event).
#'
#' @return The output object from the 'pec::pec' function, containing Brier scores over time.
#' @importFrom pec pec
#' @importFrom survival Surv
#' @export
BrierScore<-function(predsurv, predsurvtimes, obstimes, obsevents){
  # Ensure event is numeric 0/1
  obsevents_numeric <- as.numeric(obsevents == 1)
  # Create data frame for pec
  datforpec<-data.frame(time=obstimes, event=obsevents_numeric)

  # Ensure times are sorted and unique, >= 0
  eval_times <- sort(unique(predsurvtimes[predsurvtimes >= 0]))
   if (length(eval_times) == 0) {
      stop("No valid evaluation times (predsurvtimes >= 0) provided.")
  }
   # Ensure predsurv matches eval_times (subset rows if necessary)
   time_indices <- match(eval_times, predsurvtimes)
   predsurv_eval <- predsurv[time_indices, , drop = FALSE]


  # Calculate Brier score using pec::pec
  # pec::pec expects predictions matrix with rows=observations, cols=times
  # Input predsurv_eval is rows=times, cols=observations, so transpose it.
  BrierTest<-pec::pec(object=list(model=t(predsurv_eval)), # Needs to be a named list
                      formula=survival::Surv(time, event)~1,
                      data=datforpec,
                      times = eval_times,
                      exact = FALSE, # Use approximation for speed
                      reference = FALSE) # Don't compute reference model Brier score
  BrierTest
}



#' @title integratedBrier
#'
#' @description Calculate the integrated Brier score (IBS) over a specified time range.
#'
#' @param times The times at which the Brier scores were calculated.
#' @param scores Time-dependent Brier scores corresponding to 'times'.
#' @param minmax Numeric vector of length 2 specifying the integration limits [min, max].
#'
#' @return A scalar value representing the integrated Brier score, scaled by the interval length.
#'
#' @export
integratedBrier<-function(times, scores, minmax=c(1,35)){
  # Assuming Integrator function is loaded/available
  AUCMean <- Integrator(times = times, scores = scores, minmax = minmax, scale = TRUE)
  AUCMean
}
