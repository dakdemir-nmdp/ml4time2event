#' @title timedepConcordanceCR
#'
#' @description Get time dependent concordance for competing risk outcomes (only event 1)
#' @param predCIF Predicted matrix of probabilities (CIFs, rows=observations, cols=times)
#' @param predCIFtimes The times for which the probabilities are predicted
#' @param obstimes Observed times vector
#' @param obsevents Observed event indicator vector (0=censored, 1=event of interest, 2=competing event)
#' @param TestMat test dataset (optional, used if formula needs predictors)
#'
#' @return output from the "pec::cindex' function.
#' @importFrom pec cindex
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export
timedepConcordanceCR<-function(predCIF, predCIFtimes, obstimes, obsevents, TestMat=NULL){
  # Ensure obsevents is numeric (0, 1, 2)
  obsevents_numeric <- as.numeric(obsevents)
  # Create data frame for pec, ensuring correct column names
  datforpec<-data.frame(time=obstimes, event=obsevents_numeric)
  # Add predictors from TestMat if provided and needed by formula
  if (!is.null(TestMat)) {
      datforpec <- cbind(datforpec, TestMat)
      formula_str <- "Hist(time, event)~." # Use all predictors in TestMat
  } else {
      formula_str <- "Hist(time, event)~1" # Use intercept only if no predictors
  }

  # Calculate c-index using pec::cindex
  # Note: pec::cindex expects predictions matrix with rows=times, cols=observations
  # The input predCIF is rows=observations, cols=times, so we transpose it.
  cindexTest<-pec::cindex(object=t(predCIF),
                          formula=stats::as.formula(formula_str),
                          data=datforpec,
                          eval.times=predCIFtimes,
                          cause = 1) # Specify cause of interest for CR
  cindexTest
}
