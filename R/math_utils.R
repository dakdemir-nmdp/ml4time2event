#' @title Integrator
#' @description Integrate a curve (scores) over a specified time range using the trapezoidal rule.
#' @param times numeric vector of time points.
#' @param scores numeric vector of scores corresponding to time points.
#' @param minmax numeric vector of length 2 specifying the integration limits [min, max].
#' @param scale logical, if TRUE, scale the result by the length of the integration interval.
#' @return numeric value representing the integrated score.
#' @importFrom pracma trapz
#' @noRd
Integrator<-function(times, scores, minmax=c(1,35), scale=FALSE){
  # Ensure times and scores have the same length
  if (length(times) != length(scores)) {
    stop("Length of 'times' and 'scores' must be equal.")
  }
  # Ensure minmax is valid
  if (length(minmax) != 2 || minmax[1] >= minmax[2]) {
      stop("'minmax' must be a numeric vector of length 2 with minmax[1] < minmax[2].")
  }

  # Create a mask for times within the specified range
  mask<-(times>=minmax[1]) & (times<=minmax[2])

  # Filter times and scores based on the mask
  timesn<-times[mask]
  scoressn<-scores[mask]

  # Handle cases where no points fall within the range
  if (length(timesn) < 2) {
      warning("Less than 2 points found within the specified minmax range. Returning 0.")
      return(0)
  }

  # Sort points by time to ensure correct integration order
  order_idx <- order(timesn)
  timesn <- timesn[order_idx]
  scoressn <- scoressn[order_idx]

  # Calculate the area under the curve using the trapezoidal rule
  AUCsuperlearnMean = pracma::trapz(timesn, scoressn)

  # Scale the result if requested
  if (scale){
      interval_length <- minmax[2] - minmax[1]
      if (interval_length > 0) {
          AUCsuperlearnMean <- AUCsuperlearnMean / interval_length
      } else {
          warning("Integration interval length is zero or negative. Cannot scale.")
          # Return unscaled AUC or 0 depending on desired behavior
          # return(AUCsuperlearnMean)
          return(0) # Return 0 if interval is invalid for scaling
      }
  }
  AUCsuperlearnMean
}
