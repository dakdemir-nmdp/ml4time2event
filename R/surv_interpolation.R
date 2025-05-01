#' @title survivalProbsInterpolator
#' @description Interpolate survival probabilities for new times.
#' @param x new times for interpolation
#' @param probs vector of survival probabilities
#' @param times vector of times corresponding to probs
#' @return interpolated survival probability values at times x
#' @importFrom stats approxfun
#' @noRd
survivalProbsInterpolator<-function(x, probs, times){
  # Ensure times and probs are sorted by time
  order_idx <- order(times)
  times <- times[order_idx]
  probs <- probs[order_idx]

  # Create an interpolation function
  # Use constant interpolation (f=0) for survival curves (right-continuous)
  # yleft=1 (survival starts at 1), yright=min(probs) (survival at max time)
  f<-stats::approxfun(times, probs, method = "constant", f = 0, yleft = 1, yright = min(probs, na.rm=TRUE), rule = 2)
  sapply(x, function(xi)f(xi))
}

#' @title survprobMatInterpolator
#' @description Interpolate a matrix of survival probabilities for new times.
#' @param probsMat matrix of survival probabilities (rows=observations, cols=times)
#' @param times vector of times corresponding to columns of probsMat
#' @param newtimes vector of new times for interpolation
#' @return matrix of interpolated survival probabilities (rows=newtimes, cols=observations)
#' @noRd
survprobMatInterpolator<-function(probsMat, times, newtimes){
  # Define a helper function to interpolate a single row (one observation's survival curve)
  interpolate1<-function(probs_row){
    # Add time 0 with prob 1 if not present
    if (!0 %in% times) {
        times_aug <- c(0, times)
        probs_row_aug <- c(1, probs_row)
    } else {
        times_aug <- times
        probs_row_aug <- probs_row
    }
    # Interpolate using the single-vector function
    y <- survivalProbsInterpolator(newtimes, probs_row_aug, times_aug)
    y
  }
  # Apply the interpolation function to each row of the probability matrix
  # Note: apply returns matrix with rows=newtimes, cols=observations
  probsMat1<-apply(probsMat, 1, interpolate1)

  # Ensure monotonicity (survival probability should be non-increasing)
  # This step replaces values that increase with the previous minimum value.
  # Apply this correction column-wise (for each observation)
  probsMat2<-apply(probsMat1, 2, function(col_probs){
      # Find the first index where the probability increases compared to the cumulative minimum
      first_increase_idx <- which(col_probs > cummin(col_probs))[1]
      # If an increase is found, replace all values up to that point with the minimum value found so far
      # Also replace any NA values with the minimum
      if (!is.na(first_increase_idx)) {
          replace(col_probs, (seq_along(col_probs) <= first_increase_idx) | is.na(col_probs), min(col_probs[1:first_increase_idx], na.rm = TRUE))
      } else {
          # If no increase, just replace NAs if any exist (e.g., from extrapolation)
          replace(col_probs, is.na(col_probs), min(col_probs, na.rm = TRUE))
      }
  })
  # Ensure the result is a matrix, especially if only one newtime is requested
  if (!is.matrix(probsMat2)) {
      probsMat2 <- matrix(probsMat2, nrow = length(newtimes), ncol = ncol(probsMat1))
  }
  probsMat2
}


#' @title survprobMatListAveraging
#' @description Average a list of survival probability matrices on the cumulative hazard scale.
#' @param listprobsMat list of survival probability matrices (each matrix: rows=newtimes, cols=observations)
#' @return averaged survival probability matrix (rows=newtimes, cols=observations)
#' @importFrom stats na.omit
#' @noRd
survprobMatListAveraging<-function(listprobsMat){
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Check dimensions consistency
  dims <- lapply(listprobsMat, dim)
  if (length(unique(sapply(dims, paste, collapse="x"))) > 1) {
      stop("Matrices in listprobsMat must have the same dimensions.")
  }

  # Create an array to hold cumulative hazards
  HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
  for (i in 1:length(listprobsMat)){
    # Calculate cumulative hazard: -log(S(t))
    # Add small epsilon to avoid log(0)
    HazzardArray[,,i]<--log(listprobsMat[[i]] + 1e-10)
  }
  # Calculate the mean cumulative hazard across models
  MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(na.omit(x)))) # na.omit might be risky if all are NA
  # Convert mean cumulative hazard back to survival probability: exp(-H)
  NewProbs<-exp(-MeanHazzard)
  NewProbs
}
