#' @title cifInterpolator
#' @description Interpolate Cumulative Incidence Function (CIF) for new times.
#' @param x new times for interpolation
#' @param probs vector of CIF probabilities
#' @param times vector of times corresponding to probs
#' @return interpolated CIF values at times x
#' @importFrom stats approxfun
#' @noRd
cifInterpolator<-function(x, probs, times){
  # Create an interpolation function based on existing times and probabilities
  f<-stats::approxfun(times, probs, method = "linear", yleft = 0, yright = max(probs, na.rm = TRUE), rule = 2) # Ensure rule=2 to extrapolate using nearest value
  # Apply the interpolation function to the new times
  sapply(x, function(xi)f(xi))
}

#' @title cifMatInterpolaltor
#' @description Interpolate a matrix of CIFs for new times.
#' @param probsMat matrix of CIFs (rows=observations, cols=times)
#' @param times vector of times corresponding to columns of probsMat
#' @param newtimes vector of new times for interpolation
#' @return matrix of interpolated CIFs (rows=newtimes, cols=observations)
#' @noRd
cifMatInterpolaltor<-function(probsMat, times,newtimes){
  # Define a helper function to interpolate a single row (one observation's CIF curve)
  interpolate1<-function(probs_row){
    # Add time 0 with CIF 0 if not present
    if (!0 %in% times) {
        times_aug <- c(0, times)
        probs_row_aug <- c(0, probs_row)
    } else {
        times_aug <- times
        probs_row_aug <- probs_row
    }
    # Interpolate using the single-vector function
    cifInterpolator(newtimes, probs_row_aug, times_aug)
  }
  # Apply the interpolation function to each row of the probability matrix
  # Note: apply returns matrix with rows=newtimes, cols=observations
  probsMat1<-apply(probsMat, 1, interpolate1)

  # Ensure monotonicity (CIF should be non-decreasing)
  # This step replaces values that decrease with the previous maximum value.
  # Apply this correction column-wise (for each observation)
  probsMat2<-apply(probsMat1, 2, function(col_probs){
      # Find the first index where the probability decreases compared to the cumulative maximum
      first_decrease_idx <- which(col_probs < cummax(col_probs))[1]
      # If a decrease is found, replace all values up to that point with the maximum value found so far
      # Also replace any NA values with the maximum
      if (!is.na(first_decrease_idx)) {
          replace(col_probs, (seq_along(col_probs) <= first_decrease_idx) | is.na(col_probs), max(col_probs[1:first_decrease_idx], na.rm = TRUE))
      } else {
          # If no decrease, just replace NAs if any exist (e.g., from extrapolation)
          replace(col_probs, is.na(col_probs), max(col_probs, na.rm = TRUE))
      }
  })
  # Ensure the result is a matrix, especially if only one newtime is requested
  if (!is.matrix(probsMat2)) {
      probsMat2 <- matrix(probsMat2, nrow = length(newtimes), ncol = ncol(probsMat1))
  }
  probsMat2
}



#' @title cifMatListAveraging
#' @description Average a list of CIF matrices, either on probability or cumulative hazard scale.
#' @param listprobsMat list of CIF matrices (each matrix: rows=newtimes, cols=observations)
#' @param type character, either "CumHaz" (average on cumulative hazard scale) or "prob" (average on probability scale)
#' @return averaged CIF matrix (rows=newtimes, cols=observations)
#' @noRd
cifMatListAveraging<-function(listprobsMat, type="CumHaz"){
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Check dimensions consistency
  dims <- lapply(listprobsMat, dim)
  if (length(unique(sapply(dims, paste, collapse="x"))) > 1) {
      stop("Matrices in listprobsMat must have the same dimensions.")
  }

  if (type=="CumHaz"){
    # Create an array to hold cumulative hazards
    HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
    for (i in 1:length(listprobsMat)){
      # Calculate cumulative hazard: -log(1 - P)
      # Add small epsilon to avoid log(0)
      HazzardArray[,,i]<--log(1 - listprobsMat[[i]] + 1e-10)
    }
    # Calculate the mean cumulative hazard across models
    MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(x, na.rm = TRUE))) # Use na.rm=TRUE
    # Convert mean cumulative hazard back to probability: 1 - exp(-H)
    NewProbs<-1-exp(-MeanHazzard)
  } else if (type=="prob"){
    # Create an array to hold probabilities
    ProbsArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
    for (i in 1:length(listprobsMat)){
      ProbsArray[,,i]<-listprobsMat[[i]]
    }
    # Calculate the mean probability across models
    NewProbs<-apply(ProbsArray, c(1,2),function(x)(mean(x, na.rm = TRUE))) # Use na.rm=TRUE
    # Ensure probabilities are bounded between 0 and 1
    NewProbs <- pmax(pmin(NewProbs, 1.0), 0.0)
  } else {
      stop("Type must be either 'CumHaz' or 'prob'")
  }
  NewProbs
}
