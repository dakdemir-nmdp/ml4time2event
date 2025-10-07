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
  # Use linear interpolation for smooth survival curves (not step functions)
  # yleft=1 (survival starts at 1), yright=min(probs) (survival at max time)
  f<-stats::approxfun(times, probs, method = "linear", yleft = 1, yright = min(probs, na.rm=TRUE), rule = 2)
  sapply(x, function(xi)f(xi))
}

#' @title survprobMatInterpolator
#' @description Interpolate a matrix of survival probabilities for new times.
#' @param probsMat matrix of survival probabilities (rows=times, cols=observations)
#' @param times vector of times corresponding to rows of probsMat
#' @param newtimes vector of new times for interpolation
#' @return matrix of interpolated survival probabilities (rows=newtimes, cols=observations)
#' @noRd
survprobMatInterpolator <- function(probsMat, times, newtimes) {
  # Input: probsMat with rows=times, cols=observations
  # Output: matrix with rows=newtimes, cols=observations

  # Ensure matrix format
  if (!is.matrix(probsMat)) {
    probsMat <- as.matrix(probsMat)
  }

  # Add time 0 if not present
  if (!0 %in% times) {
    times <- c(0, times)
    # Add row of 1s at the beginning (S(0) = 1 for all)
    probsMat <- rbind(rep(1, ncol(probsMat)), probsMat)
  }

  # Interpolate each observation's survival curve (each column)
  n_obs <- ncol(probsMat)
  probs_interp <- matrix(NA, nrow = length(newtimes), ncol = n_obs)

  for (i in seq_len(n_obs)) {
    # Get survival curve for observation i
    surv_curve <- probsMat[, i]
    # Interpolate to new times
    probs_interp[, i] <- survivalProbsInterpolator(newtimes, surv_curve, times)
  }

  # Ensure monotonicity (survival should be non-increasing)
  # Apply to each column (each observation)
  for (i in seq_len(ncol(probs_interp))) {
    probs_interp[, i] <- cummin(probs_interp[, i])
  }

  probs_interp
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
