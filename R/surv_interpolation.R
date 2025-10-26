#' @title survivalProbsInterpolator
#' @description Interpolate survival probabilities for new times using step-function (constant) interpolation.
#'   This is appropriate for survival curves which are typically non-increasing step functions.
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

  # Handle single point case - need to respect yleft=1 for times before the point
  if (length(times) == 1) {
    result <- numeric(length(x))
    for (i in seq_along(x)) {
      if (x[i] < times[1]) {
        result[i] <- 1  # yleft = 1 (survival starts at 1)
      } else {
        result[i] <- probs[1]  # yright = last (and only) prob value
      }
    }
    return(result)
  }

  # Create an interpolation function using constant (step) interpolation
  # yleft=1 (survival starts at 1), yright=last survival probability
  # f=0 means left-continuous step function (use value just before the breakpoint)
  f <- stats::approxfun(times, probs, method = "constant", yleft = 1,
                        yright = probs[length(probs)], rule = 2, f = 0)
  sapply(x, function(xi) f(xi))
}

#' @title survprobMatInterpolator
#' @description Interpolate a matrix of survival probabilities for new times.
#' @param probsMat matrix of survival probabilities (rows=times, cols=observations)
#' @param times vector of times corresponding to rows of probsMat
#' @param new_times vector of new times for interpolation
#' @return matrix of interpolated survival probabilities (rows=new_times, cols=observations)
#' @noRd
survprobMatInterpolator <- function(probsMat, times, new_times) {
  # Input: probsMat with rows=times, cols=observations
  # Output: matrix with rows=new_times, cols=observations

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
  # Ensure new_times is sorted
  new_times_order <- order(new_times)
  new_times_sorted <- new_times[new_times_order]
  
  probs_interp_sorted <- matrix(NA, nrow = length(new_times_sorted), ncol = n_obs)

  for (i in seq_len(n_obs)) {
    # Get survival curve for observation i
    surv_curve <- probsMat[, i]
    # Interpolate to new times
    probs_interp_sorted[, i] <- survivalProbsInterpolator(new_times_sorted, surv_curve, times)
  }

  # Ensure monotonicity (survival probabilities should be non-increasing)
  # This is a safety measure in case the input data has any non-monotonic artifacts
  for (i in seq_len(ncol(probs_interp_sorted))) {
    probs_interp_sorted[, i] <- cummin(probs_interp_sorted[, i])
  }

  # Revert to original order of new_times
  probs_interp <- probs_interp_sorted[order(new_times_order), , drop = FALSE]


  probs_interp
}


#' @title survprobMatListAveraging
#' @description Average a list of survival probability matrices on the cumulative hazard scale.
#' @param listprobsMat list of survival probability matrices (each matrix: rows=new_times, cols=observations)
#' @param na.rm logical, whether to remove NAs when averaging.
#' @return averaged survival probability matrix (rows=new_times, cols=observations)
#' @importFrom stats na.omit
#' @noRd
survprobMatListAveraging<-function(listprobsMat, na.rm = FALSE){
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Filter out NULL entries
  listprobsMat <- Filter(Negate(is.null), listprobsMat)

  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Check dimensions consistency
  dims <- lapply(listprobsMat, dim)
  
  # Handle cases where some predictions might be vectors not matrices
  is_matrix <- sapply(dims, function(d) !is.null(d) && length(d) == 2)
  if (!all(is_matrix)) {
      warning("Some predictions are not matrices and will be excluded.")
      listprobsMat <- listprobsMat[is_matrix]
      dims <- dims[is_matrix]
      if (length(listprobsMat) <= 1) return(if(length(listprobsMat) == 1) listprobsMat[[1]] else NULL)
  }
  
  dim_strings <- sapply(dims, paste, collapse="x")

  # If dimensions are inconsistent, throw error
  if (length(unique(dim_strings)) > 1) {
    stop("All matrices in listprobsMat must have the same dimensions. Found dimensions: ",
         paste(unique(dim_strings), collapse = ", "))
  }

  # Create an array to hold cumulative hazards
  HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
  for (i in seq_along(listprobsMat)){
    # Calculate cumulative hazard: -log(S(t))
    # Add small epsilon to avoid log(0)
    HazzardArray[,,i]<--log(listprobsMat[[i]] + 1e-10)
  }
  # Calculate the mean cumulative hazard across models
  MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(x, na.rm = na.rm)))
  # Convert mean cumulative hazard back to survival probability: exp(-H)
  NewProbs<-exp(-MeanHazzard)

  # Ensure probabilities are bounded [0, 1]
  NewProbs <- pmax(0, pmin(NewProbs, 1))
  
  # Preserve matrix dimensions (apply can drop dims with single row/col)
  if (!is.matrix(NewProbs)) {
    dim(NewProbs) <- dim(listprobsMat[[1]])
  }

  NewProbs
}
