#' @title Integrator
#' @description Integrate a curve (scores) over a specified time range using the trapezoidal rule.
#' @param times numeric vector of time points.
#' @param scores numeric vector of scores corresponding to time points.
#' @param minmax numeric vector of length 2 specifying the integration limits (min, max).
#' @param scale logical, if TRUE, scale the result by the length of the integration interval.
#' @return numeric value representing the integrated score.
#' @importFrom pracma trapz
#' @noRd
Integrator <- function(times, scores, minmax = c(1, 35), scale = FALSE) {
  # Ensure times and scores have the same length
  if (length(times) != length(scores)) {
    stop("Length of 'times' and 'scores' must be equal.")
  }
  # Ensure minmax is valid
  if (length(minmax) != 2 || minmax[1] >= minmax[2]) {
    stop("'minmax' must be a numeric vector of length 2 with minmax[1] < minmax[2].")
  }

  # Create a mask for times within the specified range
  mask <- (times >= minmax[1]) & (times <= minmax[2])

  # Filter times and scores based on the mask
  timesn <- times[mask]
  scoressn <- scores[mask]

  # Handle cases where no points fall within the range
  if (length(timesn) < 2) {
    if (minmax[1] == minmax[2]) {
      return(0)
    }
    return(0)
  }

  # Sort points by time to ensure correct integration order
  order_idx <- order(timesn)
  timesn <- timesn[order_idx]
  scoressn <- scoressn[order_idx]

  # Calculate the area under the curve using the trapezoidal rule
  AUCsuperlearnMean <- pracma::trapz(timesn, scoressn)

  # Scale the result if requested
  if (scale) {
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

#' @title calculate_expected_time_lost
#' @description Calculate Expected Time Lost (ETL) by integrating event probabilities over time.
#' For survival models: ETL = ∫₀^T (1 - S(t)) dt
#' For competing risks: ETL = ∫₀^T CIF(t) dt
#' @param times numeric vector of time points.
#' @param event_probs numeric vector or matrix of event probabilities (1-S(t) for survival, CIF for CR).
#'   If matrix, each column represents one observation.
#' @param upper_limit numeric upper limit for integration (default: max(times)).
#' @param lower_limit numeric lower limit for integration (default: 0).
#' @return numeric vector of ETL values, one per observation.
#' @noRd
calculate_expected_time_lost <- function(times, event_probs, upper_limit = NULL, lower_limit = 0) {
  # Input validation
  if (!is.numeric(times)) {
    stop("'times' must be numeric")
  }
  if (!is.numeric(event_probs) && !is.matrix(event_probs)) {
    stop("'event_probs' must be numeric vector or matrix")
  }
  if (is.matrix(event_probs) && length(times) != nrow(event_probs)) {
    stop("'times' length must match number of rows in 'event_probs' matrix")
  }
  if (is.vector(event_probs) && length(times) != length(event_probs)) {
    stop("'times' and 'event_probs' must have equal length")
  }

  if (is.null(upper_limit)) {
    upper_limit <- max(times, na.rm = TRUE)
  }

  if (!is.numeric(upper_limit) || length(upper_limit) != 1 || upper_limit <= lower_limit) {
    stop("'upper_limit' must be a single numeric value greater than 'lower_limit'")
  }
  if (!is.numeric(lower_limit) || length(lower_limit) != 1 || lower_limit < 0) {
    stop("'lower_limit' must be a single non-negative numeric value")
  }

  # Ensure times are sorted
  if (!all(diff(times) >= 0)) {
    stop("'times' must be sorted in ascending order")
  }

  # Handle vector input (single observation)
  if (is.vector(event_probs)) {
    event_probs <- matrix(event_probs, nrow = length(times), ncol = 1)
  }

  # Calculate ETL for each observation
  etl_vector <- apply(event_probs, 2, function(obs_probs) {
    Integrator(
      times = times,
      scores = obs_probs,
      minmax = c(lower_limit, upper_limit),
      scale = FALSE
    )
  })

  # Ensure names are preserved for matrix inputs
  if (!is.null(colnames(event_probs))) {
    names(etl_vector) <- colnames(event_probs)
  }

  return(etl_vector)
}
