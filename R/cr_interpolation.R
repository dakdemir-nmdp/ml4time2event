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


#' @title aalenJohansenCIF
#' @description Compute proper CIF using Aalen-Johansen estimator from cause-specific hazards.
#'
#' This function corrects the CIF predictions from cause-specific models by properly
#' accounting for competing risks. Instead of CIF = 1 - S(t), it uses:
#' CIF_j(t) = integral_0^t h_j(u) * S_overall(u) du
#' where S_overall(t) is the overall survival accounting for all causes.
#'
#' @param cause_specific_survs list of survival probability matrices from cause-specific models
#'   Each element should be a matrix with rows=times, cols=observations
#'   The list should have names corresponding to event codes (e.g., "1", "2")
#' @param times vector of time points corresponding to the rows of survival matrices
#' @param event_of_interest character or numeric, the event code of interest
#'
#' @return matrix of calibrated CIF values (rows=times, cols=observations)
#'
#' @details
#' The Aalen-Johansen estimator properly handles competing risks by:
#' 1. Computing cause-specific hazards from each cause-specific survival curve
#' 2. Computing overall survival as the product across all causes
#' 3. Integrating the cause-specific hazard weighted by overall survival
#'
#' This corrects the bias in the naive CIF = 1 - S(t) approach which assumes
#' competing events cannot happen (treating them as censored).
#'
#' @importFrom stats approx
#' @noRd
aalenJohansenCIF <- function(cause_specific_survs, times, event_of_interest) {

  # Input validation
  if (!is.list(cause_specific_survs) || length(cause_specific_survs) == 0) {
    stop("'cause_specific_survs' must be a non-empty list of survival matrices")
  }

  event_of_interest <- as.character(event_of_interest)

  if (!event_of_interest %in% names(cause_specific_survs)) {
    stop("'event_of_interest' (", event_of_interest, ") not found in cause_specific_survs names: ",
         paste(names(cause_specific_survs), collapse = ", "))
  }

  # Get dimensions
  n_times <- length(times)
  n_obs <- ncol(cause_specific_survs[[1]])

  # Validate all matrices have same dimensions
  for (k in names(cause_specific_survs)) {
    if (nrow(cause_specific_survs[[k]]) != n_times) {
      stop("All survival matrices must have ", n_times, " rows (one per time point)")
    }
    if (ncol(cause_specific_survs[[k]]) != n_obs) {
      stop("All survival matrices must have ", n_obs, " columns (one per observation)")
    }
  }

  # Step 1: Calculate cause-specific hazards for each cause
  # hazard_k(t) = - [S_k(t) - S_k(t-1)] / S_k(t-1)
  cause_specific_hazards <- vector("list", length(cause_specific_survs))
  names(cause_specific_hazards) <- names(cause_specific_survs)

  for (k in names(cause_specific_survs)) {
    S_k <- cause_specific_survs[[k]]  # [times, observations]

    # Initialize hazard matrix
    hazard_k <- matrix(0, nrow = n_times, ncol = n_obs)

    # Calculate incremental hazards
    # For t > 1: h(t) = -[S(t) - S(t-1)] / S(t-1)
    for (i in 2:n_times) {
      S_prev <- S_k[i-1, ]
      S_curr <- S_k[i, ]

      # Avoid division by zero
      S_prev[S_prev <= 0] <- 1e-10

      # Calculate hazard increment
      hazard_increment <- -(S_curr - S_prev) / S_prev

      # Ensure non-negative and bounded
      hazard_increment <- pmax(0, pmin(hazard_increment, 1 - 1e-10))

      hazard_k[i, ] <- hazard_increment
    }

    cause_specific_hazards[[k]] <- hazard_k
  }

  # Step 2: Calculate overall survival S_overall(t)
  # S_overall(t) = product over all causes k of [1 - h_k(t)]
  S_overall <- matrix(1, nrow = n_times, ncol = n_obs)

  for (i in 1:n_times) {
    for (k in names(cause_specific_hazards)) {
      # Multiply by (1 - hazard) for each cause
      S_overall[i, ] <- S_overall[i, ] * (1 - cause_specific_hazards[[k]][i, ])
    }
  }

  # Ensure S_overall is non-increasing
  for (j in 1:n_obs) {
    S_overall[, j] <- cummin(S_overall[, j])
  }

  # Step 3: Calculate CIF using Aalen-Johansen formula
  # CIF_j(t) = sum_{s <= t} h_j(s) * S_overall(s-)
  # where S_overall(s-) is S_overall just before time s

  CIF <- matrix(0, nrow = n_times, ncol = n_obs)
  h_j <- cause_specific_hazards[[event_of_interest]]

  # Ensure h_j is a matrix
  if (!is.matrix(h_j)) {
    h_j <- matrix(h_j, ncol = 1)
  }

  # Ensure S_overall is a matrix
  if (!is.matrix(S_overall)) {
    S_overall <- matrix(S_overall, ncol = 1)
  }

  for (i in 1:n_times) {
    if (i == 1) {
      # At first time point, use h_j(t) * 1 (everyone at risk initially)
      if (n_obs == 1) {
        CIF[i, 1] <- h_j[i, 1]
      } else {
        CIF[i, ] <- h_j[i, ]
      }
    } else {
      # CIF(t) = CIF(t-1) + h_j(t) * S_overall(t-1)
      if (n_obs == 1) {
        CIF[i, 1] <- CIF[i-1, 1] + h_j[i, 1] * S_overall[i-1, 1]
      } else {
        CIF[i, ] <- CIF[i-1, ] + h_j[i, ] * S_overall[i-1, ]
      }
    }
  }

  # Ensure CIF is bounded [0, 1] and non-decreasing
  CIF <- pmax(0, pmin(CIF, 1))

  # Ensure it's still a matrix after pmax/pmin operations
  if (!is.matrix(CIF)) {
    CIF <- matrix(CIF, nrow = n_times, ncol = n_obs)
  }

  # Ensure monotonicity
  for (j in 1:n_obs) {
    CIF[, j] <- cummax(CIF[, j])
  }

  return(CIF)
}
