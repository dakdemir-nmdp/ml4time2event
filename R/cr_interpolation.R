#' @title cifInterpolator
#' @description Interpolate Cumulative Incidence Function (CIF) for new times using linear interpolation.
#'   Properly propagates NAs: if interpolating between points where any endpoint is NA, result is NA.
#' @param x new times for interpolation
#' @param probs vector of CIF probabilities (may contain NAs)
#' @param times vector of times corresponding to probs
#' @return interpolated CIF values at times x
#' @importFrom stats approx
#' @noRd
cifInterpolator<-function(x, probs, times){
  # Sort times and probs together to handle unsorted input
  if (length(times) > 1) {
    sort_order <- order(times)
    times <- times[sort_order]
    probs <- probs[sort_order]
  }

  # Handle single time/prob case by adding (0, 0) point if needed
  if (length(times) == 1 && length(probs) == 1) {
    if (times[1] == 0) {
      # If the single point is at time 0, extrapolate with constant value
      return(rep(probs[1], length(x)))
    } else {
      # Add (0, 0) point for proper linear interpolation
      times <- c(0, times)
      probs <- c(0, probs)
    }
  }

  # Use stats::approx with ties = "ordered" and rule = 2 for proper handling
  # This will interpolate but stats::approx doesn't propagate NAs correctly
  # We need to post-process to handle NA propagation
  result <- stats::approx(x = times, y = probs, xout = x,
                          method = "linear", yleft = 0,
                          yright = probs[length(probs)],
                          rule = 2, ties = "ordered")$y

  # Post-process to propagate NAs: if interpolating between any NA endpoint, set to NA
  for (i in seq_along(x)) {
    xi <- x[i]
    # Find the two surrounding time points
    if (xi <= times[1]) {
      # Before first point - already handled by yleft
      next
    } else if (xi >= times[length(times)]) {
      # After last point - check if last point is NA
      if (is.na(probs[length(probs)])) {
        result[i] <- NA
      }
    } else {
      # Between points - find surrounding indices
      idx_after <- which(times > xi)[1]
      idx_before <- idx_after - 1
      # If either endpoint is NA, result should be NA
      if (is.na(probs[idx_before]) || is.na(probs[idx_after])) {
        result[i] <- NA
      }
    }
  }

  result
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
    # Interpolate using the single-vector function (now uses linear interpolation)
    cifInterpolator(newtimes, probs_row_aug, times_aug)
  }
  # Apply the interpolation function to each row of the probability matrix
  # Note: apply returns matrix with rows=newtimes, cols=observations
  probsMat1<-apply(probsMat, 1, interpolate1)

  # Ensure probsMat1 is always a matrix (when newtimes has length 1, apply returns a vector)
  if (!is.matrix(probsMat1)) {
    probsMat1 <- matrix(probsMat1, nrow = length(newtimes), ncol = nrow(probsMat))
  }

  # Enforce monotonicity (CIF should be non-decreasing)
  # This is crucial for correctness, even if models should ideally produce this.
  for (i in 1:ncol(probsMat1)) {
    probsMat1[, i] <- cummax(probsMat1[, i])
  }

  probsMat1
}



#' @title cifMatListAveraging
#' @description Average a list of CIF matrices, either on probability or cumulative hazard scale.
#' @param listprobsMat list of CIF matrices (each matrix: rows=newtimes, cols=observations)
#' @param type character, either "CumHaz" (average on cumulative hazard scale) or "prob" (average on probability scale)
#' @param na.rm logical, whether to remove NAs when averaging.
#' @return averaged CIF matrix (rows=newtimes, cols=observations)
#' @noRd
cifMatListAveraging<-function(listprobsMat, type="CumHaz", na.rm = FALSE){
  # Validate type parameter first, regardless of list length
  if (!type %in% c("CumHaz", "prob")) {
    stop("Type must be either 'CumHaz' or 'prob'")
  }

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

  if (type=="CumHaz"){
    # Create an array to hold cumulative hazards
    HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
    for (i in 1:length(listprobsMat)){
      # Calculate cumulative hazard: -log(1 - P)
      # Add small epsilon to avoid log(0)
      HazzardArray[,,i]<--log(1 - listprobsMat[[i]] + 1e-10)
    }
    # Calculate the mean cumulative hazard across models
    MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(x, na.rm = na.rm)))
    # Convert mean cumulative hazard back to probability: 1 - exp(-H)
    NewProbs<-1-exp(-MeanHazzard)
  } else if (type=="prob"){
    # Create an array to hold probabilities
    ProbsArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
    for (i in 1:length(listprobsMat)){
      ProbsArray[,,i]<-listprobsMat[[i]]
    }
    # Calculate the mean probability across models
    NewProbs<-apply(ProbsArray, c(1,2),function(x)(mean(x, na.rm = na.rm)))
    # Ensure probabilities are bounded between 0 and 1
    NewProbs <- pmax(pmin(NewProbs, 1.0), 0.0)
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


#' @title aalenJohansenFromCoxModels
#' @description Calculate CIF using proper Aalen-Johansen estimator from cause-specific Cox models.
#' @param cox_models list of fitted Cox models for each cause
#' @param newdata data frame of new observations
#' @param times vector of time points
#' @param event_of_interest the event type to calculate CIF for
#' @return matrix of CIF values (rows=times, cols=observations)
#' @noRd
aalenJohansenFromCoxModels <- function(cox_models, newdata, times, event_of_interest) {
  
  event_of_interest <- as.character(event_of_interest)
  n_times <- length(times)
  n_obs <- nrow(newdata)
  
  # Initialize CIF matrix
  cif_matrix <- matrix(0, nrow = n_times, ncol = n_obs)
  

  
  # For each observation
  for (i in 1:n_obs) {
    obs_data <- newdata[i, , drop = FALSE]

    if (anyNA(obs_data)) {
      cif_matrix[, i] <- NA_real_
      next
    }
    
    # Calculate cause-specific hazards at each time point
    cause_hazards <- vector("list", length(cox_models))
    names(cause_hazards) <- names(cox_models)
    
    for (cause in names(cox_models)) {
      cox_model <- cox_models[[cause]]
      
      # Get baseline cumulative hazard. This will now error out if newdata is problematic.
      base_surv <- survival::survfit(cox_model, newdata = obs_data)
      
      # Extract hazard increments (Nelson-Aalen style)
      if (length(base_surv$time) > 0 && !any(is.na(base_surv$surv))) {
        # Interpolate baseline cumulative hazard to our time grid
        cum_base_haz <- stats::approx(
          x = c(0, base_surv$time), 
          y = c(0, -log(pmax(base_surv$surv, 1e-10))), # Avoid log(0)
          xout = times, 
          method = "constant", 
          f = 0, 
          rule = 2
        )$y
        
        # Convert to hazard increments
        haz_increments <- diff(c(0, cum_base_haz))
      } else {
        haz_increments <- rep(0, n_times)
      }
      
      cause_hazards[[cause]] <- haz_increments
    }
    
    # Calculate overall survival using all cause-specific hazards
    overall_surv <- rep(1, n_times)
    cum_overall_haz <- 0
    
    for (t in 1:n_times) {
      # Add hazard increment from all causes
      for (cause in names(cause_hazards)) {
        cum_overall_haz <- cum_overall_haz + cause_hazards[[cause]][t]
      }
      overall_surv[t] <- exp(-cum_overall_haz)
    }
    

    
    # Calculate CIF using Aalen-Johansen formula
    # CIF_j(t) = ∫_0^t S(s-) * λ_j(s) ds
    # Discrete version: CIF_j(t) = Σ_{s≤t} S(s-) * Δλ_j(s)
    
    if (event_of_interest %in% names(cause_hazards)) {
      target_hazards <- cause_hazards[[event_of_interest]]
      
      for (t in 1:n_times) {
        if (t == 1) {
          # At first time point, S(0-) = 1
          cif_matrix[t, i] <- 1.0 * target_hazards[t]
        } else {
          # CIF(t) = CIF(t-1) + S(t-1) * Δλ_j(t)
          cif_matrix[t, i] <- cif_matrix[t-1, i] + overall_surv[t-1] * target_hazards[t]
        }
      }
    }
  }
  
  # Ensure CIF is bounded [0, 1] and monotonic
  # Use matrix() to preserve dimensions after pmax/pmin operations
  original_dims <- dim(cif_matrix)
  cif_matrix <- matrix(pmax(0, pmin(as.vector(cif_matrix), 1)), 
                       nrow = original_dims[1], ncol = original_dims[2])
  
  # Ensure monotonicity
  if (n_obs == 1) {
    cif_matrix[, 1] <- cummax(cif_matrix[, 1])
  } else {
    for (j in 1:n_obs) {
      cif_matrix[, j] <- cummax(cif_matrix[, j])
    }
  }
  
  return(cif_matrix)
}
