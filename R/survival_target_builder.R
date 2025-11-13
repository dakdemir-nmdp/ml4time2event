# ==============================================================================
# Survival Target Construction for Super Learner
# ==============================================================================

buildObservedSurvivalMatrix <- function(data, timevar, eventvar, eval_times) {
  if (missing(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (missing(timevar) || !timevar %in% colnames(data)) {
    stop("'timevar' must be a column in data")
  }
  if (missing(eventvar) || !eventvar %in% colnames(data)) {
    stop("'eventvar' must be a column in data")
  }
  if (missing(eval_times) || !is.numeric(eval_times)) {
    stop("'eval_times' must be numeric")
  }
  
  # Ensure eval_times includes 0 and is sorted
  eval_times <- sort(unique(c(0, eval_times)))
  
  n_obs <- nrow(data)
  n_times <- length(eval_times)
  
  # Initialize observed survival matrix
  obs_surv <- matrix(NA, nrow = n_times, ncol = n_obs)
  
  # For each observation, determine survival status at each eval time
  for (i in seq_len(n_obs)) {
    obs_time <- data[[timevar]][i]
    obs_event <- data[[eventvar]][i]
    
    for (t_idx in seq_len(n_times)) {
      eval_t <- eval_times[t_idx]
      
      if (eval_t == 0) {
        # Everyone starts alive
        obs_surv[t_idx, i] <- 1
      } else if (eval_t <= obs_time) {
        # Still alive/uncensored at this time
        obs_surv[t_idx, i] <- 1
      } else if (obs_event == 1) {
        # Event occurred before eval_t
        obs_surv[t_idx, i] <- 0
      } else {
        # Censored before eval_t - use Kaplan-Meier estimate
        # For individual observations that are censored, we can't know the exact survival
        # We'll use the population KM estimate at this time point
        km_fit <- tryCatch(
          survival::survfit(survival::Surv(data[[timevar]], data[[eventvar]]) ~ 1),
          error = function(e) NULL
        )
        
        if (!is.null(km_fit)) {
          # Interpolate KM estimate at eval_t
          km_times <- c(0, km_fit$time)
          km_surv <- c(1, km_fit$surv)
          
          # Find survival probability at eval_t
          if (eval_t <= max(km_times)) {
            surv_at_t <- survivalProbsInterpolator(eval_t, km_surv, km_times)
            obs_surv[t_idx, i] <- surv_at_t
          } else {
            # Beyond follow-up, use last available estimate
            obs_surv[t_idx, i] <- tail(km_surv, 1)
          }
        } else {
          # Fallback if KM fitting fails
          obs_surv[t_idx, i] <- 0.5
        }
      }
    }
  }
  
  # Ensure monotonicity for each individual (survival should be non-increasing)
  for (i in seq_len(n_obs)) {
    obs_surv[, i] <- cummin(obs_surv[, i])
  }
  
  # Bound probabilities to [0, 1]
  obs_surv <- pmax(0, pmin(obs_surv, 1))
  
  # Ensure return is a matrix
  dim(obs_surv) <- c(n_times, n_obs)
  
  obs_surv
}

buildObservedCIFMatrix <- function(data, timevar, eventvar, cause_of_interest, eval_times) {
  if (missing(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (missing(timevar) || !timevar %in% colnames(data)) {
    stop("'timevar' must be a column in data")
  }
  if (missing(eventvar) || !eventvar %in% colnames(data)) {
    stop("'eventvar' must be a column in data")
  }
  if (missing(cause_of_interest) || !is.numeric(cause_of_interest)) {
    stop("'cause_of_interest' must be numeric")
  }
  if (missing(eval_times) || !is.numeric(eval_times)) {
    stop("'eval_times' must be numeric")
  }
  
  # Ensure eval_times includes 0 and is sorted
  eval_times <- sort(unique(c(0, eval_times)))
  
  n_obs <- nrow(data)
  n_times <- length(eval_times)
  
  # Initialize observed CIF matrix
  obs_cif <- matrix(NA, nrow = n_times, ncol = n_obs)
  
  # Build population-level CIF estimate for reference
  if (requireNamespace("cmprsk", quietly = TRUE)) {
    cif_fit <- tryCatch({
      cmprsk::cuminc(data[[timevar]], data[[eventvar]], cencode = 0)
    }, error = function(e) NULL)
  } else {
    cif_fit <- NULL
  }
  
  # For each observation, determine CIF status at each eval time
  for (i in seq_len(n_obs)) {
    obs_time <- data[[timevar]][i]
    obs_event <- data[[eventvar]][i]
    
    for (t_idx in seq_len(n_times)) {
      eval_t <- eval_times[t_idx]
      
      if (eval_t == 0) {
        # No events at time 0
        obs_cif[t_idx, i] <- 0
      } else if (eval_t <= obs_time) {
        if (obs_event == cause_of_interest) {
          # Event of interest occurred exactly at obs_time
          if (eval_t == obs_time) {
            obs_cif[t_idx, i] <- 1
          } else {
            obs_cif[t_idx, i] <- 0
          }
        } else {
          # Still at risk or competing event hasn't occurred yet
          obs_cif[t_idx, i] <- 0
        }
      } else if (obs_event == cause_of_interest) {
        # Event of interest occurred before eval_t
        obs_cif[t_idx, i] <- 1
      } else if (obs_event > 0 && obs_event != cause_of_interest) {
        # Competing event occurred - can never experience cause of interest
        obs_cif[t_idx, i] <- 0
      } else {
        # Censored before eval_t - use population CIF estimate
        if (!is.null(cif_fit)) {
          cause_key <- paste(cause_of_interest, cause_of_interest, sep = " ")
          if (cause_key %in% names(cif_fit)) {
            cif_data <- cif_fit[[cause_key]]
            cif_times <- c(0, cif_data$time)
            cif_probs <- c(0, cif_data$est)
            
            # Interpolate CIF at eval_t
            if (eval_t <= max(cif_times)) {
              cif_at_t <- survivalProbsInterpolator(eval_t, cif_probs, cif_times)
              obs_cif[t_idx, i] <- cif_at_t
            } else {
              # Beyond follow-up, use last available estimate
              obs_cif[t_idx, i] <- tail(cif_probs, 1)
            }
          } else {
            # Fallback if cause not found
            obs_cif[t_idx, i] <- 0.1
          }
        } else {
          # Fallback if CIF fitting fails
          obs_cif[t_idx, i] <- 0.1
        }
      }
    }
  }
  
  # Ensure monotonicity for each individual (CIF should be non-decreasing)
  for (i in seq_len(n_obs)) {
    obs_cif[, i] <- cummax(obs_cif[, i])
  }
  
  # Bound probabilities to [0, 1]
  obs_cif <- pmax(0, pmin(obs_cif, 1))
  
  # Ensure return is a matrix
  dim(obs_cif) <- c(n_times, n_obs)
  
  obs_cif
}