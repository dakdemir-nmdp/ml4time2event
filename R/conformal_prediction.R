#' Fit Censoring Model for Conformal Prediction
#'
#' Fits a Kaplan-Meier estimator for the censoring distribution.
#'
#' @param data A data frame.
#' @param time_col Character, name of time column.
#' @param event_col Character, name of event column.
#'
#' @return A `survfit` object.
#' @keywords internal
ml4t2e_fit_censoring <- function(data, time_col, event_col) {
  # Censoring event: status == 0
  # We assume standard coding: 0 = censored, 1+ = event.
  # So for censoring distribution, 'event' is status == 0.
  
  # Using KM for simplicity and robustness.
  # Note: This assumes censoring is independent of covariates (random censoring).
  # If censoring depends on X, a Cox model would be better.
  # Given "computationally efficient and simple", KM is the right starting point.
  
  f <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, " == 0) ~ 1"))
  survival::survfit(f, data = data)
}

#' Predict Censoring Probabilities
#'
#' @param model A `survfit` object (censoring model).
#' @param times Numeric vector of times to predict.
#'
#' @return Numeric vector of probabilities P(C > t).
#' @keywords internal
ml4t2e_predict_censoring <- function(model, times) {
  # summary(model, times = times)$surv gives P(C > t)
  # We need to handle times beyond the max follow-up carefully.
  # survival::summary.survfit extends to the requested times.
  
  # Note: summary.survfit might fail if times are outside range or if times is empty.
  if (length(times) == 0) return(numeric(0))
  
  # Use stepfun for robust interpolation/extrapolation (constant after last event)
  # survfit object has $time and $surv (and $n.risk etc)
  # We want P(C > t).
  # $surv is P(C > t) at $time.
  # Initial survival is 1 at t=0.
  
  sf <- stats::stepfun(model$time, c(1, model$surv))
  sf(times)
}

#' Compute Conformal Scores
#'
#' @param pred_matrix Matrix of predicted survival/CIF probabilities (rows=obs, cols=times).
#' @param data Calibration data frame.
#' @param time_grid Numeric vector of times corresponding to columns of pred_matrix.
#' @param censoring_model Fitted censoring model.
#' @param task The task object (to get column names).
#' @param event_of_interest For CR, the cause of interest. NULL for survival.
#'
#' @return A list with `scores` (matrix) and `weights` (matrix).
#' @keywords internal
.compute_scores_core <- function(pred_matrix, data, time_grid, censoring_model, task, event_of_interest = NULL) {
  n <- nrow(data)
  m <- length(time_grid)
  
  time_col <- task$time_col
  event_col <- task$event_col
  
  T_i <- data[[time_col]]
  Delta_i <- data[[event_col]] # 0=cens, 1..k=event
  
  # Predict G(t) for all t in time_grid
  G_t <- ml4t2e_predict_censoring(censoring_model, time_grid)
  
  # Predict G(T_i) for all i
  G_Yi <- ml4t2e_predict_censoring(censoring_model, T_i)
  
  # Avoid division by zero (small epsilon)
  G_Yi <- pmax(G_Yi, 1e-5)
  G_t <- pmax(G_t, 1e-5)
  
  scores_mat <- matrix(NA, nrow=n, ncol=m)
  weights_mat <- matrix(0, nrow=n, ncol=m)
  
  is_cr <- !is.null(event_of_interest)
  
  for (j in 1:m) {
    t_val <- time_grid[j]
    G_val <- G_t[j]
    
    # Indices
    idx_gt <- T_i > t_val
    idx_le <- T_i <= t_val
    
    # Weights
    # For Y > t: w = 1/G(t)
    weights_mat[idx_gt, j] <- 1 / G_val
    
    # For Y <= t, Event: w = 1/G(Y_i)
    # Note: We only weight observed events. Censored before t are 0 weight.
    # For Survival: event is Delta == 1 (or >0 if simple survival)
    # For CR: event is Delta != 0
    
    idx_event_observed <- if (is_cr) (Delta_i != 0) else (Delta_i == 1)
    
    idx_obs_event <- idx_le & idx_event_observed
    weights_mat[idx_obs_event, j] <- 1 / G_Yi[idx_obs_event]
    
    # Targets
    if (!is_cr) {
      # Survival
      # Target is 1(T > t)
      # Y > t -> 1
      # Y <= t, Event -> 0
      target <- rep(NA, n)
      target[idx_gt] <- 1
      target[idx_obs_event] <- 0
      
      # Score s = |Target - Pred|
      valid <- weights_mat[, j] > 0
      scores_mat[valid, j] <- abs(target[valid] - pred_matrix[valid, j])
      
    } else {
      # CIF for cause k
      # Target is 1(T <= t, Cause = k)
      # Y > t -> 0
      # Y <= t, Cause=k -> 1
      # Y <= t, Cause!=k, Event -> 0
      
      target <- rep(NA, n)
      target[idx_gt] <- 0
      
      idx_cause <- idx_le & (Delta_i == event_of_interest)
      target[idx_cause] <- 1
      
      idx_other <- idx_le & (Delta_i != 0) & (Delta_i != event_of_interest)
      target[idx_other] <- 0
      
      valid <- weights_mat[, j] > 0
      scores_mat[valid, j] <- abs(target[valid] - pred_matrix[valid, j])
    }
  }
  
  list(scores = scores_mat, weights = weights_mat)
}

#' Weighted Quantile
#'
#' Computes the weighted quantile of a distribution.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of weights.
#' @param alpha Significance level (e.g. 0.1 for 90% coverage).
#'
#' @return The weighted (1-alpha) quantile.
#' @keywords internal
ml4t2e_weighted_quantile <- function(x, w, alpha) {
  # Filter NAs and zero weights
  valid <- !is.na(x) & !is.na(w) & w > 0
  x <- x[valid]
  w <- w[valid]
  
  if (length(x) == 0) return(NA_real_)
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w)
  sw <- cw[length(cw)]
  
  # Target cumulative weight
  # We want P(Score <= q) >= 1 - alpha
  target <- (1 - alpha) * sw
  
  # Find first index where cw >= target
  idx <- which(cw >= target)[1]
  
  if (is.na(idx)) return(max(x)) # Should not happen if sw > 0
  
  x[idx]
}

#' Reshape Predictions to Matrix
#'
#' Helper to convert tidy predictions to matrix for score computation.
#'
#' @param preds Tidy prediction data frame.
#' @param data Data frame (to align rows).
#' @param time_grid Numeric vector of times.
#' @param type "surv" or "cif".
#'
#' @return Matrix (n x m).
#' @keywords internal
.reshape_preds_to_matrix <- function(preds, data, time_grid, type, id_col = NULL) {
  # Filter preds to only times in time_grid
  preds_sub <- preds[preds$time %in% time_grid, ]
  
  # Pivot using tidyr
  if (type == "surv") {
    wide <- preds_sub %>%
      dplyr::select(id, time, surv) %>%
      tidyr::pivot_wider(names_from = time, values_from = surv)
  } else {
    wide <- preds_sub %>%
      dplyr::select(id, time, cif) %>%
      tidyr::pivot_wider(names_from = time, values_from = cif)
  }
  
  # Align rows with data
  if (!is.null(id_col) && id_col %in% colnames(data)) {
    # Match wide$id with data[[id_col]]
    # Ensure we match correctly even if types differ slightly (e.g. factor vs char)
    target_ids <- as.character(data[[id_col]])
    wide_ids <- as.character(wide$id)
    
    # Check if we have all IDs
    if (!all(target_ids %in% wide_ids)) {
      rlang::warn("Some calibration data IDs are missing from predictions.")
    }
    
    idx <- match(target_ids, wide_ids)
    wide <- wide[idx, ]
  } else {
    # If no id_col provided or not in data, we assume predictions are already aligned 
    # or we can't safely reorder. 
    # However, pivot_wider might change order.
    # We should warn if we can't align.
    if (nrow(wide) != nrow(data)) {
      rlang::warn("Prediction matrix rows do not match data rows and no ID column provided for alignment.")
    }
  }
  
  # Select time columns in order
  # wide columns (except id) are times.
  # We need to ensure we select them in the order of time_grid.
  # pivot_wider names are strings.
  col_names <- as.character(time_grid)
  
  # Check if all time columns exist
  missing_times <- setdiff(col_names, colnames(wide))
  if (length(missing_times) > 0) {
    rlang::abort("Some time points are missing from pivoted predictions.")
  }
  
  mat <- as.matrix(wide[, col_names])
  mat
}
