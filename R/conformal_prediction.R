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
  if (length(times) == 0) {
    return(numeric(0))
  }

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

  scores_mat <- matrix(NA, nrow = n, ncol = m)
  weights_mat <- matrix(0, nrow = n, ncol = m)

  is_cr <- !is.null(event_of_interest)

  # Prepare logical matrices (n x m)
  # Rows = observations, Cols = time points
  # Using outer() creates n x m logical matrix
  mask_gt <- outer(T_i, time_grid, ">")
  mask_le <- !mask_gt

  # Weights computation
  # 1. Case T_i > t: weight = 1 / G(t)
  # Replicate G_t (m) to match rows (n) -> repeats G_t n times row-wise? No.
  # matrix(rep(G_t, each=n), nrow=n) makes columns identical to G_t values.
  G_t_mat <- matrix(rep(G_t, each = n), nrow = n, ncol = m)
  trm1 <- mask_gt / G_t_mat
  # Where mask_gt is FALSE, trm1 is 0.

  # 2. Case T_i <= t AND Observed Event: weight = 1 / G(T_i)
  G_Yi_mat <- matrix(rep(G_Yi, m), nrow = n, ncol = m)

  if (is_cr) {
    event_observed <- (Delta_i != 0)
  } else {
    event_observed <- (Delta_i == 1)
  }

  # Logical AND with vector recycling for event_observed (length n)
  mask_obs_event <- mask_le & event_observed

  trm2 <- matrix(0, nrow = n, ncol = m)
  trm2[mask_obs_event] <- (1 / G_Yi_mat)[mask_obs_event]

  weights_mat <- trm1 + trm2

  # Targets computation
  if (!is_cr) {
    # Survival
    # Target 1 if T > t
    # Target 0 if T <= t AND Event
    # (Ignored if censored)

    # Start with matrix of 0s or NAs?
    # Since we only use valid weights, we can default to 0.
    target_mat <- matrix(0, nrow = n, ncol = m)

    # Set 1s where T > t
    target_mat[mask_gt] <- 1

    # 0s already set.
  } else {
    # Competing Risks
    # Target 1 if T <= t AND Cause == k
    # Target 0 if T > t
    # Target 0 if T <= t AND Cause != k AND Event

    target_mat <- matrix(0, nrow = n, ncol = m)

    # Set 1s where T <= t AND Cause == k
    idx_cause_k <- (Delta_i == event_of_interest)
    mask_cause_k <- mask_le & idx_cause_k

    target_mat[mask_cause_k] <- 1

    # 0s already set.
  }

  # Scores computation
  # s = |Target - Pred| where weight > 0
  scores_mat <- abs(target_mat - pred_matrix)

  # Set scores to NA where weight is 0 (censored before t) to be safe/clean
  # logical matrix of valid weights
  mask_valid <- weights_mat > 0

  # Assign NA to invalid scores to signify they shouldn't be used (though weighting handles it usually)
  # But the ml4t2e_weighted_quantile function handles explicit weights.
  # For consistency with previous loop which set scores only for valid cases:
  scores_mat[!mask_valid] <- NA

  list(scores = scores_mat, weights = weights_mat)
}

#' Weighted Quantile
#'
#' Computes the weighted quantile of a distribution.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of weights.
#' @param alpha Significance level (e.g. 0.1 for 90 percent coverage).
#'
#' @return Numeric value representing the weighted (1-alpha) quantile.
#' @export
ml4t2e_weighted_quantile <- function(x, w, alpha) {
  # Filter NAs and zero weights
  valid <- !is.na(x) & !is.na(w) & w > 0
  x <- x[valid]
  w <- w[valid]

  if (length(x) == 0) {
    return(NA_real_)
  }

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

  if (is.na(idx)) {
    return(max(x))
  } # Should not happen if sw > 0

  x[idx]
}



#' Calibrate a fitted model for conformal prediction
#'
#' Computes conformal scores on a calibration set and attaches them to the fit object.
#'
#' @param object A `t2e_fit` object.
#' @param data A data frame to be used for calibration. Must contain the same time/event columns as the training data.
#' @param split Optional numeric (0 < split < 1). If provided, and `data` is NULL, the original training data (stored in `object$task$data`)
#'   will be split, with `split` fraction used for calibration. *Note: this requires re-fitting if the model was already trained on all data,
#'   which this function does NOT do. This argument is reserved for future pipelines where fitting and calibration are coupled.*
#'   Currently, you must provide `data`.
#'
#' @return A `t2e_fit` object with a `$conformal` slot containing scores and weights.
#' @export
ml4t2e_calibrate <- function(object, data) {
  if (!inherits(object, "t2e_fit")) {
    rlang::abort("`object` must be a `t2e_fit`.")
  }

  task <- object$task
  # Ensure data has necessary columns
  required_cols <- c(task$time_col, task$event_col, task$features)
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    rlang::abort(paste0("Calibration `data` is missing columns: ", paste(missing, collapse = ", ")))
  }

  # Ensure we have a valid ID col if it exists in task
  if (!is.null(task$id_col) && !(task$id_col %in% colnames(data))) {
    # If ID col is missing but used in task, we might need it for matching.
    # But if we treat data as new rows, we might generate IDs?
    # For now, require it if it was in train.
    rlang::warn("ID column from task not found in calibration data. Generated IDs may simplify matching.")
    data[[task$id_col]] <- seq_len(nrow(data))
  }

  # Fit Censoring Model on Calibration Data?
  # Standard Conformal: Censoring model fits on "Proper Training Set".
  # Here we accept a fitted object. We assume 'object' is trained on 'Proper Training Set'.
  # The censoring model G should ideally be fitted on the same training set (or full set?).
  # Usually G is estimated on the training set.
  # But we don't have the training data easily accessible except via task$data.
  # Let's fit G on the TRAINING data (task$data) as per standard IPCW approaches.

  train_data <- task$data
  censoring_model <- ml4t2e_fit_censoring(
    data = train_data,
    time_col = task$time_col,
    event_col = task$event_col
  )

  time_grid <- object$time_grid
  outcome_type <- object$outcome_type
  conformal_scores <- list()

  for (engine in object$model_names) {
    if (outcome_type == "survival") {
      # We need predictions on calibration data
      preds <- predict(object, newdata = data, times = time_grid, type = "survival", include = engine)
      pred_mat <- ml4t2e_reshape_preds_to_matrix(preds, data, time_grid, "surv", id_col = task$id_col)
      scores <- .compute_scores_core(pred_mat, data, time_grid, censoring_model, task, event_of_interest = NULL)
    } else {
      # Competing risks
      preds <- predict(object, newdata = data, times = time_grid, type = "cif", include = engine)
      causes <- unique(preds$cause)
      scores_list <- list()
      for (cause_val in causes) {
        cause_preds <- preds[preds$cause == cause_val, ]
        pred_mat <- ml4t2e_reshape_preds_to_matrix(cause_preds, data, time_grid, "cif", id_col = task$id_col)
        scores_list[[as.character(cause_val)]] <- .compute_scores_core(
          pred_mat, data, time_grid, censoring_model, task, cause_val
        )
      }
      scores <- list(scores = scores_list) # Wrapper to match structure
    }
    conformal_scores[[engine]] <- scores
  }

  object$conformal <- list(
    scores = conformal_scores,
    censoring_model = censoring_model,
    cal_data_size = nrow(data)
  )

  object
}
