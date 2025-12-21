ml4t2e_align_surv_predictions <- function(predsurv, pred_times, obstimes, context = "predsurv") {
  if (is.null(predsurv)) {
    rlang::abort(sprintf("'%s' cannot be NULL.", context))
  }

  pred_matrix <- as.matrix(predsurv)
  if (!is.numeric(pred_matrix)) {
    rlang::abort(sprintf("'%s' must be numeric.", context))
  }
  if (any(!is.finite(pred_matrix))) {
    rlang::abort(sprintf("'%s' contains non-finite values.", context))
  }

  if (missing(pred_times) || is.null(pred_times)) {
    rlang::abort("'pred_times' must be supplied and cannot be NULL.")
  }
  if (!is.numeric(pred_times)) {
    rlang::abort("'pred_times' must be numeric.")
  }
  if (length(pred_times) == 0) {
    rlang::abort("'pred_times' cannot be empty.")
  }

  obstimes <- as.numeric(obstimes)
  if (!all(is.finite(obstimes))) {
    rlang::abort("'obstimes' must be numeric and finite.")
  }
  n_obs <- length(obstimes)

  n_times <- length(pred_times)
  dims <- dim(pred_matrix)

  orientation_fixed <- FALSE
  if (nrow(pred_matrix) == n_times && ncol(pred_matrix) == n_obs) {
    orientation_fixed <- TRUE
  } else if (ncol(pred_matrix) == n_times && nrow(pred_matrix) == n_obs) {
    pred_matrix <- t(pred_matrix)
    orientation_fixed <- TRUE
  }

  if (!orientation_fixed) {
    rlang::abort(sprintf(
      "Unable to align predictions. Expected matrix with %d time points and %d observations, got %dx%d.",
      n_times, n_obs, dims[1], dims[2]
    ))
  }

  # Ensure times are strictly increasing for interpolation
  if (any(is.na(pred_times))) {
    rlang::abort("'pred_times' contains NA values.")
  }
  order_idx <- order(pred_times)
  pred_times_sorted <- pred_times[order_idx]
  if (any(diff(pred_times_sorted) == 0)) {
    keep_idx <- !duplicated(pred_times_sorted)
    pred_times_sorted <- pred_times_sorted[keep_idx]
    pred_matrix <- pred_matrix[order_idx, , drop = FALSE][keep_idx, , drop = FALSE]
  } else {
    pred_matrix <- pred_matrix[order_idx, , drop = FALSE]
  }

  rownames(pred_matrix) <- NULL

  list(
    matrix = pred_matrix,
    times = pred_times_sorted,
    obstimes = obstimes,
    n_obs = n_obs
  )
}

ml4t2e_interpolate_survival <- function(pred_matrix, pred_times, eval_time) {
  if (length(eval_time) != 1L || !is.finite(eval_time)) {
    rlang::abort("'eval_time' must be a single finite numeric value.")
  }

  if (eval_time <= pred_times[1]) {
    return(pred_matrix[1, ])
  }
  if (eval_time >= pred_times[length(pred_times)]) {
    return(pred_matrix[nrow(pred_matrix), ])
  }

  apply(pred_matrix, 2, function(curve) {
    if (all(is.na(curve))) {
      return(NA_real_)
    }
    stats::approx(
      x = pred_times,
      y = curve,
      xout = eval_time,
      method = "linear",
      rule = 2
    )$y
  })
}

ml4t2e_validate_events <- function(obsevents) {
  if (is.null(obsevents)) {
    rlang::abort("'obsevents' cannot be NULL.")
  }
  if (!is.numeric(obsevents) && !is.logical(obsevents)) {
    rlang::abort("'obsevents' must be numeric or logical.")
  }
  ev <- as.numeric(obsevents)
  ev[is.na(ev)] <- NA_real_
  unique_vals <- unique(ev[!is.na(ev)])
  if (!all(unique_vals %in% c(0, 1))) {
    warning("'obsevents' contains values other than 0/1. Treating non-1 values as censored.")
    ev[!is.na(ev) & ev != 1] <- 0
  }
  ev
}

ml4t2e_cindex_at_time <- function(pred_surv, eval_time, obstimes, obsevents, G_Ti, G_t) {
  # Define Cases and Controls
  # Cases: Event occurred by eval_time
  idx_cases <- which(obsevents == 1 & obstimes <= eval_time)

  # Controls: Survived past eval_time
  idx_controls <- which(obstimes > eval_time)

  if (length(idx_cases) == 0 || length(idx_controls) == 0) {
    return(NA_real_)
  }

  weights_cases <- 1 / G_Ti[idx_cases]
  weights_controls <- rep(1 / G_t, length(idx_controls))

  # Risk scores (Higher risk = Lower survival)
  risk <- 1 - pred_surv

  risk_cases <- risk[idx_cases]
  risk_controls <- risk[idx_controls]

  # Calculate Weighted AUC
  # Numerator: Sum_{i in cases} Sum_{j in controls} W_i * W_j * I(Risk_i > Risk_j)
  # Denominator: Sum_{i in cases} Sum_{j in controls} W_i * W_j

  # Explicit double loop (O(N_cases * N_controls)) - same complexity as before roughly
  # Vectorized in R is better but takes memory O(N*M).
  # If N=2000, N*M = 4e6, manageable.

  numerator <- 0

  # Vectorized comparison:
  # Outer product comparison
  # comp_mat[i, j] = 1 if risk_cases[i] > risk_controls[j]
  #                  0.5 if tie
  #                  0 if <

  comp_mat <- outer(risk_cases, risk_controls, function(r_c, r_ctrl) {
    ifelse(r_c > r_ctrl, 1, ifelse(r_c == r_ctrl, 0.5, 0))
  })

  # Weights matrix
  w_mat <- outer(weights_cases, weights_controls, "*")

  numerator <- sum(comp_mat * w_mat, na.rm = TRUE)
  denominator <- sum(w_mat, na.rm = TRUE)

  if (denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

#' Time-Dependent Concordance Index for Survival Analysis
#'
#' Computes Harrell's C-index at each time point in the prediction grid.
#' This assesses how well predicted survival probabilities discriminate between
#' pairs of observations where one had an earlier event.
#'
#' **Formula**:
#' At time \eqn{t}, \eqn{C(t) = \frac{\sum_{i,j} I(t_i < t_j) \cdot I(\hat{S}_i(t) > \hat{S}_j(t)) \cdot \delta_j}{\sum_{i,j} I(t_i < t_j) \cdot \delta_j}}
#'
#' **Interpretation**:
#' - Range: 0.5 (random predictions) to 1.0 (perfect discrimination)
#' - 0.5 = worse than random; 1.0 = perfect ranking
#' - Values > 0.7 typically considered good discrimination
#'
#' @param predsurv Matrix of predicted survival probabilities (rows = times, cols = observations)
#' @param pred_times Numeric vector of time points at which predictions were made
#' @param obstimes Observed event/censoring times
#' @param obsevents Binary event indicator (1 = event, 0 = censored)
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return List of class `ml4time2event_cindex` containing:
#'   - `AppCindex$matrix`: Vector of C-index values at each time point
#'   - `AppCindex$time`: Corresponding time points
#'   - `time`: Time grid
#'   - `call`: Original function call
#'
#' @references
#' Harrell, F. E., Lee, K. L., & Mark, D. B. (1996).
#' "Multivariable prognostic models." *Statistics in Medicine*, 15(4), 361–387.
#'
#' @keywords internal
#'

timedepConcordance <- function(predsurv, pred_times, obstimes, obsevents, TestMat = NULL, cens_model = NULL) {
  alignment <- ml4t2e_align_surv_predictions(predsurv, pred_times, obstimes, context = "predsurv")
  obsevents_numeric <- ml4t2e_validate_events(obsevents)
  obstimes_sorted <- alignment$obstimes

  # Fit censoring model if needed
  if (is.null(cens_model)) {
    cens_model <- .fit_censoring_km(obstimes_sorted, obsevents_numeric)
  }

  # Pre-calculate G(T_i) for all subjects
  G_Ti <- .get_censoring_probs(cens_model, obstimes_sorted)
  G_Ti <- pmax(G_Ti, 1e-5)

  c_values <- vapply(
    seq_along(alignment$times),
    function(idx) {
      t <- alignment$times[idx]
      surv_at_time <- alignment$matrix[idx, ]

      # G(t) for controls
      G_t <- .get_censoring_probs(cens_model, t)
      G_t <- max(G_t, 1e-5)

      ml4t2e_cindex_at_time(
        pred_surv = surv_at_time,
        eval_time = t,
        obstimes = obstimes_sorted,
        obsevents = obsevents_numeric,
        G_Ti = G_Ti,
        G_t = G_t
      )
    },
    numeric(1)
  )

  result <- list(
    AppCindex = list(
      matrix = c_values,
      time = alignment$times
    ),
    time = alignment$times,
    call = match.call()
  )
  class(result) <- c("ml4time2event_cindex", "pecCindex")
  result
}



# Internal helpers for IPCW
#' @noRd
.fit_censoring_km <- function(times, events) {
  # Fit KM for censoring (status=0 is event)
  # 0 = censored for survival, so status = 0 is the event of interest for censoring
  # We suppress warnings about potential singular computations if data is small
  suppressWarnings(survival::survfit(survival::Surv(times, events == 0) ~ 1))
}

.get_censoring_probs <- function(km_fit, times) {
  # Predict P(C > t)
  # Use stepfun for robust interpolation (constant survival between events)

  if (length(km_fit$time) == 0) {
    # No events/censoring? (e.g. empty fit or all 0 time?)
    return(rep(1, length(times)))
  }

  # stepfun(x, y). y must be length(x) + 1.
  # y[1] is value < x[1]. y[2] is value >= x[1].
  # We want P(C > t).
  # Before first time, prob is 1.
  # After first time (if event), prob drops.

  stats::stepfun(km_fit$time, c(1, km_fit$surv), right = FALSE)(times)
}

#' Brier Score for Survival Predictions (IPCW Corrected)
#'
#' Computes the Inverse Probability of Censoring Weighted (IPCW) Brier Score.
#' Measures mean squared error between predicted survival probabilities and
#' observed binary event status, weighted to account for censoring.
#'
#' **Formula**:
#' \eqn{BS(t) = \frac{1}{N} \sum_{i=1}^N \hat{W}_i(t) (Y_i(t) - \hat{S}(t|X_i))^2}
#'
#' Weights \eqn{\hat{W}_i(t)} definition:
#' - If \eqn{T_i \le t} and \eqn{\delta_i = 1} (Observed event): \eqn{1/\hat{G}(T_i-)}
#' - If \eqn{T_i > t} (Observed survivor): \eqn{1/\hat{G}(t)}
#' - Otherwise (Censored before \eqn{t}): \eqn{0}
#'
#' Where \eqn{\hat{G}} is the Kaplan-Meier estimate of the censoring distribution.
#'
#' @param predsurv Matrix of predicted survival probabilities (rows = times, cols = observations)
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param eval_times Optional evaluation time points. If NULL, computed at all `pred_times`.
#' @param TestMat Optional test matrix (currently unused)
#' @param cens_model Optional pre-fitted censoring model (survfit object). If NULL, computed from provided data.
#'
#' @return List of class `ml4time2event_brier` containing Brier scores and time grid
#'
#' @references
#' Graf, E., Schmoor, C., Sauerbrei, W., & Schumacher, M. (1999).
#' "Assessment and comparison of prognostic classification schemes for survival data."
#' *Statistics in Medicine*, 18(17‐18), 2529-2545.
#'
#' @keywords internal
#'
BrierScore <- function(predsurv, pred_times, obstimes, obsevents,
                       eval_times = NULL, TestMat = NULL, cens_model = NULL) {
  alignment <- ml4t2e_align_surv_predictions(predsurv, pred_times, obstimes, context = "predsurv")
  if (is.null(eval_times)) {
    eval_times <- alignment$times
  } else {
    if (!is.numeric(eval_times)) {
      rlang::abort("'eval_times' must be numeric.")
    }
    if (any(!is.finite(eval_times))) {
      rlang::abort("'eval_times' must contain only finite values.")
    }
  }

  obsevents_numeric <- ml4t2e_validate_events(obsevents)
  obstimes_sorted <- alignment$obstimes # Aligned obstimes
  n_obs <- length(obstimes_sorted)

  # fit or use censoring model
  if (is.null(cens_model)) {
    cens_model <- .fit_censoring_km(obstimes_sorted, obsevents_numeric)
  }

  # Pre-calculate G(Ti) for all subjects
  # We use a small epsilon to avoid division by zero if G(t)=0
  # (though usually G(t) > 0 for valid follow-up range)
  G_Ti <- .get_censoring_probs(cens_model, obstimes_sorted)
  G_Ti <- pmax(G_Ti, 1e-5)

  # Pre-calculate G(t) for all eval_times
  G_t_all <- .get_censoring_probs(cens_model, eval_times)
  G_t_all <- pmax(G_t_all, 1e-5)

  brier_values <- vapply(
    seq_along(eval_times),
    function(j) {
      t_eval <- eval_times[j]
      G_t <- G_t_all[j]

      pred_surv <- ml4t2e_interpolate_survival(alignment$matrix, alignment$times, t_eval)

      # Determine weights and contributions
      # Case 1: T_i <= t, Delta_i = 1 -> Weight 1/G(Ti)
      # Case 2: T_i > t -> Weight 1/G(t)
      # Case 3: T_i <= t, Delta_i = 0 -> Weight 0

      weights <- numeric(n_obs)

      # Case 1
      idx_event <- obstimes_sorted <= t_eval & obsevents_numeric == 1
      weights[idx_event] <- 1 / G_Ti[idx_event]

      # Case 2
      idx_surv <- obstimes_sorted > t_eval
      weights[idx_surv] <- 1 / G_t

      # Observations outcomes Y_i(t)
      # Y_i(t) = 1 if T_i > t (Survivor)
      # Y_i(t) = 0 if T_i <= t (Event)
      # (Note: Standard Brier is (Y - S)^2. If S is Prob(T>t), then Y should be I(T>t))

      obs_Y <- as.numeric(idx_surv)

      # Squared errors
      sq_err <- (obs_Y - pred_surv)^2

      # Weighted Mean Square Error
      # Exclude those with 0 weight from the mean?
      # Graf et al definition uses 1/N sum(Weights * SquaredError).
      # So we sum and divide by N (total sample size), NOT sum of weights.

      # Handle NAs in prediction
      valid_pred <- !is.na(pred_surv)
      if (sum(valid_pred) == 0) {
        return(NA_real_)
      }

      # If predictions have NAs, we should theoretically ignore them.
      # If we ignore `k` predictions, we divide by `N-k`.

      score_contrib <- weights[valid_pred] * sq_err[valid_pred]
      mean(score_contrib)
    },
    numeric(1)
  )

  result <- list(
    AppErr = list(
      model = brier_values,
      time = eval_times
    ),
    time = eval_times,
    call = match.call()
  )
  class(result) <- c("ml4time2event_brier", "pec")
  result
}


#' Integrated Brier Score (IBS)
#'
#' Integrates Brier scores over the time range using trapezoidal rule.
#' Provides a single summary measure of prediction error across all time points.
#'
#' **Formula**:
#' \eqn{IBS = \frac{1}{t_{max} - t_{min}} \int_0^{t_{max}} BS(t) dt}
#'
#' Computed numerically via trapezoidal rule over the time grid.
#'
#' **Interpretation**:
#' - Range: 0 (perfect) to 1 (worst)
#' - Combines calibration across entire follow-up
#' - Weights later times equally (unlike summary C-index)
#'
#' @param predsurv Matrix of predicted survival probabilities
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param eval_times Optional evaluation times. If NULL, uses `pred_times`.
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return Scalar integrated Brier score value
#'
#' @references
#' Mogensen, U. B., Ishwaran, H., & Gerds, T. A. (2012).
#' "Evaluating random forests for survival analysis using prediction error curves."
#' *Journal of Statistical Software*, 50(11), 1–23.
#'
#' @keywords internal
#'
integratedBrier <- function(predsurv, pred_times, obstimes, obsevents,
                            eval_times = NULL, TestMat = NULL) {
  # Get Brier scores at all time points
  brier_obj <- BrierScore(
    predsurv = predsurv,
    pred_times = pred_times,
    obstimes = obstimes,
    obsevents = obsevents,
    eval_times = eval_times,
    TestMat = TestMat
  )

  # Extract the Brier score values
  # pec objects have AppErr slot with the prediction error
  if (is.null(brier_obj$AppErr)) {
    rlang::abort("BrierScore output is missing the 'AppErr' component.")
  }
  brier_values <- brier_obj$AppErr$model
  eval_times <- brier_obj$time

  # Integrate using trapezoidal rule
  if (length(eval_times) < 2) {
    return(brier_values[1])
  }

  # Calculate integrated Brier score (area under curve)
  time_diffs <- diff(eval_times)
  avg_brier <- (brier_values[-1] + brier_values[-length(brier_values)]) / 2
  integrated_bs <- sum(time_diffs * avg_brier) / (max(eval_times) - min(eval_times))

  return(integrated_bs)
}


#' Integrated Concordance Index
#'
#' Averages time-dependent C-index values over the time grid.
#' Provides a summary measure of discrimination across follow-up.
#'
#' **Interpretation**:
#' - Average of all C-index values at supplied time points
#' - Single summary value (unlike time-dependent C-index)
#' - Range: 0.5 to 1.0
#'
#' @param predsurv Matrix of predicted survival probabilities
#' @param pred_times Numeric vector of time points
#' @param obstimes Observed times
#' @param obsevents Binary event indicator
#' @param TestMat Optional test matrix (currently unused)
#'
#' @return Scalar integrated C-index value
#'
#' @keywords internal
#'
integratedC <- function(predsurv, pred_times, obstimes, obsevents, TestMat = NULL) {
  # Get time-dependent concordance
  cindex_obj <- timedepConcordance(
    predsurv = predsurv,
    pred_times = pred_times,
    obstimes = obstimes,
    obsevents = obsevents,
    TestMat = TestMat
  )

  # Extract concordance values at each time point
  if (is.null(cindex_obj$AppCindex)) {
    rlang::abort("timedepConcordance output is missing the 'AppCindex' component.")
  }
  cindex_values <- cindex_obj$AppCindex$matrix

  # Calculate mean concordance over time (integrated C-index)
  integrated_c <- mean(cindex_values, na.rm = TRUE)

  return(integrated_c)
}
