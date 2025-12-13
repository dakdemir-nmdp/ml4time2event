.cr_prepare_prediction_matrix <- function(Predictions,
                                          n_obs,
                                          pred_times = NULL,
                                          default_time = NULL,
                                          context = "Predictions") {
  original_obj <- Predictions

  if (is.data.frame(Predictions)) {
    Predictions <- as.matrix(Predictions)
  } else if (is.vector(Predictions)) {
    if (length(Predictions) != n_obs) {
      rlang::abort(sprintf(
        "%s vector must have length equal to the number of observations (%d).",
        context, n_obs
      ))
    }
    Predictions <- matrix(Predictions, nrow = 1)
  } else {
    Predictions <- as.matrix(Predictions)
  }

  if (!is.matrix(Predictions)) {
    rlang::abort(sprintf("%s must be coercible to a matrix.", context))
  }

  mat <- Predictions

  if (ncol(mat) == n_obs) {
    # Already time-by-observation
  } else {
    rlang::abort(sprintf(
      "%s must have rows=times and columns=observations. Got %dx%d, expected nrow=%d (times), ncol=%d (observations).",
      context, nrow(mat), ncol(mat), length(pred_times) %||% nrow(mat), n_obs
    ))
  }

  times <- pred_times
  if (is.null(times)) {
    times <- attr(original_obj, "Times")
  }
  if (is.null(times)) {
    times <- attr(original_obj, "times")
  }
  if (is.null(times)) {
    rownames_numeric <- suppressWarnings(as.numeric(rownames(mat)))
    if (!any(is.na(rownames_numeric))) {
      times <- rownames_numeric
    }
  }
  if (is.null(times)) {
    if (nrow(mat) == 1 && !is.null(default_time)) {
      times <- default_time
    }
  }
  if (is.null(times)) {
    rlang::abort(sprintf(
      "Unable to determine the time grid for %s. Provide 'pred_times' or supply numeric rownames/'Times' attribute.",
      tolower(context)
    ))
  }

  times <- suppressWarnings(as.numeric(times))
  if (any(is.na(times))) {
    rlang::abort(sprintf("'pred_times' must be numeric for %s.", tolower(context)))
  }

  if (length(times) != nrow(mat)) {
    rlang::abort(sprintf(
      "'pred_times' must have length %d to match the number of rows in %s (rows represent time points).",
      nrow(mat), tolower(context)
    ))
  }

  list(matrix = mat, times = times)
}

.cr_select_prediction_row <- function(pred_matrix,
                                      pred_times,
                                      eval_time,
                                      context = "Predictions") {
  if (is.null(pred_times) || length(pred_times) != nrow(pred_matrix)) {
    rlang::abort(sprintf(
      "Prediction times for %s must be provided and match the number of matrix rows.",
      tolower(context)
    ))
  }

  idx <- which.min(abs(pred_times - eval_time))
  list(
    pred = as.numeric(pred_matrix[idx, ]),
    matched_time = pred_times[idx],
    times = pred_times
  )
}

#' Time-Dependent Concordance Index for Competing Risks
#'
#' Computes cause-specific concordance index for competing risks data.
#' Measures how well predicted cumulative incidence functions rank pairs
#' of observations where one had the competing event of interest.
#'
#' **Formula**:
#' For cause \eqn{j} at time \eqn{t}:
#' \eqn{C_j(t) = \frac{\sum_{i,j} I(t_i < t_j) \cdot I(\hat{CIF}_{j,i}(t) > \hat{CIF}_{j,j}(t)) \cdot I(\delta_j = j)}{\sum_{i,j} I(t_i < t_j) \cdot I(\delta_j = j)}}
#'
#' **Interpretation**:
#' - Range: 0.5 (random) to 1.0 (perfect)
#' - Computed only among pairs where the event of interest occurred
#' - Competing events exclude observations from risk set
#'
#' @param SurvObj Surv object with structure: `Surv(time, status=cause)`
#' @param Predictions Matrix of predicted CIFs (rows = times, cols = observations)
#' @param time Scalar evaluation time
#' @param cause Integer cause of interest (default: 1)
#' @param TestMat Optional test matrix (currently unused)
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar C-index value at specified time, or NA if no pairs comparable
#'
#' @references
#' Fine, J. P., & Gray, R. J. (1999). "A proportional hazards model for the
#' subdistribution of a competing risk." *Journal of the American Statistical Association*,
#' 94(446), 496–509.
#'
#' @keywords internal
#'
timedepConcordanceCR <- function(SurvObj, Predictions, eval_times = NULL, cause, TestMat = NULL, pred_times = NULL, time = NULL, cens_model = NULL) {
  if (!is.null(time) && is.null(eval_times)) {
    message("Parameter 'time' is deprecated. Please use 'eval_times' instead.")
    eval_times <- time
  }
  if (is.null(eval_times)) rlang::abort("'eval_times' must be specified")

  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(eval_times)) {
    rlang::abort("'eval_times' must be numeric")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  obstimes <- SurvObj[, "time"]
  obsevents <- SurvObj[, "status"]
  n_obs <- length(obstimes)

  # Fit Censoring Model
  if (is.null(cens_model)) {
    cens_model <- .fit_censoring_km(obstimes, obsevents)
  }

  # G(T_i) for all subjects
  G_Ti <- .get_censoring_probs(cens_model, obstimes)
  G_Ti <- pmax(G_Ti, 1e-5)

  # Function for single time point
  calc_auc_single <- function(t) {
    # Get Predictions at t
    pred_prep <- .cr_prepare_prediction_matrix(
      Predictions = Predictions,
      n_obs = n_obs,
      pred_times = pred_times,
      default_time = t,
      context = "Predictions"
    )
    pred_row <- .cr_select_prediction_row(
      pred_matrix = pred_prep$matrix,
      pred_times = pred_prep$times,
      eval_time = t,
      context = "Predictions"
    )
    CIF_hat <- pred_row$pred

    # Define Cases (Event 1 <= t)
    idx_cases <- which(obsevents == cause & obstimes <= t)

    # Define Controls (Not Case, but Observed)
    # 1. Survivors > t
    # 2. Competing Events <= t (Event != cause & Event != 0)
    # Assumes obsevents contain 0=Censored, 1..k=Events.

    idx_surv <- which(obstimes > t)
    idx_comp <- which(obsevents != 0 & obsevents != cause & obstimes <= t)
    idx_controls <- c(idx_surv, idx_comp)

    if (length(idx_cases) == 0 || length(idx_controls) == 0) {
      return(NA_real_)
    }

    # Weights
    w_cases <- 1 / G_Ti[idx_cases]

    w_controls <- numeric(length(idx_controls))

    # Fill weights for controls
    # Survivors: 1/G(t)
    G_t <- .get_censoring_probs(cens_model, t)
    G_t <- max(G_t, 1e-5)

    # Map back to indices for vectorized fill? easier to just iterate or subset
    # Let's rebuild weight vector aligned with idx_controls

    # Create temp dataframe or vectors
    ctrl_times <- obstimes[idx_controls]
    # ctrl_status <- obsevents[idx_controls] # Not strictly needed for weights, but good for context
    ctrl_weights <- numeric(length(idx_controls))

    # For Survivors (> t): Use G(t)
    is_surv <- ctrl_times > t
    ctrl_weights[is_surv] <- 1 / G_t

    # For Competing (<= t): Use G(T)
    # They are observed events (censoring status 1 for C).
    is_comp <- !is_surv
    ctrl_weights[is_comp] <- 1 / G_Ti[idx_controls[is_comp]]

    w_controls <- ctrl_weights

    # Comparison
    # Cases have HIGH CIF (Prediction). Controls have LOW CIF.
    # Concordant if CIF_case > CIF_control.

    risk_cases <- CIF_hat[idx_cases]
    risk_controls <- CIF_hat[idx_controls]

    # Vectorized AUC
    comp_mat <- outer(risk_cases, risk_controls, function(r_c, r_ctrl) {
      ifelse(r_c > r_ctrl, 1, ifelse(r_c == r_ctrl, 0.5, 0))
    })

    w_mat <- outer(w_cases, w_controls, "*")

    num <- sum(comp_mat * w_mat, na.rm = TRUE)
    den <- sum(w_mat, na.rm = TRUE)

    if (den == 0) {
      return(NA_real_)
    }
    num / den
  }

  if (length(eval_times) > 1) {
    sapply(eval_times, calc_auc_single)
  } else {
    calc_auc_single(eval_times[1])
  }
}


#' Brier Score for Competing Risks Predictions
#'
#' Computes mean squared error between predicted cumulative incidence functions
#' and observed event status for a specified cause. Accounts for competing events.
#'
#' **Formula**:
#' \eqn{BS_j(t) = E[(CIF_j(t|X) - Y_j(t))^2]}
#'
#' where \eqn{CIF_j(t|X)} is predicted CIF for cause \eqn{j} and \eqn{Y_j(t)}
#' is binary indicator of cause \eqn{j} occurring by time \eqn{t}.
#'
#' **Interpretation**:
#' - Range: 0 (perfect) to 1 (worst)
#' - Combines calibration and discrimination
#' - Computed only for observations at risk (event of interest or censored/competing)
#'
#' @param SurvObj Surv object with competing risks structure
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Scalar or vector of evaluation time(s)
#' @param cause Integer cause of interest (default: 1)
#' @param TestMat Optional test matrix (currently unused)
#' @param pred_times Numeric vector of time grid from predictions
#' @param time Deprecated parameter for backward compatibility; use `eval_times`
#'
#' @return Numeric Brier score at specified time(s)
#'
#' @keywords internal
#'
BrierScoreCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL, time = NULL, cens_model = NULL) {
  # Handle backward compatibility: support both 'time' and 'eval_times'
  if (!is.null(time) && is.null(eval_times)) {
    # Warn user about deprecation
    message("Parameter 'time' is deprecated. Please use 'eval_times' instead.")
    eval_times <- time
  }

  if (is.null(eval_times)) {
    rlang::abort("'eval_times' must be specified (or 'time' for backward compatibility)")
  }

  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(eval_times)) {
    rlang::abort("'eval_times' must be numeric")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]
  obsevents <- SurvObj[, "status"]

  n_obs <- length(obstimes)

  # Fit or use censoring model
  if (is.null(cens_model)) {
    # 0 = censored, 1..k = events.
    # .fit_censoring_km treats events==0 as the event of interest (censoring).
    # This works for CR data too.
    cens_model <- .fit_censoring_km(obstimes, obsevents)
  }

  # Pre-calculate weigths
  G_Ti <- .get_censoring_probs(cens_model, obstimes)
  G_Ti <- pmax(G_Ti, 1e-5)

  # Vectorized Brier Calculation
  # We loop over eval_times

  calc_brier_single <- function(t) {
    pred_prep <- .cr_prepare_prediction_matrix(
      Predictions = Predictions,
      n_obs = n_obs,
      pred_times = pred_times,
      default_time = t,
      context = "Predictions"
    )
    pred_row <- .cr_select_prediction_row(
      pred_matrix = pred_prep$matrix,
      pred_times = pred_prep$times,
      eval_time = t,
      context = "Predictions"
    )
    CIF_hat <- pred_row$pred

    # G(t)
    G_t <- .get_censoring_probs(cens_model, t)
    G_t <- pmax(G_t, 1e-5)

    # Weights and Targets
    # Case 1: T_i <= t, Delta_i != 0 (Any event observed) -> Weight 1/G(Ti)
    # Case 2: T_i > t (Survivor) -> Weight 1/G(t)
    # Case 3: T_i <= t, Delta_i == 0 (Censored) -> Weight 0

    weights <- numeric(n_obs)

    # Case 1
    idx_event_any <- obstimes <= t & obsevents != 0
    weights[idx_event_any] <- 1 / G_Ti[idx_event_any]

    # Case 2
    idx_surv <- obstimes > t
    weights[idx_surv] <- 1 / G_t

    # Target: I(T_i <= t, Delta_i == cause)
    # 1 if event of interest happened by t
    # 0 otherwise (survived past t, OR competing event happened by t)

    target <- numeric(n_obs)
    target[obstimes <= t & obsevents == cause] <- 1

    # Squared Error
    sq_err <- (target - CIF_hat)^2

    # Weighted Score (exclude NAs in prediction)
    valid_pred <- !is.na(CIF_hat)
    if (sum(valid_pred) == 0) {
      return(NA_real_)
    }

    score_contrib <- weights[valid_pred] * sq_err[valid_pred]
    mean(score_contrib)
  }

  if (length(eval_times) > 1) {
    sapply(eval_times, calc_brier_single)
  } else {
    calc_brier_single(eval_times[1])
  }
}


#' Integrated Concordance Index for Competing Risks
#'
#' Averages time-dependent concordance across evaluation times.
#' Provides a single summary measure of discrimination over follow-up.
#'
#' **Interpretation**:
#' - Mean of C-index values across time points
#' - Range: 0.5 (random) to 1.0 (perfect)
#' - Accounts for competing events in all comparisons
#'
#' @param SurvObj Surv object with competing risks
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Vector of time points at which to evaluate
#' @param cause Integer cause of interest
#' @param TestMat Optional test matrix
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar integrated C-index value
#'
#' @keywords internal
#'
integratedConcordanceCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]

  n_obs <- length(obstimes)
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    context = "Predictions"
  )
  pred_matrix <- pred_prep$matrix
  pred_times <- pred_prep$times

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }

  # Calculate concordance at each time point
  concordance_values <- sapply(eval_times, function(t) {
    timedepConcordanceCR(
      SurvObj = SurvObj,
      Predictions = pred_matrix,
      time = t,
      cause = cause,
      TestMat = TestMat,
      pred_times = pred_times
    )
  })

  # Remove NA values
  valid_idx <- !is.na(concordance_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  concordance_values <- concordance_values[valid_idx]
  # eval_times_valid <- eval_times[valid_idx]  # Not used

  # Calculate integrated concordance (mean over time)
  integrated_c <- mean(concordance_values, na.rm = TRUE)

  return(integrated_c)
}


#' Integrated Brier Score for Competing Risks
#'
#' Integrates cause-specific Brier scores over time using trapezoidal rule.
#' Provides a single summary measure of calibration and discrimination.
#'
#' **Formula**:
#' \eqn{IBS_j = \frac{1}{t_{max} - t_{min}} \int_0^{t_{max}} BS_j(t) dt}
#'
#' @param SurvObj Surv object with competing risks
#' @param Predictions Matrix of predicted CIFs
#' @param eval_times Vector of evaluation times. If NULL, uses time grid from predictions.
#' @param cause Integer cause of interest
#' @param TestMat Optional test matrix
#' @param pred_times Numeric vector of time grid from predictions
#'
#' @return Scalar integrated Brier score
#'
#' @keywords internal
#'
integratedBrierCR <- function(SurvObj, Predictions, eval_times = NULL, cause = 1, TestMat = NULL, pred_times = NULL) {
  # Input validation
  if (!inherits(SurvObj, "Surv")) {
    rlang::abort("'SurvObj' must be a Surv object")
  }
  if (!is.numeric(cause) || length(cause) != 1) {
    rlang::abort("'cause' must be a single numeric value")
  }

  # Extract time and event from Surv object
  obstimes <- SurvObj[, "time"]

  n_obs <- length(obstimes)
  pred_prep <- .cr_prepare_prediction_matrix(
    Predictions = Predictions,
    n_obs = n_obs,
    pred_times = pred_times,
    context = "Predictions"
  )
  pred_matrix <- pred_prep$matrix
  pred_times <- pred_prep$times

  if (is.null(eval_times)) {
    eval_times <- pred_times
  }

  # Calculate Brier score at each time point
  # Use the non-deprecated 'eval_times' argument to avoid warnings
  brier_values <- sapply(eval_times, function(t) {
    BrierScoreCR(
      SurvObj = SurvObj,
      Predictions = pred_matrix,
      eval_times = t,
      cause = cause,
      TestMat = TestMat,
      pred_times = pred_times
    )
  })

  # Remove NA values
  valid_idx <- !is.na(brier_values)
  if (sum(valid_idx) == 0) {
    return(NA)
  }

  brier_values <- brier_values[valid_idx]
  eval_times_valid <- eval_times[valid_idx]

  # Calculate integrated Brier score (area under curve using trapezoidal rule)
  if (length(eval_times_valid) < 2) {
    return(brier_values[1])
  }

  time_diffs <- diff(eval_times_valid)
  avg_brier <- (brier_values[-1] + brier_values[-length(brier_values)]) / 2
  integrated_bs <- sum(time_diffs * avg_brier) / (max(eval_times_valid) - min(eval_times_valid))

  return(integrated_bs)
}


#'
#'
restrictedMeanTimeLostCR <- function(Predictions, times, UL, LL = 0) {
  # Input validation
  if (!is.matrix(Predictions)) {
    rlang::abort("'Predictions' must be a matrix")
  }
  if (!is.numeric(times)) {
    rlang::abort("'times' must be numeric")
  }
  if (length(times) != nrow(Predictions)) {
    rlang::abort("'times' must match the number of rows in Predictions (rows represent time points)")
  }
  if (!is.numeric(UL) || length(UL) != 1 || UL <= LL) {
    rlang::abort("'UL' must be a single numeric value greater than 'LL'")
  }
  if (!is.numeric(LL) || length(LL) != 1 || LL < 0) {
    rlang::abort("'LL' must be a single non-negative numeric value")
  }

  # Ensure times are sorted
  if (!all(diff(times) >= 0)) {
    rlang::abort("'times' must be sorted in ascending order")
  }

  # Calculate ETL using unified function
  time_lost_vector <- calculate_expected_time_lost(
    times = times,
    event_probs = Predictions,
    upper_limit = UL,
    lower_limit = LL
  )

  return(time_lost_vector)
}
