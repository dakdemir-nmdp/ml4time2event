#' @title SurvModel_TTAH
#'
#' @description Fit a lightweight discrete-time time-varying additive hazard (TTAH) model
#'   for survival outcomes. The implementation uses spline-augmented feature encoders,
#'   low-rank time interactions, and ridge-stabilised logistic regression for a
#'   fully CPU-friendly workflow.
#'
#' @param data data.frame containing predictors and outcome columns
#' @param expvars character vector of predictor column names
#' @param timevar character name of time-to-event column
#' @param eventvar character name of event indicator (0/1)
#' @param time_grid optional numeric vector of discrete time points; if `NULL`, a grid
#'   is constructed from the observed data using quantiles
#' @param n_time integer, target number of time points when `time_grid` is not supplied
#' @param spline_knots integer, degrees of freedom for numeric spline bases
#' @param latent_dim integer, number of latent interaction factors
#' @param time_basis_df integer, degrees of freedom for the B-spline time basis
#' @param lambda numeric, ridge penalty applied to the logistic fit (excludes intercept)
#' @param engine character, one of `c("auto", "base")`
#' @param verbose logical; if `TRUE`, prints training diagnostics
#'
#' @return list with class `ml4t2e_surv_ttah` containing the fitted model and metadata
#' @export
SurvModel_TTAH <- function(data, expvars, timevar, eventvar,
                           time_grid = NULL, n_time = 50,
                           spline_knots = 5, latent_dim = 8,
                           time_basis_df = 4, lambda = 1e-3,
                           engine = c("auto", "base"),
                           verbose = FALSE) {

  engine <- match.arg(engine)
  if (engine == "auto") {
    engine <- "base"
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!all(c(timevar, eventvar) %in% colnames(data))) {
    stop("Both 'timevar' and 'eventvar' must be present in 'data'.")
  }
  if (!all(expvars %in% colnames(data))) {
    missing <- setdiff(expvars, colnames(data))
    stop("The following predictors are missing from 'data': ", paste(missing, collapse = ", "))
  }

  work_data <- data[, c(timevar, eventvar, expvars), drop = FALSE]
  work_data[[eventvar]] <- as.numeric(work_data[[eventvar]] == 1)

  complete_idx <- complete.cases(work_data)
  if (!all(complete_idx)) {
    removed <- sum(!complete_idx)
    if (verbose) {
      message("Removing ", removed, " rows with missing values for TTAH training.")
    }
    work_data <- work_data[complete_idx, , drop = FALSE]
  }
  if (nrow(work_data) < 20) {
    stop("At least 20 complete observations are required to fit SurvModel_TTAH.")
  }

  varprof <- VariableProfile(work_data, expvars)

  factor_levels <- list()
  for (var in expvars) {
    if (is.factor(work_data[[var]]) || is.character(work_data[[var]])) {
      factor_levels[[var]] <- levels(as.factor(work_data[[var]]))
    }
  }

  obs_times <- work_data[[timevar]]
  obs_events <- work_data[[eventvar]]

  grid <- ttah_build_time_grid(obs_times, time_grid = time_grid, n_time = n_time)

  prep <- ttah_prepare_features(
    data = work_data[, expvars, drop = FALSE],
    expvars = expvars,
    factor_levels = factor_levels,
    basis_specs = NULL,
    spline_knots = spline_knots
  )
  phi <- prep$phi
  basis_specs <- prep$basis_specs
  factor_levels <- prep$factor_levels

  latent_projection <- ttah_compute_latent_projection(phi, latent_dim = latent_dim)
  phi_latent <- if (ncol(latent_projection) > 0) {
    tmp <- phi %*% latent_projection
    colnames(tmp) <- colnames(latent_projection)
    tmp
  } else {
    NULL
  }

  time_basis <- ttah_time_basis(grid, df = time_basis_df)
  time_basis_matrix <- time_basis$matrix
  if (ncol(time_basis_matrix) > 0) {
    colnames(time_basis_matrix) <- paste0("time_basis", seq_len(ncol(time_basis_matrix)))
  }
  time_basis_names <- colnames(time_basis_matrix)
  interval_index <- ttah_assign_intervals(obs_times, grid)

  long_design <- ttah_build_long_design(
    phi = phi,
    phi_latent = phi_latent,
    time_basis = time_basis_matrix,
    interval_index = interval_index,
    event = obs_events
  )

  design_matrix <- cbind(Intercept = 1, long_design$X)

  fit <- ttah_ridge_glm(
    X = design_matrix,
    y = long_design$y,
    lambda = lambda
  )
  if (verbose && !fit$converged) {
    warning("TTAH logistic fit did not converge within the iteration limit.")
  }

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  offset <- 1
  coef_add <- coef[offset + seq_len(ncol(phi))]
  offset <- offset + ncol(phi)
  coef_time <- coef[offset + seq_len(ncol(time_basis_matrix))]
  offset <- offset + ncol(time_basis_matrix)
  latent_coef <- numeric(0)
  if (ncol(latent_projection) > 0) {
    latent_len <- ncol(latent_projection) * ncol(time_basis_matrix)
    if (latent_len > 0) {
      latent_coef <- coef[offset + seq_len(latent_len)]
      offset <- offset + latent_len
    }
  }

  internal_model <- list(
    intercept = coef[1],
    coef_add = coef_add,
    coef_time = coef_time,
    coef_inter = latent_coef,
    feature_names = colnames(phi),
    time_basis_names = time_basis_names,
    time_basis_specs = time_basis$specs,
    design_stats = list(
      lambda = lambda,
      deviance = fit$deviance,
      iterations = fit$iterations,
      converged = fit$converged
    )
  )

  result <- list(
    model = internal_model,
    times = grid,
    varprof = varprof,
    expvars = expvars,
    factor_levels = factor_levels,
    time_grid = grid,
    basis_specs = basis_specs,
    latent_projection = latent_projection,
    engine = engine
  )
  class(result) <- c("ml4t2e_surv_ttah", "list")
  result
}

#' @title Predict_SurvModel_TTAH
#'
#' @description Generate survival probability curves from a fitted TTAH survival model.
#'
#' @param modelout fitted model returned by `SurvModel_TTAH`
#' @param newdata data frame with predictor columns matching training predictors
#' @param new_times optional numeric vector of times for interpolation
#'
#' @return list with components:
#'   \item{Probs}{matrix of survival probabilities (rows = times, cols = observations)}
#'   \item{Times}{vector of times corresponding to the rows of `Probs`}
#' @export
Predict_SurvModel_TTAH <- function(modelout, newdata, new_times = NULL) {
  if (missing(modelout)) stop("'modelout' is required.")
  if (!inherits(modelout, "ml4t2e_surv_ttah")) {
    stop("modelout must be an object returned by SurvModel_TTAH.")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame.")
  }
  if (!all(modelout$expvars %in% colnames(newdata))) {
    missing <- setdiff(modelout$expvars, colnames(newdata))
    stop("The following predictors are missing from newdata: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(new_times)) {
    if (!is.numeric(new_times)) stop("'new_times' must be numeric.")
    if (any(new_times < 0)) stop("'new_times' must be non-negative.")
  }

  newdata_prepared <- newdata[, modelout$expvars, drop = FALSE]

  for (var in modelout$expvars) {
    if (var %in% names(modelout$factor_levels)) {
      levels_target <- modelout$factor_levels[[var]]
      newdata_prepared[[var]] <- factor(newdata_prepared[[var]], levels = levels_target)
      unseen <- setdiff(unique(newdata[[var]]), levels_target)
      if (length(unseen) > 0) {
        warning("Predictor '", var, "' has unseen levels: ",
                paste(unseen, collapse = ", "), ". They are treated as missing.")
      }
    }
  }

  if (anyNA(newdata_prepared)) {
    stop("Missing values detected in predictors after alignment; impute before prediction.")
  }

  prep_new <- ttah_prepare_features(
    data = newdata_prepared,
    expvars = modelout$expvars,
    factor_levels = modelout$factor_levels,
    basis_specs = modelout$basis_specs
  )

  phi_new <- prep_new$phi
  if (!all(modelout$model$feature_names %in% colnames(phi_new))) {
    missing_cols <- setdiff(modelout$model$feature_names, colnames(phi_new))
    stop("Basis columns missing during prediction: ", paste(missing_cols, collapse = ", "))
  }
  phi_new <- phi_new[, modelout$model$feature_names, drop = FALSE]

  latent_proj <- modelout$latent_projection
  latent_dim <- if (is.null(latent_proj)) 0 else ncol(latent_proj)
  phi_latent_new <- if (!is.null(latent_proj) && latent_dim > 0) {
    phi_new %*% latent_proj
  } else {
    NULL
  }

  time_grid <- modelout$time_grid
  time_basis <- ttah_eval_time_basis(time_grid, modelout$model$time_basis_specs)
  if (!is.null(modelout$model$time_basis_names) &&
      length(modelout$model$time_basis_names) == ncol(time_basis)) {
    colnames(time_basis) <- modelout$model$time_basis_names
  }

  n_obs <- nrow(phi_new)
  K <- length(time_grid)

  eta_matrix <- matrix(modelout$model$intercept, nrow = K, ncol = n_obs)

  main_effect <- as.vector(phi_new %*% modelout$model$coef_add)
  eta_matrix <- eta_matrix + matrix(main_effect, nrow = K, ncol = n_obs, byrow = TRUE)

  time_linear <- drop(time_basis %*% modelout$model$coef_time)
  eta_matrix <- eta_matrix + matrix(time_linear, nrow = K, ncol = n_obs, byrow = FALSE)

  if (!is.null(phi_latent_new) && length(modelout$model$coef_inter) > 0) {
    theta_matrix <- matrix(
      modelout$model$coef_inter,
      nrow = ncol(phi_latent_new),
      ncol = ncol(time_basis)
    )
    interaction_term <- phi_latent_new %*% (theta_matrix %*% t(time_basis)) # n_obs x K
    eta_matrix <- eta_matrix + t(interaction_term)
  }

  hazard_matrix <- ttah_sigmoid(eta_matrix)
  survival_body <- matrix(0, nrow = K, ncol = n_obs)
  for (j in seq_len(n_obs)) {
    survival_body[, j] <- cumprod(1 - hazard_matrix[, j])
  }

  survival_body <- pmin(pmax(survival_body, 0), 1)
  Probs <- rbind(rep(1, n_obs), survival_body)
  Times <- c(0, time_grid)

  if (!is.null(new_times)) {
    Probs <- survprobMatInterpolator(Probs, Times, new_times)
    Times <- new_times
  }

  list(
    Probs = Probs,
    Times = Times
  )
}
