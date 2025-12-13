#'
#'
#'
CRModel_TTAH <- function(data, expvars, timevar, eventvar,
                         event_codes = NULL,
                         time_grid = NULL, n_time = 50,
                         spline_knots = 5, latent_dim = 6,
                         time_basis_df = 4, lambda = 1e-3, maxit = 200,
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
    stop("'timevar' and 'eventvar' must be present in 'data'.")
  }
  if (!all(expvars %in% colnames(data))) {
    missing <- setdiff(expvars, colnames(data))
    stop("Predictors missing from 'data': ", paste(missing, collapse = ", "))
  }

  available_events <- sort(unique(data[[eventvar]][data[[eventvar]] != 0]))
  if (length(available_events) == 0) {
    stop("No competing events observed in the training data.")
  }

  if (is.null(event_codes)) {
    cause_codes <- as.character(available_events)
  } else {
    cause_codes <- as.character(event_codes)
    missing_codes <- setdiff(cause_codes, as.character(available_events))
    if (length(missing_codes) > 0) {
      stop("Requested event codes not found in data: ", paste(missing_codes, collapse = ", "))
    }
  }

  work_data <- data[, c(timevar, eventvar, expvars), drop = FALSE]
  complete_idx <- complete.cases(work_data)
  if (!all(complete_idx)) {
    removed <- sum(!complete_idx)
    if (verbose) {
      message("Removing ", removed, " rows with missing values prior to CR TTAH training.")
    }
    work_data <- work_data[complete_idx, , drop = FALSE]
  }
  if (nrow(work_data) < 30) {
    stop("At least 30 complete observations are required to fit CRModel_TTAH.")
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

  interval_index <- ttah_assign_intervals(obs_times, grid)
  multiclass_design <- ttah_build_multiclass_design(
    phi = phi,
    phi_latent = phi_latent,
    time_basis = time_basis_matrix,
    interval_index = interval_index,
    event = obs_events,
    cause_codes = cause_codes
  )

  if (length(unique(multiclass_design$target)) < 2) {
    stop("Insufficient event diversity to train multinomial hazards.")
  }

  design_df <- as.data.frame(multiclass_design$X)
  design_df$target <- multiclass_design$target

  fit <- nnet::multinom(
    target ~ .,
    data = design_df,
    decay = lambda,
    maxit = maxit,
    trace = verbose,
    MaxNWts = 20000
  )

  internal_model <- list(
    multinom = fit,
    feature_columns = setdiff(colnames(design_df), "target"),
    time_basis_specs = time_basis$specs,
    time_basis_names = colnames(time_basis_matrix),
    cause_levels = multiclass_design$levels,
    lambda = lambda,
    maxit = maxit
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
    cause_codes = cause_codes,
    engine = engine
  )
  class(result) <- c("ml4t2e_cr_ttah", "list")
  result
}

#'
#'
#'
Predict_CRModel_TTAH <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {
  if (missing(modelout) || !inherits(modelout, "ml4t2e_cr_ttah")) {
    stop("modelout must be an object returned by CRModel_TTAH.")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame.")
  }
  if (!all(modelout$expvars %in% colnames(newdata))) {
    missing <- setdiff(modelout$expvars, colnames(newdata))
    stop("Predictors missing from newdata: ", paste(missing, collapse = ", "))
  }
  if (!is.null(new_times)) {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values.")
    }
  }

  cause_codes <- modelout$cause_codes
  if (is.null(event_of_interest)) {
    event_of_interest <- cause_codes[1]
  } else {
    event_of_interest <- as.character(event_of_interest)
    if (!event_of_interest %in% cause_codes) {
      stop("Requested event code ", event_of_interest,
           " not present in fitted model. Available: ",
           paste(cause_codes, collapse = ", "))
    }
  }

  newdata_prepared <- newdata[, modelout$expvars, drop = FALSE]
  for (var in modelout$expvars) {
    if (var %in% names(modelout$factor_levels)) {
      newdata_prepared[[var]] <- factor(
        newdata_prepared[[var]],
        levels = modelout$factor_levels[[var]]
      )
      unseen <- setdiff(unique(newdata[[var]]), modelout$factor_levels[[var]])
      if (length(unseen) > 0) {
        warning("Predictor '", var, "' has unseen levels: ",
                paste(unseen, collapse = ", "), ".")
      }
    }
  }

  if (anyNA(newdata_prepared)) {
    stop("Missing values detected in predictors after level alignment; impute before prediction.")
  }

  prep_new <- ttah_prepare_features(
    data = newdata_prepared,
    expvars = modelout$expvars,
    factor_levels = modelout$factor_levels,
    basis_specs = modelout$basis_specs
  )
  phi_new <- prep_new$phi

  latent_projection <- modelout$latent_projection
  phi_latent_new <- if (!is.null(latent_projection) && ncol(latent_projection) > 0) {
    tmp <- phi_new %*% latent_projection
    colnames(tmp) <- colnames(latent_projection)
    tmp
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

  pred_design <- ttah_build_long_design(
    phi = phi_new,
    phi_latent = phi_latent_new,
    time_basis = time_basis,
    interval_index = rep(K, n_obs),
    event = rep(0, n_obs)
  )

  design_df <- as.data.frame(pred_design$X)
  feature_cols <- modelout$model$feature_columns
  missing_cols <- setdiff(feature_cols, colnames(design_df))
  if (length(missing_cols) > 0) {
    zeros <- matrix(0, nrow = nrow(design_df), ncol = length(missing_cols))
    colnames(zeros) <- missing_cols
    design_df <- cbind(design_df, zeros)
  }
  design_df <- design_df[, feature_cols, drop = FALSE]

  probs <- predict(modelout$model$multinom, newdata = design_df, type = "probs")
  prob_cols <- setdiff(modelout$model$cause_levels, "no_event")
  if (is.null(dim(probs))) {
    probs <- matrix(probs, ncol = length(prob_cols))
  }
  if (is.null(colnames(probs))) {
    if (ncol(probs) != length(prob_cols)) {
      stop("Unable to align probability columns returned by multinom().")
    }
    colnames(probs) <- prob_cols
  }
  prob_matrix <- matrix(0, nrow = nrow(design_df), ncol = length(modelout$model$cause_levels))
  colnames(prob_matrix) <- modelout$model$cause_levels

  common_cols <- intersect(colnames(probs), colnames(prob_matrix))
  prob_matrix[, common_cols] <- probs[, common_cols, drop = FALSE]
  non_baseline <- setdiff(modelout$model$cause_levels, "no_event")
  prob_matrix[, "no_event"] <- pmax(
    1 - rowSums(prob_matrix[, non_baseline, drop = FALSE]),
    1e-8
  )

  hazards_array <- array(
    0,
    dim = c(K, n_obs, length(cause_codes)),
    dimnames = list(
      interval = seq_len(K),
      observation = seq_len(n_obs),
      cause = cause_codes
    )
  )
  stay_matrix <- matrix(0, nrow = K, ncol = n_obs)

  for (row_idx in seq_len(nrow(prob_matrix))) {
    obs_id <- pred_design$obs[row_idx]
    time_id <- pred_design$time[row_idx]
    stay_matrix[time_id, obs_id] <- prob_matrix[row_idx, "no_event"]
    for (cause in cause_codes) {
      hazards_array[time_id, obs_id, cause] <- prob_matrix[row_idx, paste0("cause_", cause)]
    }
  }

  Times <- c(0, time_grid)
  TotalSurvival <- matrix(1, nrow = K + 1, ncol = n_obs)
  CauseCIF_array <- array(
    0,
    dim = c(K + 1, n_obs, length(cause_codes)),
    dimnames = list(
      time = Times,
      observation = seq_len(n_obs),
      cause = cause_codes
    )
  )

  for (j in seq_len(n_obs)) {
    for (k in seq_len(K)) {
      hazard_vec <- hazards_array[k, j, ]
      CauseCIF_array[k + 1, j, ] <- CauseCIF_array[k, j, ] + TotalSurvival[k, j] * hazard_vec
      TotalSurvival[k + 1, j] <- TotalSurvival[k, j] * stay_matrix[k, j]
    }
  }

  event_idx <- match(event_of_interest, cause_codes)
  CauseSpecificHazard <- hazards_array[, , event_idx, drop = TRUE]
  CauseSpecificCIF <- CauseCIF_array[, , event_idx, drop = TRUE]

  if (!is.matrix(CauseSpecificHazard)) {
    CauseSpecificHazard <- matrix(CauseSpecificHazard, nrow = K, ncol = n_obs)
  }
  if (!is.matrix(CauseSpecificCIF)) {
    CauseSpecificCIF <- matrix(CauseSpecificCIF, nrow = K + 1, ncol = n_obs)
  }

  CauseSpecificCumHaz <- rbind(
    rep(0, n_obs),
    apply(CauseSpecificHazard, 2, cumsum)
  )
  CauseSpecificSurvival <- rbind(
    rep(1, n_obs),
    apply(CauseSpecificHazard, 2, function(col) cumprod(1 - col))
  )

  output <- list(
    CauseSpecificHazard = CauseSpecificHazard,
    CauseSpecificCIF = CauseSpecificCIF,
    CauseSpecificCumHaz = CauseSpecificCumHaz,
    CauseSpecificSurvival = CauseSpecificSurvival,
    TotalSurvival = TotalSurvival,
    StayProbability = stay_matrix,
    AllCauseHazard = hazards_array,
    AllCauseCIF = CauseCIF_array,
    Times = Times
  )

  if (!is.null(new_times)) {
    SortedTimes <- sort(unique(new_times))
    cif_interp_list <- lapply(seq_len(dim(CauseCIF_array)[3]), function(idx) {
      cifMatInterpolator(
        probsMat = CauseCIF_array[, , idx],
        times = Times,
        new_times = SortedTimes
      )
    })

    survival_interp <- survprobMatInterpolator(
      probsMat = TotalSurvival,
      times = Times,
      new_times = SortedTimes
    )

    hazard_recalc <- array(
      0,
      dim = c(length(SortedTimes) - 1, n_obs, length(cause_codes)),
      dimnames = list(
        interval = seq_len(length(SortedTimes) - 1),
        observation = seq_len(n_obs),
        cause = cause_codes
      )
    )
    stay_interp <- matrix(0, nrow = length(SortedTimes) - 1, ncol = n_obs)

    for (cause_idx in seq_along(cause_codes)) {
      cif_mat <- cif_interp_list[[cause_idx]]
      for (j in seq_len(n_obs)) {
        surv_prev <- survival_interp[-nrow(survival_interp), j]
        surv_prev <- pmax(surv_prev, 1e-8)
        diff_cif <- diff(cif_mat[, j])
        hazard_recalc[, j, cause_idx] <- pmax(diff_cif / surv_prev, 0)
      }
    }

    stay_interp <- pmax(
      survival_interp[-1, , drop = FALSE] / survival_interp[-nrow(survival_interp), , drop = FALSE],
      0
    )

    output$CauseSpecificHazard <- hazard_recalc[, , event_idx, drop = TRUE]
    if (!is.matrix(output$CauseSpecificHazard)) {
      output$CauseSpecificHazard <- matrix(output$CauseSpecificHazard,
                                           nrow = length(SortedTimes) - 1,
                                           ncol = n_obs)
    }
    output$CauseSpecificCIF <- cif_interp_list[[event_idx]]
    output$CauseSpecificCumHaz <- rbind(
      rep(0, n_obs),
      apply(output$CauseSpecificHazard, 2, cumsum)
    )
    output$CauseSpecificSurvival <- rbind(
      rep(1, n_obs),
      apply(output$CauseSpecificHazard, 2, function(col) cumprod(1 - col))
    )
    output$TotalSurvival <- survival_interp
    output$StayProbability <- stay_interp
    output$AllCauseHazard <- hazard_recalc
    output$AllCauseCIF <- array(
      NA_real_,
      dim = c(length(SortedTimes), n_obs, length(cause_codes)),
      dimnames = list(
        time = SortedTimes,
        observation = seq_len(n_obs),
        cause = cause_codes
      )
    )
    for (cause_idx in seq_along(cause_codes)) {
      output$AllCauseCIF[, , cause_idx] <- cif_interp_list[[cause_idx]]
    }

    output$Times <- SortedTimes
  }

  output
}
