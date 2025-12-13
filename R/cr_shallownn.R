# --- Internal Helper Functions for Custom NNet Training ---
initialize_weights <- function(n_in, n_hidden) {
  # Use a common initialization scheme (e.g., random uniform)
  # The nnet package uses a range of [-0.7, 0.7]
  rand_range <- 0.7
  W1 <- matrix(runif(n_in * n_hidden, -rand_range, rand_range), nrow = n_in, ncol = n_hidden)
  b1 <- runif(n_hidden, -rand_range, rand_range)
  W2 <- matrix(runif(n_hidden, -rand_range, rand_range), nrow = n_hidden, ncol = 1)
  b2 <- runif(1, -rand_range, rand_range)

  list(W1 = W1, b1 = b1, W2 = W2, b2 = b2)
}

forward_pass <- function(X, weights) {
  # Sigmoid activation function (as in nnet package)
  sigmoid <- function(z) 1 / (1 + exp(-z))

  # Z1 = X %*% W1 + b1 (broadcast b1)
  Z1 <- sweep(X %*% weights$W1, 2, weights$b1, "+")
  A1 <- sigmoid(Z1)
  # Z2 = A1 %*% W2 + b2 (linear output)
  Z2 <- sweep(A1 %*% weights$W2, 2, weights$b2, "+")

  list(output = Z2, hidden_activations = A1)
}

prepare_newdata_for_model <- function(modelout, newdata) {
  prepared <- newdata[, modelout$expvars, drop = FALSE]

  for (var in modelout$expvars) {
    if (var %in% names(modelout$factor_levels)) {
      training_levels <- modelout$factor_levels[[var]]
      if (is.factor(prepared[[var]])) {
        prepared[[var]] <- factor(prepared[[var]], levels = training_levels)
      } else if (is.character(prepared[[var]])) {
        prepared[[var]] <- factor(prepared[[var]], levels = training_levels)
      } else {
        prepared[[var]] <- factor(prepared[[var]], levels = training_levels)
      }
    }
  }

  numeric_vars <- names(modelout$model$scaling_params$mean)
  if (length(numeric_vars) > 0) {
    prepared[, numeric_vars] <- sweep(prepared[, numeric_vars, drop = FALSE],
                                      2, modelout$model$scaling_params$mean, "-")
    prepared[, numeric_vars] <- sweep(prepared[, numeric_vars, drop = FALSE],
                                      2, modelout$model$scaling_params$sd, "/")
  }

  x_matrix <- stats::model.matrix(~ . - 1, data = prepared)

  list(data = prepared, x = x_matrix)
}

baseline_cumhaz_stepfun <- function(baseline_df) {
  stats::stepfun(baseline_df$time, c(0, baseline_df$cumhaz))
}

hazard_increment_matrix <- function(step_fun, times, risk_scores) {
  if (!is.numeric(times) || any(is.na(times))) {
    stop("`times` must be a numeric vector without NA values")
  }
  if (length(times) < 2) {
    stop("`times` must contain at least two time points (including time 0)")
  }

  baseline_cumhaz <- step_fun(times)
  baseline_increments <- diff(baseline_cumhaz)

  if (any(baseline_increments < 0)) {
    stop("Baseline cumulative hazard must be non-decreasing")
  }

  outer(baseline_increments, risk_scores)
}

predict_risk_scores <- function(modelout, x_matrix) {
  log_risk <- forward_pass(x_matrix, modelout$model$weights)$output
  log_risk <- as.vector(log_risk)
  log_risk <- pmax(pmin(log_risk, 50), -50)
  exp(log_risk)
}

unpack_weights <- function(w_vec, n_in, n_hidden) {
  W1_end <- n_in * n_hidden
  b1_end <- W1_end + n_hidden
  W2_end <- b1_end + n_hidden
  b2_end <- W2_end + 1

  list(
    W1 = matrix(w_vec[1:W1_end], nrow = n_in, ncol = n_hidden),
    b1 = w_vec[(W1_end + 1):b1_end],
    W2 = matrix(w_vec[(b1_end + 1):W2_end], nrow = n_hidden, ncol = 1),
    b2 = w_vec[b2_end]
  )
}

#'
#'
#'
#'
CRModel_ShallowNN <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                              size = 5, decay = 0.01, maxit = 1000, verbose = FALSE) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(data)) {
    stop("Input 'data' is missing")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("`expvars` must be a non-empty character vector")
  }
  if (!timevar %in% colnames(data)) {
    stop(paste0("`timevar` not found in data: ", timevar))
  }
  if (!eventvar %in% colnames(data)) {
    stop(paste0("`eventvar` not found in data: ", eventvar))
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste0("The following `expvars` not found in data: ", paste(missing_vars, collapse=", ")))
  }
  if (!is.null(event_codes) && length(event_codes) == 0) {
    stop("`event_codes` must be NULL or a non-empty vector")
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================
  if (verbose) cat("Creating variable profile...\n")
  varprof <- VariableProfile(data, expvars)

  # Store factor levels for prediction
  factor_levels <- list()
  for (var in expvars) {
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      factor_levels[[var]] <- levels(as.factor(data[[var]]))
    }
  }

  # Ensure event variable is numeric
  data[[eventvar]] <- as.numeric(data[[eventvar]])

  # Prepare data subset with complete cases
  XYTrain <- data[, c(timevar, eventvar, expvars), drop = FALSE]
  n_original <- nrow(XYTrain)
  XYTrain <- XYTrain[complete.cases(XYTrain), , drop = FALSE]
  n_complete <- nrow(XYTrain)

  if (n_complete < n_original) {
    warning(sprintf("Removed %d rows with missing values. %d rows remaining.",
                    n_original - n_complete, n_complete))
  }
  if (n_complete < 10) {
    stop("Insufficient data after removing missing values. Need at least 10 observations.")
  }

  available_events <- sort(unique(as.character(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0])))
  if (length(available_events) == 0) {
    stop("No events found in the training data.")
  }

  if (is.null(event_codes)) {
    event_codes <- available_events[1]
  }
  event_codes <- as.character(event_codes)
  if (length(event_codes) != 1) {
    stop("`event_codes` must be a single value (one event code)")
  }
  if (!event_codes %in% available_events) {
    stop(paste0("`event_codes` ", event_codes, " not present in training data. No events of type ", event_codes, "."))
  }
  event_codes_numeric <- suppressWarnings(as.numeric(event_codes))
  if (is.na(event_codes_numeric)) {
    stop(paste0("`event_codes` must be numeric or coercible to numeric. Unable to coerce '", event_codes, "' to numeric."))
  }
  primary_event_code <- event_codes[1]
  primary_event_numeric <- event_codes_numeric[1]
  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == primary_event_numeric]
  if (length(event_times) == 0) {
    stop(paste0("No events of type ", primary_event_code, " in training data. Cannot fit competing risks model."))
  }
  times <- sort(unique(event_times))

  time_range <- range(c(0, XYTrain[[timevar]][XYTrain[[eventvar]] %in% primary_event_numeric]), na.rm = TRUE)

  # ============================================================================
  # Model Fitting (Custom Shallow Neural Network for Competing Risks)
  # ============================================================================
  if (verbose) cat("Fitting shallow neural network for competing risks (event type", primary_event_code, ")...\n")

  # Standardize numeric variables
  numeric_vars <- expvars[sapply(XYTrain[, expvars, drop = FALSE], is.numeric)]
  if (length(numeric_vars) > 0) {
    scaling_params <- list(
      mean = sapply(XYTrain[, numeric_vars, drop = FALSE], mean, na.rm = TRUE),
      sd = sapply(XYTrain[, numeric_vars, drop = FALSE], sd, na.rm = TRUE)
    )
  } else {
    scaling_params <- list(mean = numeric(0), sd = numeric(0))
  }

  # Apply scaling
  XYTrain_scaled <- XYTrain
  if (length(numeric_vars) > 0) {
    XYTrain_scaled[, numeric_vars] <- scale(XYTrain[, numeric_vars, drop = FALSE])
  }

  # Create model matrix
  x_train <- stats::model.matrix(~ . - 1, data = XYTrain_scaled[, expvars, drop = FALSE])
  n_features <- ncol(x_train)

  # Sort data by descending time for efficient loss calculation
  sort_order <- order(XYTrain[[timevar]], -XYTrain[[eventvar]])
  x_train <- x_train[sort_order, ]
  time_sorted <- XYTrain[[timevar]][sort_order]
  event_sorted <- XYTrain[[eventvar]][sort_order]

  # Create cause-specific indicators for Fine-Gray style loss
  # For competing risks, we use a modified loss that accounts for censoring and competing events
  status_fg <- ifelse(event_sorted == 0, 0,  # censored
                      ifelse(event_sorted == primary_event_numeric, 1, 2))  # event of interest or competing

  # ============================================================================
  # Neural Network Training with Fine-Gray Style Loss
  # ============================================================================
  # Initialize weights
  initial_w_vec <- unlist(initialize_weights(n_features, size))

  # Define objective function (negative log-likelihood for Fine-Gray + L2 penalty)
  objective_fn <- function(w_vec) {
    weights <- unpack_weights(w_vec, n_features, size)

    # Forward pass to get log-risk (subdistribution hazard)
    log_risk <- forward_pass(x_train, weights)$output

    # Fine-Gray style loss calculation for competing risks
    # This implements a simplified version of the Fine-Gray likelihood
    risk_scores <- exp(log_risk)

    # Calculate loss similar to Fine-Gray model
    # For each event time, compute the contribution
    loss <- 0
    n_events <- 0

    for (i in seq_along(time_sorted)) {
      if (status_fg[i] == 1) {  # Event of interest
        n_events <- n_events + 1

        # Risk set at this time (all subjects still at risk)
        risk_set <- risk_scores[i:length(risk_scores)]

        # Contribution to log-likelihood
        # log(risk_i / sum(risk_set))
        if (sum(risk_set) > 0) {
          loss <- loss + log_risk[i] - log(sum(risk_set))
        }
      }
    }

    # For censored and competing events, they contribute 0 to the likelihood
    # (they're handled implicitly through the risk set calculations)

    # L2 penalty
    l2_penalty <- decay * sum(w_vec^2)

    # Return negative log-likelihood + penalty
    -loss + l2_penalty
  }

  # Define gradient function
  gradient_fn <- function(w_vec) {
    weights <- unpack_weights(w_vec, n_features, size)

    # Forward pass
    pass <- forward_pass(x_train, weights)
    log_risk <- pass$output
    A1 <- pass$hidden_activations

    # Backward pass for Fine-Gray style loss
    risk_scores <- exp(log_risk)

    d_log_risk <- rep(0, length(log_risk))

    for (i in seq_along(time_sorted)) {
      if (status_fg[i] == 1) {  # Event of interest
        # Risk set at this time
        risk_set <- risk_scores[i:length(risk_scores)]
        risk_set_sum <- sum(risk_set)

        if (risk_set_sum > 0) {
          # Gradient contribution
          d_log_risk[i] <- d_log_risk[i] + 1  # From log_risk[i] term
          d_log_risk[i:length(log_risk)] <- d_log_risk[i:length(log_risk)] -
            risk_scores[i:length(log_risk)] / risk_set_sum  # From -log(sum) term
        }
      }
    }

    # Backpropagate through output layer
    dZ2 <- d_log_risk
    dW2 <- t(A1) %*% dZ2
    db2 <- sum(dZ2)
    dA1 <- dZ2 %*% t(weights$W2)

    # Backpropagate through sigmoid
    dZ1 <- dA1 * (A1 * (1 - A1))
    dW1 <- t(x_train) %*% dZ1
    db1 <- colSums(dZ1)

    # Add L2 gradient
    l2_grad <- 2 * decay * w_vec

    # Pack gradients
    c(as.vector(dW1), db1, as.vector(dW2), db2) + l2_grad
  }

  # Optimize
  if (verbose) cat("Optimizing neural network parameters...\n")
  opt_result <- optim(
    par = initial_w_vec,
    fn = objective_fn,
    gr = gradient_fn,
    method = "BFGS",
    control = list(maxit = maxit, trace = as.numeric(verbose))
  )

  final_weights <- unpack_weights(opt_result$par, n_features, size)

  # ============================================================================
  # Calculate Baseline Subdistribution Hazard (Fine-Gray style)
  # ============================================================================
  final_log_risk <- forward_pass(x_train, final_weights)$output
  final_risk_scores <- exp(final_log_risk)

  # Calculate baseline subdistribution hazard similar to Fine-Gray
  baseline_hazard_increment <- sapply(times, function(t_i) {
    event_indices <- which(time_sorted == t_i & status_fg == 1)
    if (length(event_indices) == 0) {
      return(0)
    }

    risk_set_start <- min(event_indices)
    risk_set <- final_risk_scores[risk_set_start:length(final_risk_scores)]
    if (sum(risk_set) == 0) {
      return(0)
    }
    sum(final_risk_scores[event_indices]) / sum(risk_set)
  })

  cumulative_baseline_cumhaz <- cumsum(baseline_hazard_increment)

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Shallow neural network competing risks model fitting complete.\n")

  # Create a model object that mimics nnet structure for compatibility
  model <- list(
    weights = final_weights,
    baseline_hazard = data.frame(time = times, hazard = baseline_hazard_increment),
    baseline_cumhaz = data.frame(time = times, cumhaz = cumulative_baseline_cumhaz),
    scaling_params = scaling_params,
    n = c(n_features, size, 1),  # Mimic nnet structure
    call = match.call()
  )
  class(model) <- "nnet"

  result <- list(
    model = model,
    times = times,
    varprof = varprof,
    expvars = expvars,
    factor_levels = factor_levels,
    event_codes = primary_event_code,
    event_codes_numeric = event_codes_numeric,
    default_event_code = primary_event_code,
    model_type = "cr_shallownn",
    timevar = timevar,
    eventvar = eventvar,
    time_range = time_range
  )

  class(result) <- "ml4t2e_cr_shallownn"
  return(result)
}

#'
#'
#'
#'
Predict_CRModel_ShallowNN <- function(modelout, newdata, new_times = NULL,
                                      event_of_interest = NULL, other_models = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("`modelout` is missing")
  }
  if (!is.list(modelout) || !all(c("expvars", "default_event_code") %in% names(modelout))) {
    stop("Input 'modelout' must be output from CRModel_ShallowNN")
  }
  if (missing(newdata)) {
    stop("`newdata` is missing")
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame")
  }
  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop(paste0("The following variables are missing in `newdata`: ", paste(missing_vars, collapse = ", ")))
  }
  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    event_of_interest <- modelout$default_event_code
  }
  event_of_interest <- as.character(event_of_interest)
  if (length(event_of_interest) != 1) {
    stop("`event_of_interest` must be a single value (one event code)")
  }
  if (!identical(event_of_interest, modelout$default_event_code)) {
    stop(paste0("Shallow neural network models can only predict for the event they were trained on (event code = ",
                modelout$default_event_code, "). Requested event code: ", event_of_interest))
  }
  # Validate new_times if provided
  if (!is.null(new_times)) {
    if (!is.numeric(new_times)) {
      stop("`new_times` must be numeric")
    }
    if (any(new_times < 0)) {
      stop("`new_times` must be a numeric vector of non-negative values")
    }
  }
  
  if (!("baseline_cumhaz" %in% names(modelout$model))) {
    if ("baseline_subhaz" %in% names(modelout$model)) {
      legacy_df <- modelout$model$baseline_subhaz
      if (!all(c("time", "subhaz") %in% colnames(legacy_df))) {
        stop("Legacy model object lacks required baseline subhazard columns `time` and `subhaz`.")
      }
      cumhaz_vals <- legacy_df$subhaz
      modelout$model$baseline_cumhaz <- data.frame(time = legacy_df$time, cumhaz = cumhaz_vals)
      hazard_increments <- c(cumhaz_vals[1], diff(cumhaz_vals))
      modelout$model$baseline_hazard <- data.frame(time = legacy_df$time, hazard = hazard_increments)
    } else {
      stop(
        paste(
          "Fitted model object is missing `baseline_cumhaz`;",
          "re-fit the model using the updated ml4time2event version."
        )
      )
    }
  }

  # =========================================================================
  # Prepare newdata for primary model
  # =========================================================================
  prepared_primary <- prepare_newdata_for_model(modelout, newdata)
  x_primary <- prepared_primary$x

  if (ncol(x_primary) != modelout$model$n[1]) {
    stop("Model matrix column count does not match stored model configuration")
  }

  risk_scores_primary <- predict_risk_scores(modelout, x_primary)

  baseline_step_primary <- baseline_cumhaz_stepfun(modelout$model$baseline_cumhaz)
  primary_time_grid <- sort(unique(c(0, modelout$model$baseline_cumhaz$time)))

  # =========================================================================
  # Handle additional cause-specific models (competing events)
  # =========================================================================
  other_model_list <- list()
  if (!is.null(other_models)) {
    if (!is.list(other_models)) {
      stop("`other_models` must be a list of fitted CRModel_ShallowNN objects")
    }

    other_model_list <- other_models

    observed_event_codes <- modelout$default_event_code
    for (i in seq_along(other_model_list)) {
      model_i <- other_model_list[[i]]
      if (!inherits(model_i, "ml4t2e_cr_shallownn")) {
        stop("All elements of `other_models` must be outputs from CRModel_ShallowNN")
      }
      if (identical(model_i$default_event_code, modelout$default_event_code)) {
        stop("`other_models` contains a model for the same event code as the primary model")
      }
      if (model_i$default_event_code %in% observed_event_codes) {
        stop("Duplicate event codes detected across competing models. Ensure one model per event code.")
      }
      observed_event_codes <- c(observed_event_codes, model_i$default_event_code)
    }
  }

  # =========================================================================
  # Construct common time grid across all cause-specific models
  # =========================================================================
  additional_times <- unlist(lapply(other_model_list, function(mod) mod$model$baseline_cumhaz$time))
  full_time_grid <- sort(unique(c(primary_time_grid, additional_times)))
  if (length(full_time_grid) == 0) {
    full_time_grid <- primary_time_grid
  }

  if (!any(full_time_grid == 0)) {
    full_time_grid <- c(0, full_time_grid)
  }

  # =========================================================================
  # Compute hazard increments for primary event
  # =========================================================================
  hazard_primary <- hazard_increment_matrix(
    baseline_step_primary,
    times = full_time_grid,
    risk_scores = risk_scores_primary
  )

  # =========================================================================
  # Compute hazard increments for competing events
  # =========================================================================
  competing_hazards <- list()
  if (length(other_model_list) > 0) {
    for (i in seq_along(other_model_list)) {
      model_i <- other_model_list[[i]]
      model_name <- names(other_model_list)[i]
      if (is.null(model_name) || identical(model_name, "")) {
        model_name <- paste0("event_", model_i$default_event_code)
      }

      missing_vars_i <- setdiff(model_i$expvars, colnames(newdata))
      if (length(missing_vars_i) > 0) {
        stop(
          paste0(
            "The following variables required by competing model `",
            model_name,
            "` are missing in `newdata`: ",
            paste(missing_vars_i, collapse = ", ")
          )
        )
      }

      prepared_i <- prepare_newdata_for_model(model_i, newdata)
      risk_i <- predict_risk_scores(model_i, prepared_i$x)
      step_i <- baseline_cumhaz_stepfun(model_i$model$baseline_cumhaz)
      competing_hazards[[model_name]] <- hazard_increment_matrix(
        step_fun = step_i,
        times = full_time_grid,
        risk_scores = risk_i
      )
    }
  }

  # =========================================================================
  # Assemble cumulative hazards and survival on the native grid
  # =========================================================================
  intervals <- nrow(hazard_primary)
  n_obs <- ncol(hazard_primary)
  times_full <- full_time_grid

  cumulative_hazard_primary <- rbind(
    rep(0, n_obs),
    apply(hazard_primary, 2, cumsum)
  )

  total_hazard_matrix <- hazard_primary
  if (length(competing_hazards) > 0) {
    for (comp_mat in competing_hazards) {
      total_hazard_matrix <- total_hazard_matrix + comp_mat
    }
  }

  cumulative_total_hazard <- rbind(
    rep(0, n_obs),
    apply(total_hazard_matrix, 2, cumsum)
  )

  survival_primary <- exp(-cumulative_hazard_primary)
  survival_total <- exp(-cumulative_total_hazard)

  # =========================================================================
  # Compute CIF using Aalen-Johansen when competing models are provided
  # =========================================================================
  cif_matrix <- NULL
  if (length(competing_hazards) > 0) {
    cif_matrix <- matrix(0, nrow = intervals + 1, ncol = n_obs)
    for (k in seq_len(intervals)) {
      cif_matrix[k + 1, ] <- cif_matrix[k, ] + survival_total[k, ] * hazard_primary[k, ]
    }
    cif_matrix <- pmin(pmax(cif_matrix, 0), 1)
  }

  # =========================================================================
  # Optional interpolation to new times
  # =========================================================================
  if (!is.null(new_times)) {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    target_times <- sort(unique(new_times))
    if (!any(target_times == 0)) {
      target_times <- c(0, target_times)
    }

    step_interpolate_matrix <- function(mat, src_times, tgt_times) {
      res <- apply(mat, 2, function(col_vals) {
        stats::approx(
          x = src_times,
          y = col_vals,
          xout = tgt_times,
          method = "constant",
          rule = 2,
          f = 0
        )$y
      })
      matrix(res, nrow = length(tgt_times), ncol = ncol(mat))
    }

    cumulative_hazard_primary <- step_interpolate_matrix(
      cumulative_hazard_primary,
      src_times = times_full,
      tgt_times = target_times
    )

    survival_primary <- exp(-cumulative_hazard_primary)

    cumulative_total_hazard <- step_interpolate_matrix(
      cumulative_total_hazard,
      src_times = times_full,
      tgt_times = target_times
    )
    survival_total <- exp(-cumulative_total_hazard)

    hazard_primary <- apply(cumulative_hazard_primary, 2, function(col) c(0, diff(col)))
    hazard_primary <- matrix(
      hazard_primary,
      nrow = nrow(cumulative_hazard_primary),
      ncol = ncol(cumulative_hazard_primary)
    )

    if (!is.null(cif_matrix)) {
      cif_matrix <- cifMatInterpolator(
        probsMat = cif_matrix,
        times = times_full,
        new_times = target_times
      )
    }

    times_full <- target_times
  } else {
    hazard_primary <- rbind(rep(0, n_obs), hazard_primary)
  }

  # =========================================================================
  # Assemble results
  # =========================================================================
  result <- list(
    CIFs = cif_matrix,
    Times = times_full,
    CauseSpecificHazard = hazard_primary,
    CauseSpecificCumHaz = cumulative_hazard_primary,
    CauseSpecificSurvival = survival_primary,
    TotalSurvival = survival_total
  )

  if (is.null(cif_matrix)) {
    warning(
      paste(
        "Cumulative incidence functions require cause-specific models for all",
        "competing events. Provide them via `other_models` to obtain CIF estimates."
      )
    )
  }

  result
}
