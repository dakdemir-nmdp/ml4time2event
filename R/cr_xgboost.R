#'
#'
#'
#'
CRModel_xgboost <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                           eta = 0.01, max_depth = 5, nrounds = 100,
                           ntimes = 50, verbose = FALSE, ...) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("Input 'expvars' must be a non-empty character vector.")
  }
  if (!timevar %in% colnames(data)) {
    stop("Input 'timevar' not found in data: ", timevar)
  }
  if (!eventvar %in% colnames(data)) {
    stop("Input 'eventvar' not found in data: ", eventvar)
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following 'expvars' not found in data: ", paste(missing_vars, collapse=", "))
  }
  if (!is.null(event_codes) && length(event_codes) == 0) {
    stop("Input 'event_codes' must be NULL or a non-empty vector.")
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================
  if (verbose) cat("Creating variable profile...\n")
  varprof <- VariableProfile(data, expvars)

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

  available_events <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))
  if (length(available_events) == 0) {
    stop("No events found in the training data.")
  }

  if (is.null(event_codes)) {
    event_codes <- as.character(available_events)
  }

  event_codes <- as.character(event_codes)
  missing_event_codes <- setdiff(event_codes, as.character(available_events))
  if (length(missing_event_codes) > 0) {
    stop("The following event_codes are not present in the data: ",
         paste(missing_event_codes, collapse = ", "))
  }

  event_codes_numeric <- suppressWarnings(as.numeric(event_codes))
  if (any(is.na(event_codes_numeric))) {
    stop("XGBoost competing risks requires numeric event codes. Unable to coerce: ",
         paste(event_codes[is.na(event_codes_numeric)], collapse = ", "))
  }

  primary_event_code <- event_codes[1]
  primary_event_numeric <- event_codes_numeric[1]

  # Get unique event times for the default event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == primary_event_numeric]
  if (length(event_times) == 0) {
    stop("No events of type ", primary_event_code, " in training data. Cannot fit competing risks model.")
  }

  # Store event time range for reference
  time_range <- range(XYTrain[[timevar]][XYTrain[[eventvar]] %in% event_codes_numeric], na.rm = TRUE)

  # ============================================================================
  # Model Fitting - Cause-Specific XGBoost Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  if (verbose) {
    cat("Fitting cause-specific XGBoost models for event codes:", paste(event_codes, collapse = ", "), "\n")
  }

  # Store models for all event types in the requested set
  xgb_models_all_causes <- vector("list", length(event_codes))
  names(xgb_models_all_causes) <- event_codes
  feature_names <- NULL

  # Fit a separate XGBoost model for each event type
  for (idx in seq_along(event_codes_numeric)) {
    cause <- event_codes_numeric[idx]
    cause_char <- event_codes[idx]
    if (verbose) cat("Fitting XGBoost model for event type", cause, "...\n")

    # Create cause-specific data: event = 1 if this cause, 0 otherwise (censored or competing)
    XYTrain_cause <- XYTrain
    XYTrain_cause$status_cs <- ifelse(XYTrain[[eventvar]] == cause, 1, 0)

    if (sum(XYTrain_cause$status_cs) < 5) {
      warning("Fewer than 5 events of type ", cause_char, ". Skipping this cause.")
      xgb_models_all_causes[[cause_char]] <- NULL
      next
    }

    # ============================================================================
    # Prepare Data Matrix for XGBoost
    # ============================================================================
    # Prepare data matrix (handle factors using model.matrix)
    X <- stats::model.matrix(~-1+., XYTrain_cause[, expvars, drop = FALSE])
    if (is.null(feature_names)) {
      feature_names <- colnames(X)  # Store feature names from first model
    }

    # Prepare labels for XGBoost: negative time for censored, positive time for event
    ytimes <- XYTrain_cause[[timevar]]
    yevents <- XYTrain_cause$status_cs
    yTrain <- ytimes
    yTrain[yevents == 0] <- (-yTrain[yevents == 0])

    # ============================================================================
    # XGBoost Parameters
    # ============================================================================
    params <- list(
      objective = 'survival:cox',
      eval_metric = "cox-nloglik",
      eta = eta,
      max_depth = max_depth,
      ...
    )

    # ============================================================================
    # Train XGBoost Model
    # ============================================================================
    # Create DMatrix
    dtrain <- xgboost::xgb.DMatrix(data = X, label = yTrain)

    # Train the model
    xgb_model_cause <- tryCatch(
      xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = ifelse(verbose, 1, 0)
      ),
      error = function(e) {
        warning("Failed to fit XGBoost model for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(xgb_model_cause)) {
      next
    }

    # ============================================================================
    # Create Baseline Model for Prediction
    # ============================================================================
    # Get linear predictors for training data
    train_linear_preds <- predict(xgb_model_cause, X)

    # Create survival data for the cause-specific case
    train_survival_data <- data.frame(
      time = XYTrain_cause[[timevar]],
      event = XYTrain_cause$status_cs
    )

    # Use score2proba to create baseline hazard model
    baseline_info <- score2proba(
      datasurv = train_survival_data,
      score = train_linear_preds,
      conf.int = 0.95,
      which.est = "point"
    )

    # Store baseline model in a wrapper list
    xgb_wrapper <- list(
      model = xgb_model_cause,
      baseline_model = baseline_info$model,
      baseline_sf = baseline_info$sf
    )

    xgb_models_all_causes[[cause_char]] <- xgb_wrapper
  }

  # The main model for the event of interest
  xgb_model <- xgb_models_all_causes[[primary_event_code]]

  if (is.null(xgb_model)) {
    stop("Failed to fit XGBoost model for the event of interest (event code = ", primary_event_code, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific XGBoost model fitting complete.\n")

  result <- list(
    xgb_model = xgb_model,
    xgb_models_all_causes = xgb_models_all_causes,  # All cause-specific models for Aalen-Johansen
    event_codes = event_codes,
    event_codes_numeric = event_codes_numeric,
    default_event_code = primary_event_code,
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_xgboost",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    time_range = time_range,
    feature_names = feature_names
  )

  class(result) <- c("ml4t2e_cr_xgboost", "CRModel_xgboost")
  return(result)
}

#'
#'
#'
#'
Predict_CRModel_xgboost <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("'modelout' is missing")
  }
  if (!is.list(modelout) || !all(c("expvars", "default_event_code") %in% names(modelout))) {
    stop("must be output from CRModel_xgboost")
  }
  if (missing(newdata)) {
    stop("'newdata' is missing")
  }
  if (!is.data.frame(newdata)) {
    stop("Input 'newdata' must be a data frame.")
  }

  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop("The following 'expvars' are missing in 'newdata': ", paste(missing_vars, collapse = ", "))
  }

  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    event_of_interest <- modelout$default_event_code
  }

  event_of_interest <- as.character(event_of_interest)

  if (length(event_of_interest) != 1) {
    stop("Input 'event_of_interest' must be a single event code.")
  }

  if (!event_of_interest %in% modelout$event_codes) {
    stop("Input 'event_of_interest' was not present in training data. Available event codes: ",
         paste(modelout$event_codes, collapse = ", "))
  }

  # Generate default times if not specified
  if (is.null(new_times)) {
    # Create a reasonable default grid: include 0 and event times
    target_times <- sort(unique(c(0, modelout$times)))
  } else {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("Input 'new_times' must be a numeric vector of non-negative values.")
    }
    target_times <- sort(unique(new_times))
  }

  # ============================================================================
  # Prepare newdata
  # ============================================================================
  # Apply the same model.matrix transformation as training
  X_new <- stats::model.matrix(~-1+., newdata[, modelout$expvars, drop = FALSE])

  # Ensure X_new has the same columns as training (in case of missing factor levels)
  missing_cols <- setdiff(modelout$feature_names, colnames(X_new))
  if (length(missing_cols) > 0) {
    # Add missing columns with zeros
    for (col in missing_cols) {
      X_new <- cbind(X_new, 0)
      colnames(X_new)[ncol(X_new)] <- col
    }
  }

  # Ensure correct column order
  X_new <- X_new[, modelout$feature_names, drop = FALSE]

  # ============================================================================
  # Make Predictions using Aalen-Johansen Estimator
  # ============================================================================
  base_times <- sort(unique(c(0, modelout$times)))
  if (length(base_times) == 0) {
    base_times <- sort(unique(c(0, modelout$time_range)))
  }

  if (length(target_times) == 0) {
    target_times <- base_times
  }

  n_times <- length(base_times)
  n_obs <- nrow(X_new)
  n_causes <- length(modelout$event_codes)

  cum_hazards_all_causes <- array(NA_real_, dim = c(n_times, n_obs, n_causes))
  missing_cause_idx <- integer(0)

  for (i in seq_len(n_causes)) {
    cause_char <- modelout$event_codes[i]
    cause_wrapper <- modelout$xgb_models_all_causes[[cause_char]]

    if (is.null(cause_wrapper)) {
      missing_cause_idx <- c(missing_cause_idx, i)
      next
    }
    
    cause_model <- cause_wrapper$model

    linear_preds <- tryCatch(
      predict(cause_model, X_new),
      error = function(e) {
        warning("Prediction failed for cause ", cause_char, ": ", e$message)
        rep(NA_real_, n_obs)
      }
    )

    sf_baseline <- tryCatch(
      survival::survfit(cause_wrapper$baseline_model,
                        newdata = data.frame("score" = 0)),
      error = function(e) NULL
    )

    if (is.null(sf_baseline) || length(sf_baseline$time) == 0 || any(is.na(sf_baseline$surv))) {
      missing_cause_idx <- c(missing_cause_idx, i)
      next
    }

    baseline_cum_hazard <- stats::approx(
      x = c(0, sf_baseline$time),
      y = c(0, -log(pmax(sf_baseline$surv, 1e-10))),
      xout = base_times,
      method = "constant",
      f = 0,
      rule = 2
    )$y

    for (j in seq_len(n_obs)) {
      if (is.na(linear_preds[j])) {
        cum_hazards_all_causes[, j, i] <- NA_real_
      } else {
        cum_hazards_all_causes[, j, i] <- baseline_cum_hazard * exp(linear_preds[j])
      }
    }
  }
  if (length(missing_cause_idx) > 0) {
    warning("No model available for event code(s): ", paste(modelout$event_codes[missing_cause_idx], collapse = ", "),
            ". Assuming zero hazard contribution for these causes.")
    for (idx in missing_cause_idx) {
      if (all(is.na(cum_hazards_all_causes[, , idx]))) {
        cum_hazards_all_causes[, , idx] <- 0
      }
    }
  }

  # Step 2: Calculate overall survival S(t) = exp(-Σ_j Λ_j(t))
  overall_cum_hazard <- apply(cum_hazards_all_causes, c(1, 2), function(x) {
    if (any(is.na(x))) {
      NA_real_
    } else {
      sum(x)
    }
  })
  overall_survival <- exp(-overall_cum_hazard)

  # Step 3: Calculate CIF using Aalen-Johansen formula for requested event
  target_idx <- match(event_of_interest, modelout$event_codes)
  target_cum_hazards <- cum_hazards_all_causes[, , target_idx]
  cif_matrix <- matrix(NA_real_, nrow = n_times, ncol = n_obs)

  for (j in seq_len(n_obs)) {
    if (all(is.na(target_cum_hazards[, j])) || all(is.na(overall_survival[, j]))) {
      next
    }

    hazard_increments <- c(target_cum_hazards[1, j], diff(target_cum_hazards[, j]))
    cif_vals <- rep(NA_real_, n_times)

    for (t in seq_len(n_times)) {
      current_increment <- hazard_increments[t]

      if (is.na(current_increment)) {
        cif_vals[t] <- NA_real_
        next
      }

      if (t == 1) {
        # At time 0, CIF should be 0
        if (abs(base_times[t]) < .Machine$double.eps) {
          cif_vals[t] <- 0
        } else {
          # First time point is not 0, use Aalen-Johansen formula
          cif_vals[t] <- 1.0 * current_increment
        }
      } else {
        prev_cif <- cif_vals[t - 1]
        prev_surv <- overall_survival[t - 1, j]

        if (is.na(prev_cif) || is.na(prev_surv)) {
          cif_vals[t] <- NA_real_
        } else {
          cif_vals[t] <- prev_cif + prev_surv * current_increment
        }
      }
    }

    cif_matrix[, j] <- cif_vals
  }

  # ============================================================================
  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  requested_times <- target_times
  if (length(requested_times) == length(base_times) && all(requested_times == base_times)) {
    result_cifs <- cif_matrix
    result_times <- base_times
  } else {
    result_cifs <- cifMatInterpolaltor(
      probsMat = cif_matrix,
      times = base_times,
      new_times = requested_times
    )
    result_times <- requested_times
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    CIFs = result_cifs,
    Times = result_times
  )

  return(result)
}
