#'
#'
#'
#'

CRModel_SurvReg <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                           dist = "exponential", ntimes = 50, verbose = FALSE, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
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


  # Handle event_of_interest for consistency
  if (is.null(event_codes)) {
    event_codes <- available_events
  }
  event_codes <- as.character(event_codes)
  missing_codes <- setdiff(event_codes, available_events)
  if (length(missing_codes) > 0) {
    stop(paste0("The following `event_codes` are not present in the data: ", paste(missing_codes, collapse = ", ")))
  }
  event_codes_numeric <- suppressWarnings(as.numeric(event_codes))
  if (any(is.na(event_codes_numeric))) {
    stop(paste0("`event_codes` must be numeric or coercible to numeric. Unable to coerce: ", paste(event_codes[is.na(event_codes_numeric)], collapse = ", ")))
  }

  # If event_of_interest is provided, prioritize it as the first event code
  if (!is.null(event_of_interest)) {
    event_of_interest <- as.character(event_of_interest)
    if (!event_of_interest %in% event_codes) {
      stop(paste0("`event_of_interest` '", event_of_interest, "' is not in `event_codes`."))
    }
    # Reorder event_codes so event_of_interest is first
    event_codes <- c(event_of_interest, setdiff(event_codes, event_of_interest))
    event_codes_numeric <- suppressWarnings(as.numeric(event_codes))
  }
  primary_event_code <- event_codes[1]
  primary_event_numeric <- event_codes_numeric[1]

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == primary_event_numeric]
  if (length(event_times) == 0) {
    stop("No events of type ", primary_event_code, " in training data. Cannot fit competing risks model.")
  }

  # Store event time range for reference
  time_range <- range(c(0, XYTrain[[timevar]][XYTrain[[eventvar]] %in% event_codes_numeric]), na.rm = TRUE)

  # ============================================================================
  # Model Fitting - Cause-Specific SurvReg Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  all_event_types <- event_codes_numeric

  if (verbose) {
    cat("Fitting cause-specific SurvReg models for all event types:", paste(all_event_types, collapse = ", "), "\n")
  }

  # Store models for all event types
  survreg_models_all_causes <- vector("list", length(all_event_types))
  names(survreg_models_all_causes) <- as.character(all_event_types)

  # Fit a separate SurvReg model for each event type
  for (cause in all_event_types) {
    if (verbose) cat("Fitting SurvReg model for event type", cause, "...\n")

    # Create cause-specific data: event = 1 if this cause, 0 otherwise (censored or competing)
    XYTrain_cause <- XYTrain
    XYTrain_cause$status_cs <- ifelse(XYTrain[[eventvar]] == cause, 1, 0)

    if (sum(XYTrain_cause$status_cs) < 5) {
      warning("Fewer than 5 events of type ", cause, ". Skipping this cause.")
      survreg_models_all_causes[[as.character(cause)]] <- NULL
      next
    }

    # ============================================================================
    # Forward Selection with AIC (adapted for cause-specific modeling)
    # ============================================================================
    if (verbose) cat("Performing forward selection with AIC for event", cause, "...\n")

    selected_vars <- c()
    candidate_vars <- expvars

    # Start with intercept-only model AIC
    null_formula <- stats::as.formula(paste("survival::Surv(", timevar, ", status_cs) ~ 1"))
    null_model <- tryCatch(
      survival::survreg(null_formula, data = XYTrain_cause, dist = dist, x = FALSE, y = FALSE),
      error = function(e) {
        if (verbose) warning("Failed to fit null model for cause ", cause, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(null_model)) {
      warning("Failed to fit intercept-only model for cause ", cause, ". Skipping.")
      next
    }
    
    best_aic <- stats::AIC(null_model)
    if (verbose) print(paste("Cause", cause, "- Initial AIC (Intercept only):", round(best_aic, 2)))

    # Perform variable selection for all causes
    while (length(candidate_vars) > 0) {
        aic_values <- numeric(length(candidate_vars))
        names(aic_values) <- candidate_vars

        for (i in seq_along(candidate_vars)) {
          var <- candidate_vars[i]
          current_vars <- c(selected_vars, var)
          formula_str <- paste("survival::Surv(", timevar, ", status_cs) ~", paste(current_vars, collapse = "+"))
          formula <- stats::as.formula(formula_str)

          model <- tryCatch(
            survival::survreg(formula, data = XYTrain_cause, dist = dist, x = FALSE, y = FALSE),
            error = function(e) {
              if (verbose) warning("Failed to fit model with var ", var, ": ", e$message)
              NULL
            }
          )

          if (!is.null(model)) {
            aic_values[i] <- stats::AIC(model)
          } else {
            aic_values[i] <- Inf # Penalize models that fail to fit
          }
        }

        best_candidate_idx <- which.min(aic_values)
        best_candidate_aic <- aic_values[best_candidate_idx]
        best_candidate_var <- candidate_vars[best_candidate_idx]

        # Add variable if it improves AIC
        if (best_candidate_aic < best_aic) {
          if (verbose) print(paste("Adding", best_candidate_var, "AIC:", round(best_candidate_aic, 2),
                                  "(Improvement:", round(best_aic - best_candidate_aic, 2), ")"))
          selected_vars <- c(selected_vars, best_candidate_var)
          candidate_vars <- setdiff(candidate_vars, best_candidate_var)
          best_aic <- best_candidate_aic
        } else {
          # If no improvement, but no variable has been selected yet, pick the best single variable
          if (length(selected_vars) == 0) {
            if (verbose) print(paste("No improvement over intercept-only, but keeping best single variable:", best_candidate_var, "AIC:", round(best_candidate_aic, 2)))
            selected_vars <- best_candidate_var
          } else {
            if (verbose) print(paste("No improvement adding remaining variables. Best AIC:", round(best_aic, 2)))
          }
          break # Stop if no variable improves AIC
        }
        # End while loop
      }
      # Robustness: If after selection, selected_vars is empty, always keep the best single variable
      if (length(selected_vars) == 0 && length(candidate_vars) > 0) {
        # Defensive: pick the best single variable by AIC
        aic_values <- numeric(length(candidate_vars))
        names(aic_values) <- candidate_vars
        for (i in seq_along(candidate_vars)) {
          var <- candidate_vars[i]
          formula_str <- paste("survival::Surv(", timevar, ", status_cs) ~", var)
          model <- tryCatch(
            survival::survreg(stats::as.formula(formula_str), data = XYTrain_cause, dist = dist, x = FALSE, y = FALSE),
            error = function(e) NULL
          )
          if (!is.null(model)) {
            aic_values[i] <- stats::AIC(model)
          } else {
            aic_values[i] <- Inf
          }
        }
        best_idx <- which.min(aic_values)
        if (length(best_idx) == 1 && is.finite(aic_values[best_idx])) {
          selected_vars <- candidate_vars[best_idx]
          if (verbose) print(paste("No variable improved over intercept-only, keeping best single variable:", selected_vars, "AIC:", round(aic_values[best_idx], 2)))
        } else {
          stop("No valid single-variable model could be fit after variable selection.")
        }
      }

    # ============================================================================
    # Fit the Final Selected Model
    # ============================================================================
    if (length(selected_vars) > 0) {
      final_formula_str <- paste("survival::Surv(", timevar, ", status_cs) ~", paste(selected_vars, collapse = "+"))
      final_model_cause <- tryCatch(
        survival::survreg(stats::as.formula(final_formula_str), data = XYTrain_cause,
                         dist = dist, x = TRUE, y = TRUE),
        error = function(e) {
          warning("Failed to fit SurvReg model for cause ", cause, ": ", e$message)
          NULL
        }
      )
    } else {
      # If no variables selected, use the intercept-only model
      if (verbose) warning("No variables selected by forward selection. Using intercept-only model.")
      final_model_cause <- tryCatch(
        survival::survreg(null_formula, data = XYTrain_cause, dist = dist, x = TRUE, y = TRUE),
        error = function(e) {
          warning("Failed to fit intercept-only model for cause ", cause, ": ", e$message)
          NULL
        }
      )
    }

    if (is.null(final_model_cause)) {
      warning("Failed to fit any model for cause ", cause, ". Skipping.")
      next
    }

    # ============================================================================
    # Create Baseline Model for Prediction
    # ============================================================================
    # Get linear predictors for training data
    train_linear_preds <- predict(final_model_cause, type = "linear")

    # IMPORTANT: Convert linear predictors to risk scores (negate for proper direction)
    # In survival models, higher linear predictors mean longer survival (lower risk)
    # But for CIF, we want higher risk scores to give higher CIF
    train_risk_scores <- -train_linear_preds

    # Create survival data for the cause-specific case
    train_survival_data <- data.frame(
      time = XYTrain_cause[[timevar]],
      event = XYTrain_cause$status_cs
    )

    # Use score2proba to create baseline hazard model
    baseline_info <- score2proba(
      datasurv = train_survival_data,
      score = train_risk_scores,
      conf.int = 0.95,
      which.est = "point"
    )

    # Store baseline model in the survreg model object
    final_model_cause$baseline_model <- baseline_info$model
    final_model_cause$baseline_sf <- baseline_info$sf
    final_model_cause$selected_vars <- selected_vars  # Store selected variables for this cause

    survreg_models_all_causes[[as.character(cause)]] <- final_model_cause
  }

  # The main model for the event of interest
  final_model <- survreg_models_all_causes[[as.character(primary_event_numeric)]]

  if (is.null(final_model)) {
  stop("Failed to fit SurvReg model for the event of interest (event_code = ", primary_event_code, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific parametric survival model fitting complete.\n")

  result <- list(
    survreg_model = final_model,
    survreg_models_all_causes = survreg_models_all_causes,  # All cause-specific models for Aalen-Johansen
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_survreg",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes,
    event_codes_numeric = event_codes_numeric,
    default_event_code = primary_event_code,
    default_event_code_numeric = primary_event_numeric,
    time_range = time_range,
    dist = dist
  )

  class(result) <- c("ml4t2e_cr_survreg", "list")
  return(result)
}

#'
#'
#'
#'
Predict_CRModel_SurvReg <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("'modelout' is missing")
  }
  if (!is.list(modelout) || !all(c("expvars", "default_event_code") %in% names(modelout))) {
    stop("'modelout' must be output from CRModel_SurvReg")
  }
  if (missing(newdata)) {
    stop("'newdata' is missing")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame")
  }

  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop("The following variables missing in newdata: ",
         paste(missing_vars, collapse = ", "))
  }

  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    event_of_interest <- modelout$default_event_code
  }

  event_of_interest <- as.character(event_of_interest)
  event_idx <- match(event_of_interest, modelout$event_codes)

  if (is.na(event_idx)) {
    stop("event_of_interest ", event_of_interest, " was not present in training data. Available event codes: ",
         paste(modelout$event_codes, collapse = ", "))
  }

  target_event_numeric <- modelout$event_codes_numeric[event_idx]

  use_native_times <- is.null(new_times)
  if (!use_native_times) {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    new_times <- sort(unique(new_times))
  }

  # ============================================================================
  # Prepare newdata
  # ============================================================================
  # Ensure factor levels match training data
  newdata_prepared <- newdata[, modelout$expvars, drop = FALSE]

  for (var in modelout$expvars) {
    if (var %in% names(modelout$varprof)) {
      varinfo <- modelout$varprof[[var]]
      # Handle factors
      if (is.table(varinfo)) {
        training_levels <- names(varinfo)
        if (is.factor(newdata_prepared[[var]])) {
          # Ensure factor has same levels as training
          new_levels <- levels(newdata_prepared[[var]])
          extra_levels <- setdiff(new_levels, training_levels)
          if (length(extra_levels) > 0) {
            warning("Factor '", var, "' has new levels in newdata: ",
                    paste(extra_levels, collapse = ", "),
                    ". These will be set to NA.")
          }
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
                                            levels = training_levels)
        } else if (is.character(newdata_prepared[[var]])) {
          # Convert character to factor with training levels
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
                                            levels = training_levels)
        }
      }
    }
  }

  # Check for NAs
  if (any(is.na(newdata_prepared))) {
    warning("Missing values in newdata will result in NA predictions")
  }

  # ============================================================================
  # Make Predictions using proper parametric approach
  # ============================================================================
  n_obs <- nrow(newdata_prepared)

  # Determine time grid for prediction
  surv_times <- if (use_native_times) modelout$times else new_times
  surv_times <- sort(unique(surv_times))

  if (length(surv_times) == 0) {
    stop("No valid time points available for prediction.")
  }

  # Compute cause-specific hazards for each observation
  cause_specific_hazards <- list()

  for (cause in modelout$event_codes_numeric) {
    cause_char <- as.character(cause)
    survreg_model <- modelout$survreg_models_all_causes[[cause_char]]

    if (is.null(survreg_model)) {
      warning("No model fitted for cause ", cause_char, ". This may indicate insufficient events for this cause.")
      next
    }

    # message("Predicting for cause ", cause_char, " with ", length(survreg_model$selected_vars), " variables")

    # Use the variables that this specific model was trained on
    model_vars <- survreg_model$selected_vars
    if (is.null(model_vars) || length(model_vars) == 0) {
      model_vars <- modelout$expvars  # fallback to all variables
    }
    
    # Subset newdata to only the variables this model was trained on
    newdata_cause <- newdata_prepared[, model_vars, drop = FALSE]

    # Linear predictors for newdata
    linear_preds <- tryCatch(
      stats::predict(survreg_model, newdata = newdata_cause, type = "lp"),
      error = function(e) {
        stop("Failed to obtain linear predictors for cause ", cause_char, ": ", e$message)
      }
    )

    # Ensure vector form
    linear_preds <- as.numeric(linear_preds)

    # Compute hazard rates using the parametric form
    hazard_matrix <- matrix(NA_real_, nrow = length(surv_times), ncol = n_obs)
    for (j in seq_len(n_obs)) {
      # For exponential distribution, hazard = 1 / exp(lp)
      # For Weibull, hazard = (shape / exp(lp)) * (t / exp(lp))^(shape-1)
      if (survreg_model$dist == "exponential") {
        hazard_matrix[, j] <- rep(1 / exp(linear_preds[j]), length(surv_times))
      } else {
        # Weibull hazard
        shape <- 1 / survreg_model$scale
        scale_param <- exp(linear_preds[j])
        hazard_matrix[, j] <- shape * (surv_times / scale_param)^(shape - 1) / scale_param
      }
    }

    cause_specific_hazards[[cause_char]] <- hazard_matrix
  }

  if (length(cause_specific_hazards) == 0) {
    stop("No cause-specific models available for prediction")
  }

  # Compute CIF using exact parametric formulas for competing risks
  n_times <- length(surv_times)
  cif_matrix <- matrix(0, nrow = n_times, ncol = n_obs)

  for (j in seq_len(n_obs)) {
    # Extract hazard rates at time 0 (constant for exponential, but for Weibull we need to handle differently)
    # For exponential, hazard is constant over time
    hazard_rates <- sapply(cause_specific_hazards, function(hmat) hmat[1, j])
    
    if (anyNA(hazard_rates)) {
      cif_matrix[, j] <- NA
      next
    }
    
    total_hazard <- sum(hazard_rates)
    
    if (total_hazard == 0) {
      # No risk, CIF remains 0
      cif_matrix[, j] <- 0
      next
    }
    
    target_cause_char <- as.character(target_event_numeric)
    if (target_cause_char %in% names(hazard_rates)) {
      target_hazard <- hazard_rates[target_cause_char]
      
      # For exponential distribution, use exact CIF formula
      # CIF_j(t) = (λ_j / Σλ_k) * (1 - exp(-Σλ_k * t))
      proportion <- target_hazard / total_hazard
      cif_matrix[, j] <- proportion * (1 - exp(-total_hazard * surv_times))
    }
  }

  # Ensure bounds and monotonicity
  cif_matrix <- pmin(pmax(cif_matrix, 0), 1)
  for (j in seq_len(n_obs)) {
    cif_matrix[, j] <- cummax(cif_matrix[, j])
  }

  if (!use_native_times) {
    interpolated <- cifMatInterpolaltor(
      probsMat = cif_matrix,
      times = surv_times,
      new_times = new_times
    )
    result_cifs <- interpolated
    result_times <- new_times
  } else {
    result_cifs <- cif_matrix
    result_times <- surv_times
  }

  list(
    CIFs = result_cifs,
    Times = result_times
  )
}
