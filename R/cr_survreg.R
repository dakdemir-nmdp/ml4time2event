#' @title CRModel_SurvReg
#'
#' @description Fit a parametric survival regression model for competing risks outcomes using cause-specific modeling.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param event_codes character or numeric vector identifying the event code(s)
#'   to model. If NULL (default), all non-zero event codes observed in the data
#'   are used. The first entry defines the default event of interest.
#' @param dist distribution for the parametric model, one of "weibull", "exponential", "gaussian", "logistic","lognormal" or "loglogistic".
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{survreg_model}{the fitted cause-specific parametric survival model object}
#'   \item{times}{vector of unique event times in the training data for the event of interest}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_survreg"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{event_codes}{character vector of event codes included in the model}
#'   \item{event_codes_numeric}{numeric vector of event codes included}
#'   \item{default_event_code}{character scalar for the default event code}
#'   \item{default_event_code_numeric}{numeric scalar for the default event code}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{dist}{the distribution used for the parametric model}
#'
#' @importFrom survival survreg Surv
#' @importFrom stats AIC as.formula predict quantile complete.cases
#' @export

#' @param event_of_interest optional character or numeric scalar indicating a specific event code
#'   that should be prioritized as the primary event of interest. If provided, this
#'   event code must be one of the codes specified in 'event_codes'.
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
  if (verbose && cause == primary_event_numeric) cat("Performing forward selection with AIC for event", cause, "...\n")

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
  if (verbose && cause == primary_event_numeric) print(paste("Initial AIC (Intercept only):", round(best_aic, 2)))

    # Only do variable selection for the main event to save time
  if (cause == primary_event_numeric) {
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
    } else {
      # For other causes, use all variables to save time
      selected_vars <- expvars
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
  if (verbose && cause == primary_event_numeric) warning("No variables selected by forward selection. Using intercept-only model.")
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

#' @title Predict_CRModel_SurvReg
#'
#' @description Get predictions from a fitted cause-specific parametric survival competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_SurvReg'
#' @param newdata data frame with new observations for prediction
#' @param new_times optional numeric vector of time points for prediction.
#'   If NULL (default), uses the times from the training data.
#'   Can be any positive values - interpolation handles all time points.
#' @param event_of_interest character or numeric scalar indicating the event code
#'   to predict. If NULL (default), uses the event code stored during training.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @importFrom survival psurvreg
#' @export
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

  # Compute cause-specific survival curves for each observation
  cause_specific_survival <- list()

  for (cause in modelout$event_codes_numeric) {
    cause_char <- as.character(cause)
    survreg_model <- modelout$survreg_models_all_causes[[cause_char]]

    if (is.null(survreg_model)) {
      next
    }

    # Linear predictors for newdata
    linear_preds <- tryCatch(
      stats::predict(survreg_model, newdata = newdata_prepared, type = "lp"),
      error = function(e) {
        stop("Failed to obtain linear predictors for cause ", cause_char, ": ", e$message)
      }
    )

    # Ensure vector form
    linear_preds <- as.numeric(linear_preds)

    # Compute survival curves using the parametric form
    surv_matrix <- matrix(NA_real_, nrow = length(surv_times), ncol = n_obs)
    for (j in seq_len(n_obs)) {
      surv_vals <- 1 - survival::psurvreg(
        q = surv_times,
        mean = linear_preds[j],
        scale = survreg_model$scale,
        distribution = survreg_model$dist
      )
      surv_matrix[, j] <- pmax(pmin(surv_vals, 1), 0)
    }

    cause_specific_survival[[cause_char]] <- surv_matrix
  }

  if (length(cause_specific_survival) == 0) {
    result_times <- surv_times
    if (result_times[1] > 0) {
      result_times <- c(0, result_times)
    }
    result_cifs <- matrix(0, nrow = length(result_times), ncol = n_obs)
    return(list(CIFs = result_cifs, Times = result_times))
  }

  # Ensure time zero is included for stability
  if (surv_times[1] > 0) {
    surv_times <- c(0, surv_times)
    for (cause_char in names(cause_specific_survival)) {
      surv_matrix <- cause_specific_survival[[cause_char]]
      zero_row <- matrix(1, nrow = 1, ncol = n_obs)
      cause_specific_survival[[cause_char]] <- rbind(zero_row, surv_matrix)
    }
  } else {
    for (cause_char in names(cause_specific_survival)) {
      surv_matrix <- cause_specific_survival[[cause_char]]
      surv_matrix[1, ] <- 1
      cause_specific_survival[[cause_char]] <- surv_matrix
    }
  }

  n_times <- length(surv_times)

  # Convert survival curves to cumulative hazards
  cause_cumhaz <- lapply(
    cause_specific_survival,
    function(surv_matrix) -log(pmax(surv_matrix, 1e-12))
  )

  # Overall survival from sum of cumulative hazards
  overall_cumhaz <- Reduce(`+`, cause_cumhaz)
  overall_survival <- exp(-overall_cumhaz)

  # Compute CIF for target cause
  target_cause_char <- as.character(target_event_numeric)
  target_cumhaz <- cause_cumhaz[[target_cause_char]]

  if (is.null(target_cumhaz)) {
    stop("No model available for event_of_interest = ", target_cause_char)
  }

  cif_matrix <- matrix(0, nrow = n_times, ncol = n_obs)

  for (t in seq_len(n_times)) {
    prev_surv <- if (t == 1) rep(1, n_obs) else overall_survival[t - 1, ]
    hazard_increment <- target_cumhaz[t, ] - if (t == 1) 0 else target_cumhaz[t - 1, ]
    hazard_increment <- pmax(hazard_increment, 0)

    if (t == 1) {
      cif_matrix[t, ] <- hazard_increment
    } else {
      cif_matrix[t, ] <- cif_matrix[t - 1, ] + prev_surv * hazard_increment
    }
  }

  # Enforce bounds and monotonicity
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
