#' @title CRModel_SurvReg
#'
#' @description Fit a parametric survival regression model for competing risks outcomes using cause-specific modeling.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param failcode integer, the code for the event of interest (default: 1)
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
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{dist}{the distribution used for the parametric model}
#'
#' @importFrom survival survreg Surv
#' @importFrom stats AIC as.formula predict quantile
#' @export
CRModel_SurvReg <- function(data, expvars, timevar, eventvar, failcode = 1,
                           dist = "exponential", ntimes = 50, verbose = FALSE) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("'expvars' must be a non-empty character vector")
  }
  if (!timevar %in% colnames(data)) {
    stop("'timevar' not found in data: ", timevar)
  }
  if (!eventvar %in% colnames(data)) {
    stop("'eventvar' not found in data: ", eventvar)
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following expvars not found in data: ", paste(missing_vars, collapse=", "))
  }
  if (!is.numeric(failcode) || length(failcode) != 1 || failcode < 1) {
    stop("'failcode' must be a positive integer")
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

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == failcode]
  if (length(event_times) == 0) {
    stop("No events of type ", failcode, " in training data. Cannot fit competing risks model.")
  }

  # Store event time range for reference
  time_range <- c(0, max(event_times))

  # ============================================================================
  # Model Fitting - Cause-Specific SurvReg Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  all_event_types <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))

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
    if (verbose && cause == failcode) cat("Performing forward selection with AIC for event", cause, "...\n")

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
    if (verbose && cause == failcode) print(paste("Initial AIC (Intercept only):", round(best_aic, 2)))

    # Only do variable selection for the main event to save time
    if (cause == failcode) {
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
          if (verbose) print(paste("No improvement adding remaining variables. Best AIC:", round(best_aic, 2)))
          break # Stop if no variable improves AIC
        }
      } # End while loop
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
      if (verbose && cause == failcode) warning("No variables selected by forward selection. Using intercept-only model.")
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
  final_model <- survreg_models_all_causes[[as.character(failcode)]]

  if (is.null(final_model)) {
    stop("Failed to fit SurvReg model for the event of interest (failcode = ", failcode, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific parametric survival model fitting complete.\n")

  result <- list(
    survreg_model = final_model,
    survreg_models_all_causes = survreg_models_all_causes,  # All cause-specific models for Aalen-Johansen
    all_event_types = all_event_types,  # Event type codes
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_survreg",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
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
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the times from the training data.
#'   Can be any positive values - interpolation handles all time points.
#' @param failcode integer, the code for the event of interest for CIF prediction.
#'   If NULL (default), uses the failcode from the model training.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @export
Predict_CRModel_SurvReg <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_survreg")) {
    stop("'modelout' must be output from CRModel_SurvReg")
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

  # Handle failcode parameter
  if (is.null(failcode)) {
    failcode <- modelout$failcode  # Use the failcode from training
  } else {
    if (!is.numeric(failcode) || length(failcode) != 1 || failcode < 1) {
      stop("'failcode' must be a positive integer")
    }
    if (!failcode %in% modelout$all_event_types) {
      stop("failcode ", failcode, " was not present in training data. Available event types: ",
           paste(modelout$all_event_types, collapse = ", "))
    }
  }

  # Generate default times if not specified
  if (is.null(newtimes)) {
    # Create a reasonable default grid: include 0 and event times
    newtimes <- sort(unique(c(0, modelout$times)))
  } else {
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))
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
  # Make Predictions using Aalen-Johansen Estimator
  # ============================================================================
  # Get survival predictions from ALL cause-specific SurvReg models
  cause_specific_survs <- vector("list", length(modelout$all_event_types))
  names(cause_specific_survs) <- as.character(modelout$all_event_types)

  # Predict survival for each cause
  for (cause in modelout$all_event_types) {
    cause_char <- as.character(cause)

    if (is.null(modelout$survreg_models_all_causes[[cause_char]])) {
      # If model wasn't fitted for this cause, assume no events (S(t) = 1)
      warning("No model available for cause ", cause, ". Assuming S(t) = 1 for all times.")
      # We'll handle this in aalenJohansenCIF by providing a matrix of 1s
      next
    }

    # Get linear predictors from the SurvReg model for this cause
    linear_preds_cause <- stats::predict(modelout$survreg_models_all_causes[[cause_char]],
                                         newdata = newdata_prepared,
                                         type = "linear")

    # IMPORTANT: In survival models, higher linear predictors usually mean LONGER survival (lower risk)
    # But for competing risks CIF, we want higher risk scores to give higher CIF
    # So we need to negate the linear predictors to get proper risk scores
    risk_scores_cause <- -linear_preds_cause

    # Use the stored baseline hazard from the training fit
    sf_cause <- tryCatch(
      survival::survfit(
        modelout$survreg_models_all_causes[[cause_char]]$baseline_model,
        newdata = data.frame("score" = risk_scores_cause),
        conf.int = 0.95
      ),
      error = function(e) {
        warning("Prediction failed for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(sf_cause)) {
      next
    }

    # Extract survival probabilities
    surv_probs_cause <- sf_cause$surv
    surv_times_cause <- sf_cause$time

    # Ensure matrix format
    if (!is.matrix(surv_probs_cause)) {
      surv_probs_cause <- matrix(surv_probs_cause, ncol = 1)
    }

    # Store in the list (rows = times, cols = observations)
    cause_specific_survs[[cause_char]] <- surv_probs_cause
  }

  # Remove NULL entries
  cause_specific_survs <- cause_specific_survs[!sapply(cause_specific_survs, is.null)]

  if (length(cause_specific_survs) == 0) {
    stop("No valid cause-specific survival predictions could be generated")
  }

  # Get the time grid from the first model (they should all have similar times)
  surv_times <- sf_cause$time

  # Use Aalen-Johansen estimator to calculate proper CIF
  cif_matrix <- aalenJohansenCIF(
    cause_specific_survs = cause_specific_survs,
    times = surv_times,
    event_of_interest = failcode
  )

  # ============================================================================
  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (!is.null(newtimes)) {
    # Interpolate to new time points
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))

    # Use the standard CIF interpolation utility function
    pred_cifs <- cifMatInterpolaltor(
      probsMat = t(cif_matrix),  # cifMatInterpolaltor expects [observations, times]
      times = surv_times,
      newtimes = newtimes
    )

    # cifMatInterpolaltor returns [newtimes, observations], keep as [times, observations]
    result_cifs <- pred_cifs
    result_times <- newtimes
  } else {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times, observations]
    result_times <- surv_times
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