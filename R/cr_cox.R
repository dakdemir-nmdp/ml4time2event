#' @title CRModel_Cox
#'
#' @description Fit a cause-specific Cox model for competing risks outcomes with optional variable selection.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param varsel character string specifying variable selection method:
#'   "none" (default, no selection),
#'   "backward" (backward elimination),
#'   "forward" (forward selection),
#'   "both" (stepwise selection)
#' @param penalty character string specifying penalty criterion for stepwise methods: "AIC" (default) or "BIC"
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{cph_model}{the fitted cause-specific Cox model object}
#'   \item{times}{vector of unique event times in the training data}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_cox"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{varsel_method}{character string indicating variable selection method used}
#'   \item{time_range}{vector with min and max observed event times}
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit cause-specific Cox model without variable selection
#' model1 <- CRModel_Cox(data, expvars = c("x1", "x2"),
#'                      timevar = "time", eventvar = "event")
#'
#' # Fit with backward selection (AIC)
#' model2 <- CRModel_Cox(data, expvars = c("x1", "x2", "x3"),
#'                      timevar = "time", eventvar = "event",
#'                      varsel = "backward", penalty = "AIC")
#'
#' # Fit with forward selection (BIC)
#' model3 <- CRModel_Cox(data, expvars = c("x1", "x2", "x3"),
#'                      timevar = "time", eventvar = "event",
#'                      varsel = "forward", penalty = "BIC")
#' }
CRModel_Cox <- function(data, expvars, timevar, eventvar, failcode = 1,
                       varsel = "none", penalty = "AIC",
                       ntimes = 50, verbose = FALSE) {

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
  varsel <- match.arg(varsel, c("none", "backward", "forward", "both"))
  penalty <- match.arg(penalty, c("AIC", "BIC"))

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
  # Model Fitting - Cause-Specific Cox Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  all_event_types <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))

  if (verbose) {
    cat("Fitting cause-specific Cox models for all event types:", paste(all_event_types, collapse = ", "), "\n")
  }

  # Store models for all event types
  cph_models_all_causes <- vector("list", length(all_event_types))
  names(cph_models_all_causes) <- as.character(all_event_types)

  # Fit a separate Cox model for each event type
  for (cause in all_event_types) {
    if (verbose) cat("Fitting Cox model for event type", cause, "...\n")

    # Create cause-specific data: event = 1 if this cause, 0 otherwise (censored or competing)
    XYTrain_cause <- XYTrain
    XYTrain_cause$status_cs <- ifelse(XYTrain[[eventvar]] == cause, 1, 0)

    if (sum(XYTrain_cause$status_cs) < 5) {
      warning("Fewer than 5 events of type ", cause, ". Skipping this cause.")
      cph_models_all_causes[[as.character(cause)]] <- NULL
      next
    }

    # Define formula
    form <- stats::as.formula(
      paste0("survival::Surv(", timevar, ", status_cs) ~ ",
             paste(expvars, collapse = " + "))
    )

    # Fit initial Cox model
    cph_model_cause <- tryCatch(
      survival::coxph(form, data = XYTrain_cause, x = TRUE, y = TRUE),
      error = function(e) {
        warning("Failed to fit Cox model for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(cph_model_cause)) {
      next
    }

    # Variable selection (only for the event of interest to save time)
    if (varsel != "none" && cause == failcode) {
      if (verbose) cat("Performing variable selection (", varsel, ", ", penalty, ")...\n", sep="")

      # Set penalty parameter k for step()
      k_penalty <- switch(penalty, "AIC" = 2, "BIC" = log(nrow(XYTrain_cause)))

      # Define scope for forward/both methods
      if (varsel %in% c("forward", "both")) {
        null_formula <- stats::as.formula(
          paste0("survival::Surv(", timevar, ", status_cs) ~ 1")
        )
        null_model <- survival::coxph(null_formula, data = XYTrain_cause, x = TRUE, y = TRUE)
        lower_scope <- null_formula
        upper_scope <- form
      }

      # Perform stepwise selection
      cph_model_cause <- tryCatch(
        {
          if (varsel == "backward") {
            stats::step(cph_model_cause, direction = "backward", k = k_penalty, trace = as.numeric(verbose))
          } else if (varsel == "forward") {
            stats::step(null_model, direction = "forward", scope = list(lower = lower_scope, upper = upper_scope),
                       k = k_penalty, trace = as.numeric(verbose))
          } else if (varsel == "both") {
            stats::step(null_model, direction = "both", scope = list(lower = lower_scope, upper = upper_scope),
                       k = k_penalty, trace = as.numeric(verbose))
          }
        },
        error = function(e) {
          warning("Variable selection failed for cause ", cause, ": ", e$message, ". Using full model.")
          cph_model_cause # Return original full model
        }
      )
    }

    cph_models_all_causes[[as.character(cause)]] <- cph_model_cause
  }

  # The main model for the event of interest
  cph_model <- cph_models_all_causes[[as.character(failcode)]]

  if (is.null(cph_model)) {
    stop("Failed to fit Cox model for the event of interest (failcode = ", failcode, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific Cox model fitting complete.\n")

  result <- list(
    cph_model = cph_model,
    cph_models_all_causes = cph_models_all_causes,  # All cause-specific models for Aalen-Johansen
    all_event_types = all_event_types,  # Event type codes
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_cox",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
    varsel_method = varsel,
    time_range = time_range
  )

  class(result) <- c("ml4t2e_cr_cox", "list")
  return(result)
}

#' @title Predict_CRModel_Cox
#'
#' @description Get predictions from a fitted cause-specific Cox competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_Cox'
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
#' @importFrom survival survfit
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit model
#' model <- CRModel_Cox(data, expvars = c("x1", "x2"),
#'                     timevar = "time", eventvar = "event")
#'
#' # Predict on test data
#' preds <- Predict_CRModel_Cox(model, test_data)
#'
#' # Predict at specific times
#' preds_custom <- Predict_CRModel_Cox(model, test_data,
#'                                    newtimes = c(30, 60, 90, 180, 365))
#' }
Predict_CRModel_Cox <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_cox")) {
    stop("'modelout' must be output from CRModel_Cox")
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
  # Get survival predictions from ALL cause-specific Cox models
  cause_specific_survs <- vector("list", length(modelout$all_event_types))
  names(cause_specific_survs) <- as.character(modelout$all_event_types)

  # Predict survival for each cause
  for (cause in modelout$all_event_types) {
    cause_char <- as.character(cause)

    if (is.null(modelout$cph_models_all_causes[[cause_char]])) {
      # If model wasn't fitted for this cause, assume no events (S(t) = 1)
      warning("No model available for cause ", cause, ". Assuming S(t) = 1 for all times.")
      # We'll handle this in aalenJohansenCIF by providing a matrix of 1s
      next
    }

    cph_predict <- tryCatch(
      survival::survfit(modelout$cph_models_all_causes[[cause_char]], newdata = newdata_prepared),
      error = function(e) {
        warning("Prediction failed for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(cph_predict)) {
      next
    }

    # Extract survival probabilities
    surv_probs_cause <- cph_predict$surv
    surv_times_cause <- cph_predict$time

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
  surv_times <- cph_predict$time

  # Use Aalen-Johansen estimator to calculate proper CIF
  cif_matrix <- aalenJohansenCIF(
    cause_specific_survs = cause_specific_survs,
    times = surv_times,
    event_of_interest = failcode
  )

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(newtimes)) {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times, observations]
    result_times <- surv_times
  } else {
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
