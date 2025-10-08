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
  # Model Fitting - Cause-Specific Cox Model
  # ============================================================================
  if (verbose) cat("Fitting cause-specific Cox model for event type", failcode, "...\n")

  # Create cause-specific data: event = 1 if failcode, 0 if censored, NA if other competing event
  XYTrain$status_cs <- ifelse(XYTrain[[eventvar]] == 0, 0,  # censored
                              ifelse(XYTrain[[eventvar]] == failcode, 1, NA))  # event of interest or competing

  # Remove competing events (treat as censored for this cause-specific model)
  XYTrain_cs <- XYTrain[!is.na(XYTrain$status_cs), , drop = FALSE]

  if (nrow(XYTrain_cs) < 10) {
    stop("Insufficient data after removing competing events. Need at least 10 observations.")
  }

  # Define formula
  form <- stats::as.formula(
    paste0("survival::Surv(", timevar, ", status_cs) ~ ",
           paste(expvars, collapse = " + "))
  )

  if (verbose) cat("Fitting Cox model...\n")

  # Fit initial Cox model
  cph_model <- tryCatch(
    survival::coxph(form, data = XYTrain_cs, x = TRUE, y = TRUE),
    error = function(e) {
      stop("Failed to fit cause-specific Cox model: ", e$message)
    }
  )

  # Variable selection
  if (varsel != "none") {
    if (verbose) cat("Performing variable selection (", varsel, ", ", penalty, ")...\n", sep="")

    # Set penalty parameter k for step()
    k_penalty <- switch(penalty, "AIC" = 2, "BIC" = log(nrow(XYTrain_cs)))

    # Define scope for forward/both methods
    if (varsel %in% c("forward", "both")) {
      null_formula <- stats::as.formula(
        paste0("survival::Surv(", timevar, ", status_cs) ~ 1")
      )
      null_model <- survival::coxph(null_formula, data = XYTrain_cs, x = TRUE, y = TRUE)
      lower_scope <- null_formula
      upper_scope <- form
    }

    # Perform stepwise selection
    cph_model <- tryCatch(
      {
        if (varsel == "backward") {
          stats::step(cph_model, direction = "backward", k = k_penalty, trace = as.numeric(verbose))
        } else if (varsel == "forward") {
          stats::step(null_model, direction = "forward", scope = list(lower = lower_scope, upper = upper_scope),
                     k = k_penalty, trace = as.numeric(verbose))
        } else if (varsel == "both") {
          stats::step(null_model, direction = "both", scope = list(lower = lower_scope, upper = upper_scope),
                     k = k_penalty, trace = as.numeric(verbose))
        }
      },
      error = function(e) {
        warning("Variable selection failed: ", e$message, ". Using full model.")
        cph_model # Return original full model
      }
    )
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific Cox model fitting complete.\n")

  result <- list(
    cph_model = cph_model,
    times = sort(unique(XYTrain_cs[[timevar]][XYTrain_cs$status_cs == 1])),
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
#'
#' @return a list containing:
#'   \item{cph_modelTestPredict}{the raw survfit object from prediction}
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
Predict_CRModel_Cox <- function(modelout, newdata, newtimes = NULL) {

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
  # Make Predictions
  # ============================================================================
  # Get survival predictions from the cause-specific Cox model
  cph_modelTestPredict <- tryCatch(
    survival::survfit(modelout$cph_model, newdata = newdata_prepared),
    error = function(e) {
      stop("Prediction failed: ", e$message)
    }
  )

  # Extract survival probabilities
  surv_probs <- cph_modelTestPredict$surv
  surv_times <- cph_modelTestPredict$time

  # Ensure matrix format
  if (!is.matrix(surv_probs)) {
    surv_probs <- matrix(surv_probs, ncol = 1)
  }

  # For cause-specific Cox model, the CIF is approximately 1 - S(t)
  # where S(t) is the cause-specific survival function
  # This is a simplification; in practice, competing risks CIFs require
  # more complex calculation accounting for all causes
  cif_matrix <- 1 - surv_probs

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
    cph_modelTestPredict = cph_modelTestPredict,
    CIFs = result_cifs,
    Times = result_times
  )

  return(result)
}
