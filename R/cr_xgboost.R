#' @title CRModel_xgboost
#'
#' @description Fit an XGBoost model for competing risks outcomes using cause-specific modeling.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param eta learning rate (default: 0.01)
#' @param max_depth maximum depth of trees (default: 5)
#' @param nrounds number of boosting rounds (default: 100)
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#' @param ... additional parameters passed to xgboost
#'
#' @return a list with the following components:
#'   \item{xgb_model}{the fitted cause-specific XGBoost model object}
#'   \item{times}{vector of unique event times in the training data for the event of interest}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_xgboost"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{feature_names}{character vector of feature names used in XGBoost}
#'
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom stats model.matrix
#' @export
CRModel_xgboost <- function(data, expvars, timevar, eventvar, failcode = 1,
                           eta = 0.01, max_depth = 5, nrounds = 100,
                           ntimes = 50, verbose = FALSE, ...) {

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
  # Model Fitting - Cause-Specific XGBoost Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  all_event_types <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))

  if (verbose) {
    cat("Fitting cause-specific XGBoost models for all event types:", paste(all_event_types, collapse = ", "), "\n")
  }

  # Store models for all event types
  xgb_models_all_causes <- vector("list", length(all_event_types))
  names(xgb_models_all_causes) <- as.character(all_event_types)

  # Fit a separate XGBoost model for each event type
  for (cause in all_event_types) {
    if (verbose) cat("Fitting XGBoost model for event type", cause, "...\n")

    # Create cause-specific data: event = 1 if this cause, 0 otherwise (censored or competing)
    XYTrain_cause <- XYTrain
    XYTrain_cause$status_cs <- ifelse(XYTrain[[eventvar]] == cause, 1, 0)

    if (sum(XYTrain_cause$status_cs) < 5) {
      warning("Fewer than 5 events of type ", cause, ". Skipping this cause.")
      xgb_models_all_causes[[as.character(cause)]] <- NULL
      next
    }

    # ============================================================================
    # Prepare Data Matrix for XGBoost
    # ============================================================================
    # Prepare data matrix (handle factors using model.matrix)
    X <- stats::model.matrix(~-1+., XYTrain_cause[, expvars, drop = FALSE])
    if (cause == all_event_types[1]) {
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

    # Store baseline model in the XGBoost model object
    xgb_model_cause$baseline_model <- baseline_info$model
    xgb_model_cause$baseline_sf <- baseline_info$sf

    xgb_models_all_causes[[as.character(cause)]] <- xgb_model_cause
  }

  # The main model for the event of interest
  xgb_model <- xgb_models_all_causes[[as.character(failcode)]]

  if (is.null(xgb_model)) {
    stop("Failed to fit XGBoost model for the event of interest (failcode = ", failcode, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific XGBoost model fitting complete.\n")

  result <- list(
    xgb_model = xgb_model,
    xgb_models_all_causes = xgb_models_all_causes,  # All cause-specific models for Aalen-Johansen
    all_event_types = all_event_types,  # Event type codes
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_xgboost",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
    time_range = time_range,
    feature_names = feature_names
  )

  class(result) <- c("ml4t2e_cr_xgboost", "list")
  return(result)
}

#' @title Predict_CRModel_xgboost
#'
#' @description Get predictions from a fitted cause-specific XGBoost competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_xgboost'
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
#' @importFrom stats model.matrix
#' @importFrom xgboost xgb.DMatrix
#' @export
Predict_CRModel_xgboost <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_xgboost")) {
    stop("'modelout' must be output from CRModel_xgboost")
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
  # Get survival predictions from ALL cause-specific XGBoost models
  cause_specific_survs <- vector("list", length(modelout$all_event_types))
  names(cause_specific_survs) <- as.character(modelout$all_event_types)

  # Predict survival for each cause
  for (cause in modelout$all_event_types) {
    cause_char <- as.character(cause)

    if (is.null(modelout$xgb_models_all_causes[[cause_char]])) {
      # If model wasn't fitted for this cause, assume no events (S(t) = 1)
      warning("No model available for cause ", cause, ". Assuming S(t) = 1 for all times.")
      # We'll handle this in aalenJohansenCIF by providing a matrix of 1s
      next
    }

    # Get linear predictors from the XGBoost model for this cause
    linear_preds_cause <- predict(modelout$xgb_models_all_causes[[cause_char]], X_new)

    # Use the stored baseline hazard from the training fit
    sf_cause <- tryCatch(
      survival::survfit(
        modelout$xgb_models_all_causes[[cause_char]]$baseline_model,
        newdata = data.frame("score" = linear_preds_cause),
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