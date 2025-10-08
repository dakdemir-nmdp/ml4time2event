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

  # ============================================================================
  # Cause-Specific Data Preparation
  # ============================================================================
  # Create cause-specific outcome: 1 for event of interest, 0 for censored/other events
  XYTrain_cs <- XYTrain
  XYTrain_cs$status_cs <- ifelse(XYTrain[[eventvar]] == failcode, 1,
                                ifelse(XYTrain[[eventvar]] == 0, 0, 0))

  # Get time range for the event of interest
  event_times <- XYTrain_cs[[timevar]][XYTrain_cs$status_cs == 1]
  time_range <- c(0, max(event_times))

  # ============================================================================
  # Prepare Data Matrix for XGBoost
  # ============================================================================
  # Prepare data matrix (handle factors using model.matrix)
  X <- stats::model.matrix(~-1+., XYTrain_cs[, expvars, drop = FALSE])
  feature_names <- colnames(X)

  # Prepare labels for XGBoost: negative time for censored, positive time for event
  ytimes <- XYTrain_cs[[timevar]]
  yevents <- XYTrain_cs$status_cs
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
  if (verbose) cat("Training cause-specific XGBoost model...\n")

  # Create DMatrix
  dtrain <- xgboost::xgb.DMatrix(data = X, label = yTrain)

  # Train the model
  xgb_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = ifelse(verbose, 1, 0)
  )

  # ============================================================================
  # Create Baseline Model for Prediction
  # ============================================================================
  # Get linear predictors for training data
  train_linear_preds <- predict(xgb_model, X)

  # Create survival data for the cause-specific case
  train_survival_data <- data.frame(
    time = XYTrain_cs[[timevar]],
    event = XYTrain_cs$status_cs
  )

  # Use score2proba to create baseline hazard model
  baseline_info <- score2proba(
    datasurv = train_survival_data,
    score = train_linear_preds,
    conf.int = 0.95,
    which.est = "point"
  )

  # Store baseline model in the XGBoost model object
  xgb_model$baseline_model <- baseline_info$model
  xgb_model$baseline_sf <- baseline_info$sf

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific XGBoost model fitting complete.\n")

  result <- list(
    xgb_model = xgb_model,
    times = sort(unique(XYTrain_cs[[timevar]][XYTrain_cs$status_cs == 1])),
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
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats model.matrix
#' @importFrom xgboost xgb.DMatrix
#' @export
Predict_CRModel_xgboost <- function(modelout, newdata, newtimes = NULL) {

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
  # Make Predictions
  # ============================================================================
  # Get linear predictors from the XGBoost model
  linear_preds <- predict(modelout$xgb_model, X_new)

  # Use the stored baseline hazard from the training fit
  # and the new scores (linear_preds) to get survival probabilities for the new data
  sf <- survival::survfit(
    modelout$xgb_model$baseline_model, # Use the stored Cox model
    newdata = data.frame("score" = linear_preds),
    conf.int = 0.95
  )

  # Extract survival probabilities
  surv_probs <- sf$surv # Matrix: rows=times, cols=observations
  surv_times <- sf$time

  # Ensure matrix format
  if (!is.matrix(surv_probs)) {
    surv_probs <- matrix(surv_probs, ncol = 1)
  }

  # For cause-specific model, CIF is approximately 1 - S(t)
  # where S(t) is the cause-specific survival function
  cif_matrix <- 1 - surv_probs

  # Ensure CIFs are properly bounded [0,1]
  cif_matrix <- pmin(pmax(cif_matrix, 0), 1)

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(newtimes)) {
    # Return predictions in native time grid: [observations, times]
    result_cifs <- t(cif_matrix)  # cif_matrix is [times, observations], transpose to [observations, times]
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

    # cifMatInterpolaltor returns [newtimes, observations], transpose to [observations, newtimes]
    result_cifs <- t(pred_cifs)
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