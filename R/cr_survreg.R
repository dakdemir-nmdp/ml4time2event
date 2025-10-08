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
  # Forward Selection with AIC (adapted for cause-specific modeling)
  # ============================================================================
  if (verbose) cat("Performing forward selection with AIC...\n")

  selected_vars <- c()
  candidate_vars <- expvars

  # Start with intercept-only model AIC
  null_formula <- stats::as.formula(paste("survival::Surv(", timevar, ", status_cs) ~ 1"))
  null_model <- tryCatch(
    survival::survreg(null_formula, data = XYTrain_cs, dist = dist, x = FALSE, y = FALSE),
    error = function(e) {
      if (verbose) warning("Failed to fit null model: ", e$message)
      NULL
    }
  )
  if (is.null(null_model)) stop("Failed to fit intercept-only model.")
  best_aic <- stats::AIC(null_model)
  if (verbose) print(paste("Initial AIC (Intercept only):", round(best_aic, 2)))

  while (length(candidate_vars) > 0) {
    aic_values <- numeric(length(candidate_vars))
    names(aic_values) <- candidate_vars

    for (i in seq_along(candidate_vars)) {
      var <- candidate_vars[i]
      current_vars <- c(selected_vars, var)
      formula_str <- paste("survival::Surv(", timevar, ", status_cs) ~", paste(current_vars, collapse = "+"))
      formula <- stats::as.formula(formula_str)

      model <- tryCatch(
        survival::survreg(formula, data = XYTrain_cs, dist = dist, x = FALSE, y = FALSE),
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

  # ============================================================================
  # Fit the Final Selected Model
  # ============================================================================
  if (length(selected_vars) > 0) {
    final_formula_str <- paste("survival::Surv(", timevar, ", status_cs) ~", paste(selected_vars, collapse = "+"))
    final_model <- survival::survreg(stats::as.formula(final_formula_str), data = XYTrain_cs,
                                     dist = dist, x = TRUE, y = TRUE)
  } else {
    # If no variables selected, use the intercept-only model
    if (verbose) warning("No variables selected by forward selection. Using intercept-only model.")
    final_model <- survival::survreg(null_formula, data = XYTrain_cs, dist = dist, x = TRUE, y = TRUE)
    selected_vars <- character(0)
  }

  if (is.null(final_model)) {
    stop("No valid model could be fit even after selection. Please check the data and variables.")
  }

  # ============================================================================
  # Create Baseline Model for Prediction
  # ============================================================================
  # Get linear predictors for training data
  train_linear_preds <- predict(final_model, type = "linear")

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

  # Store baseline model in the survreg model object
  final_model$baseline_model <- baseline_info$model
  final_model$baseline_sf <- baseline_info$sf

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific parametric survival model fitting complete.\n")

  result <- list(
    survreg_model = final_model,
    times = sort(unique(XYTrain_cs[[timevar]][XYTrain_cs$status_cs == 1])),
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
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @export
Predict_CRModel_SurvReg <- function(modelout, newdata, newtimes = NULL) {

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
  # Get linear predictors from the parametric model
  linear_preds <- stats::predict(modelout$survreg_model,
                                 newdata = newdata_prepared,
                                 type = "linear")

  # Use the stored baseline hazard from the training fit
  # and the new scores (linear_preds) to get survival probabilities for the new data
  sf <- survival::survfit(
    modelout$survreg_model$baseline_model, # Use the stored Cox model
    newdata = data.frame("score" = linear_preds),
    conf.int = .95
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