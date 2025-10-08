#' @title CRModel_RF
#'
#' @description Fit a random forest model for competing risks outcomes using randomForestSRC.
#' Includes tuning of nodesize and mtry parameters.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2 where 0=censored, 1=event of interest, 2=competing event)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param ntree integer value, number of trees to grow (default: 300)
#' @param samplesize integer value, sample size for each grown tree (default: 500)
#' @param nsplit integer value, maximum number of splits for each tree (default: 5)
#' @param trace logical, trace tuning process or not (default: TRUE)
#' @param splitrule character, split rule for trees (default: "logrankCR")
#' @param nodesize_try numeric vector, nodesize values to try during tuning (default: c(1, 5, 10, 15))
#' @param ... additional parameters passed to randomForestSRC functions
#'
#' @return a list with the following components:
#'   \item{rf_model}{the fitted randomForestSRC model object}
#'   \item{times}{vector of unique event times in the training data}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_rf"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'
#' @importFrom randomForestSRC tune rfsrc
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @export
CRModel_RF <- function(data, expvars, timevar, eventvar, failcode = 1,
                      ntree = 300, samplesize = 500, nsplit = 5, trace = TRUE,
                      splitrule = "logrankCR", nodesize_try = c(1, 5, 10, 15), ...) {

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
  # Create variable profile
  varprof <- VariableProfile(data, expvars)

  # Convert character columns to factors for randomForestSRC
  for (vari in expvars) {
    if (is.character(data[[vari]])) {
      data[[vari]] <- as.factor(data[[vari]])
    }
  }

  # Define formula for randomForestSRC (competing risks)
  formRF <- stats::as.formula(paste("Surv(", timevar, ",", eventvar, ") ~ .", collapse = ""))

  # Adjust samplesize if it exceeds 70% of data
  samplesize <- min(ceiling(0.7 * nrow(data)), samplesize)

  # ============================================================================
  # Model Fitting with Tuning
  # ============================================================================
  # Tune hyperparameters (nodesize, mtry) with user-provided parameters
  o <- randomForestSRC::tune(formRF, data = data[, c(timevar, eventvar, expvars), drop = FALSE],
                             splitrule = splitrule, samptype = "swor", sampsize = samplesize,
                             trace = trace, nsplit = nsplit, stepFactor = 1.5,
                             mtryStart = 2, # Start tuning mtry from 2
                             nodesizeTry = nodesize_try, # Use user-provided nodesize values
                             ntreeTry = ntree, # Use fixed ntree for tuning speed
                             cause = c(failcode, setdiff(unique(data[[eventvar]][data[[eventvar]] != 0]), failcode)),
                             ...)

  # Fit final model with optimal parameters
  rf_model <- randomForestSRC::rfsrc(formRF, data = data[, c(timevar, eventvar, expvars), drop = FALSE],
                                     nodesize = o$optimal[[1]], ntree = ntree, mtry = o$optimal[[2]],
                                     tree.err = FALSE, importance = TRUE, statistics = TRUE,
                                     do.trace = trace, splitrule = splitrule, samptype = "swor",
                                     sampsize = samplesize, nsplit = nsplit,
                                     cause = c(failcode, setdiff(unique(data[[eventvar]][data[[eventvar]] != 0]), failcode)),
                                     ...)

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] != 0, timevar]))

  # Get time range
  time_range <- range(data[data[[eventvar]] != 0, timevar])

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    rf_model = rf_model,
    times = times,
    varprof = varprof,
    model_type = "cr_rf",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
    time_range = time_range
  )

  class(result) <- "ml4t2e_cr_rf"
  return(result)
}




#' @title Predict_CRModel_RF
#'
#' @description Get predictions from a CR random forest model for a test dataset.
#'
#' @param modelout the output from 'CRModel_RF' (a list containing model and metadata)
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the model's native time points.
#'   Can be any positive values - interpolation handles all time points.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=observations, cols=times)}
#'   \item{Times}{the times at which CIFs are calculated
#'     (always includes time 0)}
#'
#' @importFrom randomForestSRC predict.rfsrc
#' @export
Predict_CRModel_RF <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_rf")) {
    stop("'modelout' must be output from CRModel_RF")
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

  # ============================================================================
  # Make Predictions
  # ============================================================================
  # Get predictions from randomForestSRC
  pred_rf <- randomForestSRC::predict.rfsrc(modelout$rf_model, newdata = newdata_prepared)

  # Extract CIF for the event of interest (failcode)
  # randomForestSRC returns cif as [observations, times, causes]
  cif_matrix <- pred_rf$cif[, , paste0("CIF.", modelout$failcode)]

  # Add time 0 with CIF = 0
  cif_with_t0 <- cbind(0, cif_matrix)  # Add column for time 0
  times_with_t0 <- c(0, pred_rf$time.interest)

  # ============================================================================
  # Apply Interpolation if needed
  # ============================================================================
  if (is.null(newtimes)) {
    # Use model's native time points: [times, observations]
    result_cifs <- t(cif_with_t0)  # Transpose to [times, observations]
    result_times <- times_with_t0
  } else {
    # Interpolate to new time points
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))

    # Use cifMatInterpolaltor for interpolation
    result_cifs <- cifMatInterpolaltor(
      probsMat = cif_with_t0,  # [observations, times]
      times = times_with_t0,
      newtimes = newtimes
    )
    # cifMatInterpolaltor returns [newtimes, observations], keep as [times, observations]
    result_cifs <- result_cifs
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
