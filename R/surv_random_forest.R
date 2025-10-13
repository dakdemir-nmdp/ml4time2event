#' @title SurvModel_RF
#'
#' @description Fit a RF model for survival outcomes using randomForestSRC.
#' Includes tuning of nodesize and mtry.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param ntree integer value, number of trees to grow (default: 300)
#' @param samplesize integer value, sample size for each grown tree (default: 500)
#' @param nsplit integer value, maximum number of splits for each tree (default: 5)
#' @param trace logical, trace tuning process or not (default: TRUE)
#' @param splitrule character, split rule for trees (default: "bs.gradient")
#' @param nodesize_try numeric vector, nodesize values to try during tuning (default: c(1, 5, 10, 15))
#' @param ... additional parameters passed to randomForestSRC functions
#'
#' @return a list of four items: model: fitted randomForestSRC model object,
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables,
#'  expvars: the explanatory variables used.
#'
#' @importFrom randomForestSRC tune rfsrc predict.rfsrc
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @export
SurvModel_RF<-function(data, expvars, timevar, eventvar, ntree=300, samplesize=500, nsplit=5, trace=TRUE, 
                       splitrule="bs.gradient", nodesize_try=c(1, 5, 10, 15), ...){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Convert character columns to factors for randomForestSRC
  for (vari in expvars){
    if (is.character(data[[vari]])){
      data[[vari]]<-as.factor(data[[vari]])
    }
  }

  # Define formula
  formRF<-stats::as.formula(paste("Surv(",timevar, ",", eventvar,") ~ .", collapse = "")) # Removed survival::

  # Adjust samplesize if it exceeds 70% of data
  samplesize <- min(ceiling(0.7 * nrow(data)), samplesize)

  # Tune hyperparameters (nodesize, mtry) with user-provided parameters
  o <- randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                             splitrule = splitrule, samptype = "swor", sampsize = samplesize,
                             trace = trace, nsplit = nsplit, stepFactor = 1.5,
                             mtryStart = 2, # Start tuning mtry from 2
                             nodesizeTry = nodesize_try, # Use user-provided nodesize values
                             ntreeTry = ntree, # Use fixed ntree for tuning speed
                             ...)

  # Fit final model with optimal parameters
  nodesize_opt <- if (!is.null(names(o$optimal)) && "nodesize" %in% names(o$optimal)) o$optimal[["nodesize"]] else o$optimal[[1]]
  mtry_opt <- if (!is.null(names(o$optimal)) && "mtry" %in% names(o$optimal)) o$optimal[["mtry"]] else o$optimal[[2]]

  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                                   nodesize = nodesize_opt, ntree = ntree, mtry = mtry_opt,
                                   tree.err = FALSE, importance = TRUE, statistics = TRUE,
                                   do.trace = trace, splitrule = splitrule, samptype = "swor",
                                   sampsize = samplesize, nsplit = nsplit, ...)

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] == 1, timevar]))

  result <- list(model = hd.obj, times = times, varprof = varprof, expvars = expvars)
  class(result) <- c("ml4t2e_surv_rf", "list")
  return(result)
}

#' @title Predict_SurvModel_RF
#'
#' @description Get predictions from a RF survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_RF' (a list containing 'model', 'times', 'varprof', 'expvars')
#' @param newdata the data for which the predictions are to be calculated
#' @param newtimes optional vector of new time points for interpolation. If NULL, uses model's native time points.
#'
#' @return a list containing the following items:
#'  Probs: predicted Survival probability matrix (rows=times, cols=observations),
#'  Times: The times at which the probabilities are predicted.
#'
#' @importFrom randomForestSRC predict.rfsrc
#' @export
Predict_SurvModel_RF <- function(modelout, newdata, newtimes = NULL) {
  # ============================================================================
  # Input Validation
  # ============================================================================
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

  # Check for NAs and count complete cases
  has_missing <- any(is.na(newdata_prepared))
  n_complete_cases <- sum(complete.cases(newdata_prepared))
  n_expected_obs <- nrow(newdata_prepared)

  if (has_missing) {
    warning(sprintf(
      "Missing values detected in %d observation(s). ",
      n_expected_obs - n_complete_cases),
      "randomForestSRC::predict.rfsrc will drop incomplete cases."
    )
  }

  # ============================================================================
  # Make Predictions
  # ============================================================================
  predSurvsTestRF <- randomForestSRC::predict.rfsrc(modelout$model, newdata = newdata_prepared)

  # ============================================================================
  # Extract and Format Predictions
  # ============================================================================
  # Extract survival probabilities and times
  # randomForestSRC returns survival as obs x times matrix
  # We need to transpose to get times x obs and add time 0 with probability 1
  Probs <- t(predSurvsTestRF$survival)  # Transpose to times x obs
  Probs <- rbind(1, Probs)  # Add time 0 with probability 1 for all observations
  Times <- c(0, predSurvsTestRF$time.interest)

  # Validate dimensions
  n_returned_obs <- ncol(Probs)

  if (n_returned_obs != n_expected_obs) {
    # Check if the mismatch is due to missing values (expected)
    if (has_missing && n_returned_obs == n_complete_cases) {
      # This is expected - predict.rfsrc dropped incomplete cases
      # Expand predictions to match original newdata by filling with NAs
      full_pred_probs <- matrix(NA, nrow = nrow(Probs), ncol = n_expected_obs)
      complete_idx <- which(complete.cases(newdata_prepared))
      full_pred_probs[, complete_idx] <- Probs
      Probs <- full_pred_probs
    } else {
      # Unexpected mismatch
      stop(sprintf(
        "Dimension mismatch: predict.rfsrc returned %d observation(s) but expected %d.\n",
        n_returned_obs, n_expected_obs),
        sprintf("Complete cases: %d. ", n_complete_cases),
        "This indicates an unexpected issue in predict.rfsrc call or data preparation."
      )
    }
  }

  # If newtimes specified, interpolate to those times
  if (!is.null(newtimes)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, newtimes = newtimes)
    Times <- newtimes
  }

  return(list(
    Probs = Probs, # Return as rows=times, cols=observations
    Times = Times
  ))
}
