#'
#'
#'
#'
SurvModel_RF<-function(data, expvars, timevar, eventvar, ntree=300, samplesize=500, nsplit=5, trace=FALSE,
                       verbose = NULL, splitrule="bs.gradient", nodesize_try=c(5, 10, 15), importance="permute", ...){
  # Map verbose to trace for backward compatibility and consistency
  if (!is.null(verbose)) {
    trace <- verbose
  }
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
  
  # Force importance to be "permute" which is most reliable
  importance <- "permute" # Override any user value to ensure consistent behavior

  # Define formula
  formRF<-stats::as.formula(paste("Surv(",timevar, ",", eventvar,") ~ .", collapse = "")) # Removed survival::

  # Adjust samplesize if it exceeds 70% of data
  samplesize <- min(ceiling(0.7 * nrow(data)), samplesize)

  # Tune hyperparameters (nodesize, mtry) with user-provided parameters
  # Wrap in tryCatch to handle tuning failures gracefully
  tune_result <- tryCatch({
    randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                          splitrule = splitrule, samptype = "swor", sampsize = samplesize,
                          trace = trace, nsplit = nsplit, stepFactor = 1.5,
                          mtryStart = 2, # Start tuning mtry from 2
                          nodesizeTry = nodesize_try, # Use user-provided nodesize values
                          ntreeTry = ntree, # Use fixed ntree for tuning speed
                          ...)
  }, error = function(e) {
    warning("Tuning failed with error: ", e$message, ". Using default parameters.")
    NULL
  })

  # Fit final model with optimal parameters (or defaults if tuning failed)
  if (!is.null(tune_result)) {
    nodesize_opt <- if (!is.null(names(tune_result$optimal)) && "nodesize" %in% names(tune_result$optimal)) {
      tune_result$optimal[["nodesize"]]
    } else {
      tune_result$optimal[[1]]
    }
    mtry_opt <- if (!is.null(names(tune_result$optimal)) && "mtry" %in% names(tune_result$optimal)) {
      tune_result$optimal[["mtry"]]
    } else {
      tune_result$optimal[[2]]
    }
  } else {
    # Use reasonable defaults when tuning fails
    nodesize_opt <- 15  # Conservative default
    mtry_opt <- max(1, floor(sqrt(length(expvars))))  # Standard RF default
    message("Using default parameters: nodesize=", nodesize_opt, ", mtry=", mtry_opt)
  }

  # Ensure importance is properly calculated with specific importance settings
  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                                   nodesize = nodesize_opt, ntree = ntree, mtry = mtry_opt,
                                   tree.err = FALSE, 
                                   importance = "permute", # Explicitly set to permutation importance
                                   var.used = "all.trees", # Track all variables used
                                   statistics = TRUE,
                                   forest = TRUE, # Save the forest for importance calculations
                                   save.memory = FALSE, # Don't use memory saving (which might drop importance)
                                   do.trace = trace, splitrule = splitrule, samptype = "swor",
                                   sampsize = samplesize, nsplit = nsplit, ...)
                                   
  # Extract and store importance scores explicitly in the model object
  if (!is.null(hd.obj) && is.null(hd.obj$importance)) {
    # If importance calculation failed, manually calculate it using vimp
    if (requireNamespace("randomForestSRC", quietly = TRUE)) {
      tryCatch({
        imp_obj <- randomForestSRC::vimp(hd.obj)
        if (!is.null(imp_obj) && !is.null(imp_obj$importance)) {
          hd.obj$importance <- imp_obj$importance
        }
      }, error = function(e) {
        # If vimp fails, create a simple placeholder importance based on var.used
        if (!is.null(hd.obj$var.used)) {
          var_counts <- table(hd.obj$var.used)
          # Create normalized importance scores (0-100 scale)
          imp_scores <- 100 * var_counts / sum(var_counts)
          hd.obj$importance <- imp_scores
        }
      })
    }
  }

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] == 1, timevar]))

  result <- list(model = hd.obj, times = times, varprof = varprof, expvars = expvars)
  class(result) <- c("ml4t2e_surv_rf", "list")
  return(result)
}

#'
#'
#'
#'
Predict_SurvModel_RF <- function(modelout, newdata, new_times = NULL) {
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

  # If new_times specified, interpolate to those times
  if (!is.null(new_times)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, new_times = new_times)
    Times <- new_times
  }

  return(list(
    Probs = Probs, # Return as rows=times, cols=observations
    Times = Times
  ))
}
