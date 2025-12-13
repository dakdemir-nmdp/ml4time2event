#'
#'
#'
#'
CRModel_RF <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                      ntree = 300, samplesize = 500, nsplit = 5, trace = FALSE,
                      verbose = NULL, splitrule = "logrankCR", nodesize_try = c(5, 10, 15), ...) {
  # Map verbose to trace for backward compatibility and consistency
  if (!is.null(verbose)) {
    trace <- verbose
  }

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
  if (!is.null(event_codes) && length(event_codes) == 0) {
    stop("'event_codes' must be NULL or a non-empty vector")
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

  available_events <- sort(unique(data[[eventvar]][data[[eventvar]] != 0]))
  if (length(available_events) == 0) {
    stop("No events found in the training data.")
  }

  if (is.null(event_codes)) {
    event_codes <- as.character(available_events[1])
  }

  event_codes <- as.character(event_codes)

  if (length(event_codes) != 1) {
    stop("CRModel_RF supports exactly one event code. Received ", length(event_codes), ".")
  }

  if (!event_codes %in% as.character(available_events)) {
    stop("Requested event code ", event_codes, " not present in training data. Available codes: ",
         paste(as.character(available_events), collapse = ", "))
  }

  event_code_numeric <- suppressWarnings(as.numeric(event_codes))
  if (is.na(event_code_numeric)) {
    stop("Random forest competing risks requires numeric event codes. Unable to coerce '",
         event_codes, "' to numeric.")
  }

  # Define formula for randomForestSRC (competing risks)
  formRF <- stats::as.formula(paste("Surv(", timevar, ",", eventvar, ") ~ .", collapse = ""))

  # Adjust samplesize if it exceeds 70% of data
  samplesize <- min(ceiling(0.7 * nrow(data)), samplesize)

  # ============================================================================
  # Model Fitting with Tuning
  # ============================================================================
  # Tune hyperparameters (nodesize, mtry) with user-provided parameters
  # Wrap in tryCatch to handle tuning failures gracefully
  tune_result <- tryCatch({
    randomForestSRC::tune(formRF, data = data[, c(timevar, eventvar, expvars), drop = FALSE],
                          splitrule = splitrule, samptype = "swor", sampsize = samplesize,
                          trace = trace, nsplit = nsplit, stepFactor = 1.5,
                          mtryStart = 2, # Start tuning mtry from 2
                          nodesizeTry = nodesize_try, # Use user-provided nodesize values
                          ntreeTry = ntree, # Use fixed ntree for tuning speed
                          cause = c(event_code_numeric,
                                     setdiff(unique(data[[eventvar]][data[[eventvar]] != 0]), event_code_numeric)),
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

  rf_model <- randomForestSRC::rfsrc(formRF, data = data[, c(timevar, eventvar, expvars), drop = FALSE],
                                     nodesize = nodesize_opt, ntree = ntree, mtry = mtry_opt,
                                     tree.err = FALSE, importance = TRUE, statistics = TRUE,
                                     do.trace = trace, splitrule = splitrule, samptype = "swor",
                                     sampsize = samplesize, nsplit = nsplit,
                                     cause = c(event_code_numeric,
                                               setdiff(unique(data[[eventvar]][data[[eventvar]] != 0]), event_code_numeric)),
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
    event_codes = event_codes,
    event_code_numeric = event_code_numeric,
    time_range = time_range
  )

  class(result) <- c("ml4t2e_cr_rf", "CRModel_RF")
  return(result)
}




#'
#'
#'
#'
Predict_CRModel_RF <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("'modelout' is missing")
  }
  if (!is.list(modelout) || !all(c("expvars", "event_codes") %in% names(modelout))) {
    stop("'modelout' must be output from CRModel_RF")
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
    event_of_interest <- modelout$event_codes
  }

  event_of_interest <- as.character(event_of_interest)

  if (length(event_of_interest) != 1) {
    stop("Random forest competing risks predictions require a single event code")
  }

  if (!identical(event_of_interest, modelout$event_codes)) {
    stop("RF models can only predict for the event they were trained on (event code = ",
         modelout$event_codes, "). Requested event code: ", event_of_interest)
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
  cif_matrix <- pred_rf$cif[, , paste0("CIF.", modelout$event_code_numeric)]

  # Convert to time-by-observation orientation and add time 0 with CIF = 0
  cif_time_obs <- t(cif_matrix)
  cif_time_obs <- rbind(rep(0, ncol(cif_time_obs)), cif_time_obs)
  times_with_t0 <- c(0, pred_rf$time.interest)

  # ============================================================================
  # Apply Interpolation if needed
  # ============================================================================
  if (is.null(new_times)) {
    result_cifs <- cif_time_obs
    result_times <- times_with_t0
  } else {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    new_times <- sort(unique(new_times))

    result_cifs <- cifMatInterpolator(
      probsMat = cif_time_obs,
      times = times_with_t0,
      new_times = new_times
    )
    result_times <- new_times
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
