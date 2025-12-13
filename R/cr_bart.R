#'
#'
#'
#'
CRModel_BART <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                         K = 10, ntree = 100, ndpost = 1000, nskip = 100,
                         keepevery = 10, verbose = FALSE) {
  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("Input 'expvars' must be a non-empty character vector.")
  }
  if (!timevar %in% colnames(data)) {
    stop("Input 'timevar' not found in data: ", timevar)
  }
  if (!eventvar %in% colnames(data)) {
    stop("Input 'eventvar' not found in data: ", eventvar)
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following 'expvars' not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (!is.null(event_codes) && length(event_codes) == 0) {
    stop("Input 'event_codes' must be NULL or a non-empty vector.")
  }

  # Identify available event codes (exclude censoring)
  available_events <- sort(unique(as.character(data[[eventvar]][data[[eventvar]] != 0])))
  if (length(available_events) == 0) {
    stop("No events found in the training data.")
  }

  if (is.null(event_codes)) {
    event_codes <- available_events[1]
  }

  event_codes <- as.character(event_codes)

  if (length(event_codes) != 1) {
    stop("CRModel_BART supports exactly one event code. Received ", length(event_codes), ".")
  }

  if (!event_codes %in% available_events) {
    stop(
      "Requested event code ", event_codes, " not present in training data. Available codes: ",
      paste(available_events, collapse = ", ")
    )
  }

  event_code_numeric <- suppressWarnings(as.numeric(event_codes))
  if (is.na(event_code_numeric)) {
    stop(
      "BART competing risks requires numeric event codes. Unable to coerce '",
      event_codes, "' to numeric."
    )
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================
  # Create variable profile
  varprof <- VariableProfile(data, expvars)

  # Prepare data
  times_train <- data[[timevar]]
  delta_train <- as.integer(data[[eventvar]])

  # Create model matrix
  x_train <- as.matrix(stats::model.matrix(~ -1 + ., data = data[, expvars, drop = FALSE]))

  # ============================================================================
  # Model Fitting
  # ============================================================================
  if (verbose) message("Fitting Competing Risks BART model...")
  # Fit BART model with retry logic for robustness, suppressing all output unless verbose
  failcount <- 0
  bart_model <- NULL

  while (is.null(bart_model) && failcount < 4) {
    failcount <- failcount + 1
    bart_model <- tryCatch(
      {
        if (verbose) {
          suppressMessages(BART::crisk.bart(
            x.train = x_train,
            times = times_train,
            delta = delta_train,
            x.test = x_train, # Include training data for predictions
            K = K,
            ntree = ntree,
            ndpost = ndpost,
            nskip = nskip,
            keepevery = keepevery,
            numcut = 2 # Small number of cuts for speed
          ))
        } else {
          bart_fit <- NULL
          invisible(capture.output({
            bart_fit <- suppressMessages(BART::crisk.bart(
              x.train = x_train,
              times = times_train,
              delta = delta_train,
              x.test = x_train, # Include training data for predictions
              K = K,
              ntree = ntree,
              ndpost = ndpost,
              nskip = nskip,
              keepevery = keepevery,
              numcut = 2 # Small number of cuts for speed
            ))
          }))
          bart_fit
        }
      },
      error = function(e) {
        if (verbose) message("BART fitting attempt ", failcount, " failed: ", e$message)
        NULL
      }
    )
  }

  if (is.null(bart_model)) {
    stop("Failed to fit BART model after ", failcount, " attempts")
  }

  # Get unique event times from training data for the modeled event code
  times <- sort(unique(data[data[[eventvar]] == event_code_numeric, timevar]))

  if (length(times) == 0) {
    stop("No events of type ", event_codes, " in training data. Cannot fit competing risks model.")
  }

  # Get time range restricted to modeled event codes
  time_range <- range(c(0, data[data[[eventvar]] %in% event_code_numeric, timevar]), na.rm = TRUE)

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    bart_model = bart_model,
    times = times,
    varprof = varprof,
    model_type = "cr_bart",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes,
    event_codes_numeric = event_code_numeric,
    default_event_code = event_codes,
    time_range = time_range,
    x_train = x_train,
    times_train = times_train,
    delta_train = delta_train
  )

  class(result) <- "ml4t2e_cr_bart"
  return(result)
}


#'
#'
#'
#'
Predict_CRModel_BART <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {
  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(newdata)) {
    stop("Input 'newdata' must be a data frame.")
  }

  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop("The following 'expvars' are missing in 'newdata': ", paste(missing_vars, collapse = ", "))
  }

  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    event_of_interest <- modelout$default_event_code
  }

  event_of_interest <- as.character(event_of_interest)

  if (length(event_of_interest) != 1) {
    stop("Input 'event_of_interest' must be a single event code.")
  }

  if (!identical(event_of_interest, modelout$default_event_code)) {
    stop(
      "BART models can only predict for the event they were trained on (event code = ",
      modelout$default_event_code, "). Requested event code: ", event_of_interest
    )
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
            warning(
              "Factor '", var, "' has new levels in newdata: ",
              paste(extra_levels, collapse = ", "),
              ". These will be set to NA."
            )
          }
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
            levels = training_levels
          )
        } else if (is.character(newdata_prepared[[var]])) {
          # Convert character to factor with training levels
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
            levels = training_levels
          )
        }
      }
    }
  }

  # ============================================================================
  # Make Predictions
  # ============================================================================
  # Create test model matrix
  x_test <- as.matrix(stats::model.matrix(~ -1 + ., data = newdata_prepared))

  # Prepare BART prediction structure
  pre <- BART::crisk.pre.bart(
    time = modelout$times_train,
    delta = modelout$delta_train,
    x.train = BART::bartModelMatrix(modelout$x_train),
    x.test = BART::bartModelMatrix(x_test),
    x.train2 = BART::bartModelMatrix(modelout$x_train),
    x.test2 = BART::bartModelMatrix(x_test),
    K = modelout$bart_model$K
  )

  # Predict using the fitted BART model
  pred <- predict(modelout$bart_model, pre$tx.test, pre$tx.test2)

  # Extract CIF predictions
  # pred$cif.test.mean is organized as [obs1_time1, obs1_time2, ..., obs1_timeK, obs2_time1, ...]
  # We need to reshape it to [observations, times]
  N <- nrow(x_test) # Number of actual test observations
  K <- modelout$bart_model$K # Number of time points

  total_len <- length(pred$cif.test.mean)
  expected_len <- N * K
  if (total_len != expected_len) {
    stop("Unexpected length of BART CIF predictions. Expected ", expected_len, " values, got ", total_len, ".")
  }

  # Reshape to matrix [observations, times]
  cif_matrix <- matrix(pred$cif.test.mean, nrow = N, ncol = K, byrow = TRUE)


  # Convert to time-by-observation orientation and add time 0 with CIF = 0
  cif_time_obs <- t(cif_matrix)
  cif_time_obs <- rbind(rep(0, ncol(cif_time_obs)), cif_time_obs)
  times_with_t0 <- c(0, modelout$bart_model$times)

  # ============================================================================
  # Apply Interpolation if needed
  # ============================================================================
  if (is.null(new_times)) {
    result_cifs <- cif_time_obs
    result_times <- times_with_t0
  } else {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("Input 'new_times' must be a numeric vector of non-negative values.")
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
