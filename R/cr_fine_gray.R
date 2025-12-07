#'
#'
#'
#'
CRModel_FineGray <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                            ntimes = 50, verbose = FALSE) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("`expvars` must be a non-empty character vector")
  }
  if (!timevar %in% colnames(data)) {
    stop(paste0("`timevar` not found in data: ", timevar))
  }
  if (!eventvar %in% colnames(data)) {
    stop(paste0("`eventvar` not found in data: ", eventvar))
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste0("The following `expvars` not found in data: ", paste(missing_vars, collapse=", ")))
  }
  if (!is.null(event_codes) && length(event_codes) == 0) {
    stop("`event_codes` must be NULL or a non-empty vector")
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

  available_events <- sort(unique(as.character(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0])))
  if (length(available_events) == 0) {
    stop("No events found in the training data.")
  }

  if (is.null(event_codes)) {
    event_codes <- available_events[1]
  }

  event_codes <- as.character(event_codes)

  if (length(event_codes) != 1) {
    stop("`event_codes` must be a single value (one event code)")
  }
  if (!event_codes %in% available_events) {
    stop(paste0("`event_codes` ", event_codes, " not present in training data. Available codes: ", paste(available_events, collapse = ", ")))
  }
  failcode <- suppressWarnings(as.numeric(event_codes))
  if (is.na(failcode)) {
    stop(paste0("`event_codes` must be numeric or coercible to numeric. Unable to coerce '", event_codes, "' to numeric."))
  }

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == failcode]
  if (length(event_times) == 0) {
    stop(paste0("No events of type ", event_codes, " in training data. Cannot fit competing risks model."))
  }

  # Store event time range for reference
  time_range <- c(0, max(event_times))

  # ============================================================================
  # Model Fitting
  # ============================================================================
  if (verbose) cat("Fitting Fine-Gray model for event type", failcode, "...\n")

  # Create model matrix
  covmat <- stats::model.matrix(~ -1 + ., data = XYTrain[, expvars, drop = FALSE])

  # Scale covariates
  covmat_scaled <- scale(covmat, center = TRUE, scale = TRUE)
  meanTrain <- attr(covmat_scaled, "scaled:center")
  sdTrain <- attr(covmat_scaled, "scaled:scale")
  
  # Store column names for later validation
  names(meanTrain) <- colnames(covmat)
  names(sdTrain) <- colnames(covmat)

    # SVD for dimensionality reduction (keep top 20 components or fewer)
  svdcovmat <- svd(covmat_scaled)
  n_components <- min(c(20, ncol(covmat_scaled)))
  Feat <- as.data.frame((covmat_scaled %*% svdcovmat$v)[, 1:n_components])
  colnames(Feat) <- paste0("PC", 1:n_components)

  # Fit model with default elastic net parameters (lambda=0.01, alpha=0.5)
  # For simplicity, we'll use fixed parameters instead of grid search
  fg_model <- tryCatch({
    # Prepare covariate matrix as in temp_test_fastcmprsk.R
    cov <- as.matrix(Feat)
    model_data <- data.frame(
      ftime = XYTrain[[timevar]],
      fstatus = XYTrain[[eventvar]]
    )
    model_data$cov <- I(cov)  # I() prevents data.frame from splitting matrix into columns
    # Explicitly set failcode and cencode for Crisk
    fastcmprsk::fastCrrp(
      fastcmprsk::Crisk(ftime, fstatus, failcode = failcode, cencode = min(XYTrain[[eventvar]], na.rm = TRUE)) ~ cov,
      data = model_data,
      lambda = 0.01,
      alpha = 0.5,
      penalty = "ENET",
      standardize = FALSE,
      max.iter = 5000
    )
  },
  error = function(e) {
    stop("Failed to fit Fine-Gray model: ", e$message)
  })

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Fine-Gray model fitting complete.\n")

  result <- list(
    fg_model = fg_model,
    time_range = time_range,
    varprof = varprof,
    model_type = "fine_gray",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes,
    event_code_numeric = failcode,
    train_data = XYTrain,
    scaling = list(meanTrain = meanTrain, sdTrain = sdTrain),
    loadings = svdcovmat$v
  )

  class(result) <- c("ml4t2e_cr_finegray", "CRModel_FineGray")
  return(result)
}


#'
#'
#'
#'
#'
#'
#'
Predict_CRModel_FineGray <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame")
  }
  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop(paste0("The following variables are missing in `newdata`: ", paste(missing_vars, collapse = ", ")))
  }
  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    event_of_interest <- modelout$event_codes
  }
  event_of_interest <- as.character(event_of_interest)
  if (length(event_of_interest) != 1) {
    stop("`event_of_interest` must be a single value (one event code)")
  }
  if (!identical(event_of_interest, modelout$event_codes)) {
    stop(paste0("Fine-Gray models can only predict for the event they were trained on (event code = ",
                modelout$event_codes, "). Requested event code: ", event_of_interest))
  }
  failcode <- suppressWarnings(as.numeric(event_of_interest))
  if (is.na(failcode)) {
    stop(paste0("`event_of_interest` must be numeric or coercible to numeric. Unable to coerce '", event_of_interest, "' to numeric."))
  }

  # Generate default times if not specified
  if (!is.null(new_times)) {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    new_times <- sort(unique(new_times))
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
  # Create model matrix and apply same preprocessing as training
  covmat <- stats::model.matrix(~ -1 + ., data = newdata_prepared)

  # Ensure covmat has same columns as training (handle missing/extra factor levels)
  expected_cols <- names(modelout$scaling$meanTrain)
  if (!is.null(expected_cols)) {
    # Add missing columns with 0s
    missing_cols <- setdiff(expected_cols, colnames(covmat))
    if (length(missing_cols) > 0) {
      missing_matrix <- matrix(0, nrow = nrow(covmat), ncol = length(missing_cols))
      colnames(missing_matrix) <- missing_cols
      covmat <- cbind(covmat, missing_matrix)
    }
    # Remove extra columns and reorder to match training
    covmat <- covmat[, expected_cols, drop = FALSE]
  }

  # Apply scaling
  covmat_scaled <- scale(covmat,
                         center = modelout$scaling$meanTrain,
                         scale = modelout$scaling$sdTrain)

  # Apply SVD transformation
  n_components <- min(c(20, ncol(covmat_scaled)))
  Feat <- (covmat_scaled %*% modelout$loadings)[, 1:n_components, drop = FALSE]

  # Get baseline cumulative hazard from the model
  # The model provides Breslow jumps at specific times
  breslow <- modelout$fg_model$breslowJump
  if (is.null(breslow) || !is.matrix(breslow) || nrow(breslow) == 0) {
    baseline_times <- numeric(0)
    baseline_haz <- numeric(0)
  } else {
    baseline_times <- as.numeric(breslow[, 1])
    baseline_haz <- as.numeric(breslow[, 2])
  }

  # Compute CIF for each observation
  n_obs <- as.integer(nrow(Feat))
  cif_matrix <- matrix(NA_real_, nrow = as.integer(length(baseline_times) + 1L), ncol = n_obs)

  # Time 0 has CIF = 0
  cif_matrix[1, ] <- 0

  # Compute CIF at each baseline time point
  # Ensure coefficients are a numeric vector
  beta <- as.numeric(modelout$fg_model$coef)
  for (i in 1:n_obs) {
    # Linear predictor
    lp <- sum(Feat[i, ] * beta)

    # Cumulative hazard for this observation
    cumhaz <- exp(lp) * baseline_haz

    # Cumulative incidence function (subdistribution hazard)
    cif_vals <- 1 - exp(-cumsum(cumhaz))

    cif_matrix[-1, i] <- cif_vals
  }

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(new_times)) {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times+1, observations]
    result_times <- c(0, baseline_times)
  } else {
    # Interpolate to new time points
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    new_times <- sort(unique(new_times))

    result_cifs <- cifMatInterpolaltor(
      probsMat = cif_matrix,
      times = c(0, baseline_times),
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
