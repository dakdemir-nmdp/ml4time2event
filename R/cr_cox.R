#' @importFrom survival coxph Surv
NULL

#'
#'
#'
#'
#'
#'
#'
CRModel_Cox <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                       varsel = "none", penalty = "AIC",
                       alpha = 0.5, nfolds = 10,
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
    stop("`event_codes` must be NULL or a non-empty vector if provided")
  }
  varsel <- match.arg(varsel, c("none", "backward", "forward", "both", "penalized"))
  penalty <- match.arg(penalty, c("AIC", "BIC"))
  
  # Validate alpha and nfolds parameters
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a numeric value between 0 and 1")
  }
  if (!is.numeric(nfolds) || nfolds < 3) {
    stop("`nfolds` must be a numeric value >= 3")
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

  # Get unique event types, excluding censor code (0)
  unique_events <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))
  if (length(unique_events) == 0) {
    stop("No events found in the training data.")
  }

  # If event_codes is not provided, use all found events
  if (is.null(event_codes)) {
    event_codes <- as.character(unique_events)
  } else {
    # Check if specified event_codes are in the data
    event_codes <- as.character(event_codes)
    missing_event_codes <- setdiff(event_codes, as.character(unique_events))
    if (length(missing_event_codes) > 0) {
      stop("The following event_codes are not present in the data: ",
           paste(missing_event_codes, collapse = ", "))
    }
  }

  # Check for sufficient events for each cause
  for (cause in event_codes) {
    event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == as.numeric(cause)]
    if (length(event_times) == 0) {
      stop(paste0("No events of type ", cause, " in training data. Cannot fit competing risks model."))
    }
  }

  time_range <- range(XYTrain[[timevar]][XYTrain[[eventvar]] %in% as.numeric(event_codes)], na.rm = TRUE)
  times <- seq(time_range[1], time_range[2], length.out = ntimes)

  # ============================================================================
  # Fit Cause-Specific Cox Models
  # ============================================================================
  if (verbose) cat("Fitting cause-specific Cox models...\n")

  cph_models_all_causes <- list()

  for (cause in event_codes) {
    cause_val <- as.numeric(cause)

    # Create a binary event indicator for the current cause
    event_indicator <- ifelse(XYTrain[[eventvar]] == cause_val, 1, 0)

    # Create a temporary dataset for this cause
    cause_data <- XYTrain
    cause_data$event_bin <- event_indicator

    # Formula for the Cox model
    formula_str <- paste0("Surv(", timevar, ", event_bin) ~ ", paste(expvars, collapse = " + "))
    formula <- as.formula(formula_str)

    # Fit the initial Cox model
    cph_model <- tryCatch({
      coxph(formula, data = cause_data, x = TRUE)
    }, error = function(e) {
      warning(sprintf("Failed to fit initial Cox model for cause %s: %s", cause, e$message))
      return(NULL)
    })

    if (is.null(cph_model)) {
      cph_models_all_causes[[cause]] <- NULL
      next
    }

    # Variable selection if specified
    if (varsel == "penalized") {
      if (verbose) cat(sprintf("Fitting penalized Cox model for cause %s (alpha=%.2f, nfolds=%d)...\n", cause, alpha, nfolds))
      
      # Prepare data for glmnet (matrix form)
      X <- stats::model.matrix(as.formula(paste0("~", paste(expvars, collapse = " + "))), data = cause_data)[, -1, drop = FALSE]
      times_cox <- cause_data[[timevar]]
      events_cox <- cause_data$event_bin
      
      # Fit penalized Cox model via glmnet
      cph_model <- tryCatch({
        glmnet::cv.glmnet(X, Surv(times_cox, events_cox), family = "cox", 
                         alpha = alpha, nfolds = nfolds, standardize = TRUE)
      }, error = function(e) {
        warning(sprintf("Penalized Cox fitting failed for cause %s: %s. Using standard Cox model.", cause, e$message))
        return(NULL)
      })
      
      if (is.null(cph_model)) {
        # Fallback to standard Cox model
        cph_model <- tryCatch({
          coxph(formula, data = cause_data, x = TRUE)
        }, error = function(e) {
          warning(sprintf("Fallback Cox model also failed for cause %s: %s", cause, e$message))
          return(NULL)
        })
      }
      
      if (!is.null(cph_model)) {
        cph_models_all_causes[[cause]] <- cph_model
      } else {
        cph_models_all_causes[[cause]] <- NULL
      }
      
    } else if (varsel != "none") {
      if (verbose) cat(sprintf("Performing %s selection for cause %s...\n", varsel, cause))

      # Define the scope for stepwise selection
      scope <- list(
        lower = as.formula(paste0("Surv(", timevar, ", event_bin) ~ 1")),
        upper = formula
      )

      # Determine direction for step function
      direction <- switch(varsel,
                          "backward" = "backward",
                          "forward" = "forward",
                          "both" = "both")

      # Penalty for step function
      k_penalty <- ifelse(penalty == "BIC", log(nrow(cause_data)), 2)

      # Run stepwise selection
      cph_model_selected <- tryCatch({
        stats::step(cph_model, scope = scope, direction = direction,
                    trace = 0, k = k_penalty)
      }, error = function(e) {
        warning(sprintf("Stepwise selection failed for cause %s: %s. Using initial model.", cause, e$message))
        return(cph_model) # Fallback to the initial model
      })

      cph_models_all_causes[[cause]] <- cph_model_selected
    } else {
      cph_models_all_causes[[cause]] <- cph_model
    }
  }

  # ============================================================================
  # Finalize and Return Output
  # ============================================================================

  # Check if any model failed
  failed_models <- names(cph_models_all_causes)[sapply(cph_models_all_causes, is.null)]
  if (length(failed_models) > 0) {
    warning("Failed to fit Cox models for the following causes: ", paste(failed_models, collapse = ", "))
  }

  # Extract final variables from the models
  final_expvars <- unique(unlist(lapply(cph_models_all_causes, function(m) {
    if (!is.null(m)) labels(stats::terms(m)) else NULL
  })))

  if (is.null(final_expvars)) final_expvars <- character(0)


  output <- list(
    cph_models_all_causes = cph_models_all_causes,
    times = times,
    varprof = varprof,
    model_type = "cr_cox",
    expvars = final_expvars,
    timevar = timevar,
    eventvar = eventvar,
    event_codes = event_codes,
    varsel_method = varsel,
    time_range = time_range,
    alpha = alpha,
    nfolds = nfolds
  )

  class(output) <- c("ml4t2e_cr_cox", "CRModel_Cox")

  return(output)
}

#'
#'
#'
#'
#'
#'
#'
Predict_CRModel_Cox <- function(modelout, newdata, new_times = NULL, event_of_interest = NULL) {

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

  # Handle event_of_interest parameter
  if (is.null(event_of_interest)) {
    if (length(modelout$event_codes) == 0) {
      stop("Model does not have any event codes stored")
    }
    event_of_interest <- modelout$event_codes[1]
  }

  event_of_interest <- as.character(event_of_interest)

  if (length(event_of_interest) != 1) {
    stop("'event_of_interest' must be a single event code")
  }

  if (!event_of_interest %in% modelout$event_codes) {
    stop("Requested event_of_interest ", event_of_interest, " was not present in training data. Available event codes: ",
         paste(modelout$event_codes, collapse = ", "))
  }

  # Generate default times if not specified
  use_default_times <- is.null(new_times)
  if (!use_default_times) {
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
  # Make Predictions using Aalen-Johansen Estimator
  # ============================================================================
  # Use proper Aalen-Johansen approach with cause-specific hazards
  base_times <- sort(unique(c(0, modelout$times)))
  cif_matrix <- aalenJohansenFromCoxModels(
    cox_models = modelout$cph_models_all_causes,
    newdata = newdata_prepared,
    times = base_times,
    event_of_interest = event_of_interest
  )

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (use_default_times) {
    # Return results at base times (including 0)
    result_cifs <- cif_matrix
    result_times <- base_times
  } else {
    # Interpolate to new time points
    # Use the standard CIF interpolation utility function
    result_cifs <- cifMatInterpolator(
      probsMat = cif_matrix,
      times = base_times,
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
