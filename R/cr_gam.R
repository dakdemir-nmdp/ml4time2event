#' @title score2proba
#' @description internal function: from linear score to survival probabilities using Cox PH baseline.
#' @param datasurv data frame with 'time' and 'event' columns.
#' @param score numeric vector of linear predictor scores.
#' @param conf.int confidence level for survival curve.
#' @param which.est which estimate to return ("point", "lower", "upper").
#' @return a list containing the Cox model ('model') and survfit object ('sf').
#' @importFrom survival coxph Surv survfit coxph.control
#' @noRd
score2proba <-
  function(datasurv, score, conf.int=0.95, which.est=c("point", "lower", "upper")) {
    which.est <- match.arg(which.est)
    # pred <- rep(NA, length(score)) # Not used
    # names(pred) <- names(score) # Not used
    datacox<-cbind(datasurv, data.frame(score=score))
    # Ensure column names are 'time' and 'event' for Surv formula
    colnames(datacox)[1:2]<-c("time","event")
    # Fit a Cox model with score as the only predictor, fixing coefficient to 1 (offset)
    # This estimates the baseline hazard based on the provided scores
    predm <- survival::coxph(survival::Surv(time, event) ~ score, data=datacox, init=1, control=survival::coxph.control(iter.max = 0))
    # Get survival curve predictions based on the fitted baseline hazard and new scores
    sf <- survival::survfit(predm, newdata=data.frame("score"=score), conf.int=conf.int)
    return(list(model=predm,sf=sf))
  }

#' @title CRModel_GAM
#'
#' @description Fit a GAM model for competing risks outcomes using cause-specific modeling.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param event_codes character or numeric vector identifying the event code(s) to
#'   model. GAM competing risks can fit multiple causes simultaneously. If NULL
#'   (default), all non-zero event codes observed in the data are used. The first
#'   entry defines the default event of interest.
#' @param shrinkTreshold integer value, minimum number of factor levels for factor variables to be considered for shrinkage ('re' basis).
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#' @param event_of_interest optional character or numeric scalar indicating a specific event code
#'   that should be prioritized as the primary event of interest. If provided, this
#'   event code must be one of the codes specified in 'event_codes'.
#'
#' @return a list with the following components:
#'   \item{gam_model}{the fitted cause-specific GAM model object from mgcv::gam}
#'   \item{times}{vector of unique event times in the training data for the event of interest}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_gam"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{event_codes}{character vector of event codes included in the model}
#'   \item{event_codes_numeric}{numeric vector of event codes included}
#'   \item{default_event_code}{character scalar for the default event code}
#'   \\item{default_event_code_numeric}{numeric scalar for the default event code}
#'   \item{time_range}{vector with min and max observed event times}
#'
#' @importFrom mgcv gam cox.ph s
#' @importFrom stats as.formula predict
#' @export
CRModel_GAM <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                        shrinkTreshold = 10, ntimes = 50, verbose = FALSE, event_of_interest = NULL) {

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
    event_codes <- available_events
  }

  event_codes <- as.character(event_codes)

  # Validate event_of_interest AFTER event_codes is determined
  if (!is.null(event_of_interest) && !as.character(event_of_interest) %in% event_codes) {
    stop("'event_of_interest' must be one of the specified event_codes")
  }

  missing_event_codes <- setdiff(event_codes, available_events)
  if (length(missing_event_codes) > 0) {
    stop("The following event_codes are not present in the data: ",
         paste(missing_event_codes, collapse = ", "))
  }

  event_codes_numeric <- suppressWarnings(as.numeric(event_codes))
  if (any(is.na(event_codes_numeric))) {
    stop("GAM competing risks requires numeric event codes. Unable to coerce: ",
         paste(event_codes[is.na(event_codes_numeric)], collapse = ", "))
  }

  primary_event_code <- event_codes[1]
  primary_event_numeric <- event_codes_numeric[1]

  # If event_of_interest is provided, prioritize it as the primary event code
  if (!is.null(event_of_interest)) {
    primary_event_code <- as.character(event_of_interest)
    primary_event_numeric <- as.numeric(event_of_interest)
  }

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == primary_event_numeric]
  if (length(event_times) == 0) {
    stop("No events of type ", primary_event_code, " in training data. Cannot fit competing risks model.")
  }

  # Store event time range for reference
  time_range <- range(c(0, XYTrain[[timevar]][XYTrain[[eventvar]] %in% event_codes_numeric]), na.rm = TRUE)

  # ============================================================================
  # Model Fitting - Cause-Specific GAM Models for ALL Competing Events
  # ============================================================================
  if (verbose) {
    cat("Fitting cause-specific GAM models for event codes:", paste(event_codes_numeric, collapse = ", "), "\n")
  }

  # Store models for requested event codes
  gam_models_all_causes <- vector("list", length(event_codes_numeric))
  names(gam_models_all_causes) <- as.character(event_codes_numeric)

  # Fit a separate GAM model for each event type
  for (cause in event_codes_numeric) {
    if (verbose) cat("Fitting GAM model for event type", cause, "...\n")

    # Create cause-specific data: event = 1 if this cause, 0 otherwise (censored or competing)
    XYTrain_cause <- XYTrain
    XYTrain_cause$status_cs <- ifelse(XYTrain[[eventvar]] == cause, 1, 0)

    if (sum(XYTrain_cause$status_cs) < 5) {
      warning("Fewer than 5 events of type ", cause, ". Skipping this cause.")
      gam_models_all_causes[[as.character(cause)]] <- NULL
      next
    }

    # Identify variable types for GAM formula construction
    numvars <- expvars[which(sapply(as.data.frame(XYTrain_cause[, expvars, drop=FALSE]), is.numeric))]
    fctvars <- expvars[which(sapply(as.data.frame(XYTrain_cause[, expvars, drop=FALSE]), function(x) {
      (is.factor(x) | is.character(x))
    }))]

    # Convert characters to factors
    for (i in fctvars) {
      if (is.character(XYTrain_cause[[i]])) {
        XYTrain_cause[[i]] <- as.factor(XYTrain_cause[[i]])
      }
    }

    # Separate variables based on number of levels for smoothing/shrinkage
    catvarstoshrink <- if (length(fctvars) > 0) {
      fctvars[which(sapply(as.data.frame(XYTrain_cause[, fctvars, drop=FALSE]), function(x) {
        length(levels(x)) > shrinkTreshold
      }))]
    } else {
      c()
    }
    catvarsnottoshrink <- if (length(fctvars) > 0) {
      fctvars[which(sapply(as.data.frame(XYTrain_cause[, fctvars, drop=FALSE]), function(x) {
        length(levels(x)) <= shrinkTreshold
      }))]
    } else {
      c()
    }
    numvarstosmooth <- if (length(numvars) > 0) {
      numvars[which(sapply(as.data.frame(XYTrain_cause[, numvars, drop=FALSE]), function(x) {
        length(unique(x[!is.na(x)])) > shrinkTreshold
      }))]
    } else {
      c()
    }
    numvarsnottosmooth <- if (length(numvars) > 0) {
      numvars[which(sapply(as.data.frame(XYTrain_cause[, numvars, drop=FALSE]), function(x) {
        length(unique(x[!is.na(x)])) <= shrinkTreshold
      }))]
    } else {
      c()
    }

    # Build GAM formula string
    formGAM_terms <- c()
    for (vari in catvarstoshrink) {
      formGAM_terms <- c(formGAM_terms, paste0("s(", vari, ", bs='re')")) # Random effect smooth for high-cardinality factors
    }
    for (vari in numvarstosmooth) {
      formGAM_terms <- c(formGAM_terms, paste0("s(", vari, ")")) # Default smooth for numerics
    }
    # Add linear terms for remaining variables
    formGAM_terms <- c(formGAM_terms, numvarsnottosmooth, catvarsnottoshrink)

    # Combine terms into formula
    if (length(formGAM_terms) > 0) {
      formGAM_rhs <- paste(formGAM_terms, collapse = "+")
    } else {
      formGAM_rhs <- "1" # Intercept only if no predictors
    }

    formGAM <- stats::as.formula(
      paste0("survival::Surv(", timevar, ", status_cs) ~ ", formGAM_rhs)
    )
    if (verbose) print(formGAM)

    # Fit the GAM model with Cox PH family for cause-specific modeling
    gam_model_cause <- tryCatch(
      mgcv::gam(
        formGAM,
        family = mgcv::cox.ph(),
        data = XYTrain_cause,
        select = TRUE # Enable shrinkage via double penalty approach
      ),
      error = function(e) {
        warning("Failed to fit GAM model for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(gam_model_cause)) {
      next
    }

    # Create baseline model for prediction using score2proba approach
    # Get linear predictors for training data
    train_linear_preds <- stats::predict(gam_model_cause,
                                         newdata = XYTrain_cause,
                                         type = "link")

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

    # Store baseline model in the GAM model object
    gam_model_cause$baseline_model <- baseline_info$model
    gam_model_cause$baseline_sf <- baseline_info$sf

    gam_models_all_causes[[as.character(cause)]] <- gam_model_cause
  }

  # The main model for the event of interest
  gam_model <- gam_models_all_causes[[as.character(primary_event_numeric)]]

  if (is.null(gam_model)) {
    stop("Failed to fit GAM model for the event of interest (event_code = ", primary_event_code, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific GAM model fitting complete.\n")

  result <- list(
    gam_model = gam_model,
    gam_models_all_causes = gam_models_all_causes,  # All cause-specific models for Aalen-Johansen
    event_codes = event_codes,
    event_codes_numeric = event_codes_numeric,
    default_event_code = primary_event_code,
    default_event_code_numeric = primary_event_numeric,
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_gam",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    time_range = time_range
  )

  class(result) <- c("ml4t2e_cr_gam", "list")
  return(result)
}

#' @title Predict_CRModel_GAM
#'
#' @description Get predictions from a fitted cause-specific GAM competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_GAM'
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the times from the training data.
#'   Can be any positive values - interpolation handles all time points.
#' @param event_of_interest character or numeric scalar indicating which event code
#'   to predict. If NULL (default), uses the event code stored during training.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @export
Predict_CRModel_GAM <- function(modelout, newdata, newtimes = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_gam")) {
    stop("'modelout' must be output from CRModel_GAM")
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
    event_of_interest <- modelout$default_event_code
  }

  event_of_interest <- as.character(event_of_interest)
  event_idx <- match(event_of_interest, modelout$event_codes)

  if (is.na(event_idx)) {
    stop("event_of_interest ", event_of_interest, " was not present in training data. Available event codes: ",
         paste(modelout$event_codes, collapse = ", "))
  }

  target_event_numeric <- modelout$event_codes_numeric[event_idx]

  # Generate default times if not specified
  use_native_times <- is.null(newtimes)
  if (!use_native_times) {
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
  # Make Predictions using Aalen-Johansen Estimator
  # ============================================================================
  # Get survival predictions from ALL cause-specific GAM models
  cause_specific_survs <- vector("list", length(modelout$event_codes_numeric))
  names(cause_specific_survs) <- as.character(modelout$event_codes_numeric)

  # Predict survival for each cause
  surv_times <- NULL
  for (cause in modelout$event_codes_numeric) {
    cause_char <- as.character(cause)

    if (is.null(modelout$gam_models_all_causes[[cause_char]])) {
      warning("No model available for cause ", cause, ". Assuming S(t) = 1 for all times.")
      next
    }

    # Get linear predictors from the GAM model
    linear_preds <- tryCatch(
      stats::predict(modelout$gam_models_all_causes[[cause_char]],
                     newdata = newdata_prepared,
                     type = "link"),
      error = function(e) {
        warning("Prediction failed for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(linear_preds)) {
      next
    }

    # Use the stored baseline hazard to get survival probabilities
    sf <- tryCatch(
      survival::survfit(
        modelout$gam_models_all_causes[[cause_char]]$baseline_model,
        newdata = data.frame("score" = linear_preds),
        conf.int = 0.95
      ),
      error = function(e) {
        warning("survfit failed for cause ", cause, ": ", e$message)
        NULL
      }
    )

    if (is.null(sf)) {
      next
    }

    if (is.null(surv_times)) {
      surv_times <- sf$time
    }

    # Extract survival probabilities
    surv_probs_cause <- sf$surv

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

  if (is.null(surv_times)) {
    stop("Unable to derive survival time grid from cause-specific models")
  }

  # Use proper Aalen-Johansen approach by calculating overall survival first  
  # Similar to XGBoost approach since GAM also uses baseline Cox models
  n_obs <- nrow(newdata_prepared)
  n_times <- length(surv_times)
  
  # Calculate cumulative hazards for each cause and each observation
  cum_hazards_all_causes <- array(0, dim = c(n_times, n_obs, length(modelout$event_codes_numeric)))
  
  for (i in seq_along(modelout$event_codes_numeric)) {
    cause <- modelout$event_codes_numeric[i]
    cause_char <- as.character(cause)

    if (!is.null(modelout$gam_models_all_causes[[cause_char]])) {
      # Get linear predictors from GAM
      linear_preds <- tryCatch(
        as.vector(stats::predict(
          modelout$gam_models_all_causes[[cause_char]],
          newdata = newdata_prepared,
          type = "link"
        )),
        error = function(e) {
          warning("GAM prediction failed for cause ", cause, ": ", e$message)
          rep(0, n_obs)
        }
      )

      # Get baseline cumulative hazard from survfit
      sf_baseline <- survival::survfit(
        modelout$gam_models_all_causes[[cause_char]]$baseline_model,
        newdata = data.frame("score" = rep(0, 1))
      )

      # Check if sf_baseline has valid data
      if (length(sf_baseline$time) == 0 || any(is.na(sf_baseline$surv))) {
        # If no valid baseline data, assume no events for this cause
        baseline_cum_hazard <- rep(0, n_times)
      } else {
        # Interpolate baseline cumulative hazard to our time grid
        baseline_cum_hazard <- stats::approx(
          x = c(0, sf_baseline$time),
          y = c(0, -log(pmax(sf_baseline$surv, 1e-10))),
          xout = surv_times,
          method = "constant",
          f = 0,
          rule = 2
        )$y
      }

      # Apply individual risk factors: Λ_j(t|x) = Λ_0j(t) * exp(β*x)
      for (j in seq_len(n_obs)) {
        cum_hazards_all_causes[, j, i] <- baseline_cum_hazard * exp(linear_preds[j])
      }
    }
  }
  
  # Step 2: Calculate overall survival S(t) = exp(-Σ_j Λ_j(t))
  overall_cum_hazard <- apply(cum_hazards_all_causes, c(1, 2), sum)
  overall_survival <- exp(-overall_cum_hazard)
  
  # Step 3: Calculate CIF using Aalen-Johansen formula
  cif_matrix <- matrix(0, nrow = n_times, ncol = n_obs)
  event_idx_numeric <- which(modelout$event_codes_numeric == target_event_numeric)

  if (length(event_idx_numeric) > 0) {
    target_cum_hazards <- cum_hazards_all_causes[, , event_idx_numeric]
    
    for (j in seq_len(n_obs)) {
      for (t in seq_len(n_times)) {
        if (t == 1) {
          # For first time point, CIF = hazard increment
          target_hazard_increment <- target_cum_hazards[t, j] 
          cif_matrix[t, j] <- 1.0 * target_hazard_increment
        } else {
          # CIF(t) = CIF(t-1) + S(t-1) * Δλ_j(t)
          target_hazard_increment <- target_cum_hazards[t, j] - target_cum_hazards[t-1, j]
          cif_matrix[t, j] <- cif_matrix[t-1, j] + overall_survival[t-1, j] * target_hazard_increment
        }
      }
    }
  }
  
  # Ensure bounds and monotonicity
  # Use matrix() to preserve dimensions after pmax/pmin operations
  original_dims <- dim(cif_matrix)
  cif_matrix <- matrix(pmax(0, pmin(as.vector(cif_matrix), 1)), 
                       nrow = original_dims[1], ncol = original_dims[2])
  for (j in 1:n_obs) {
    cif_matrix[, j] <- cummax(cif_matrix[, j])
  }

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  # Ensure time zero included with zero CIF
  if (!any(abs(surv_times) < .Machine$double.eps)) {
    surv_times <- c(0, surv_times)
    cif_matrix <- rbind(rep(0, n_obs), cif_matrix)
  } else {
    zero_idx <- which.min(abs(surv_times))
    cif_matrix[zero_idx, ] <- 0
  }

  if (use_native_times) {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times, observations]
    result_times <- surv_times
  } else {
    # Interpolate to new time points
    # Use the standard CIF interpolation utility function
    pred_cifs <- cifMatInterpolaltor(
      probsMat = t(cif_matrix),  # cifMatInterpolaltor expects [observations, times]
      times = surv_times,
      newtimes = newtimes
    )

    # cifMatInterpolaltor returns [newtimes, observations], keep as [times, observations]
    result_cifs <- pred_cifs
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
