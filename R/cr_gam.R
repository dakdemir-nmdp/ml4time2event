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
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param shrinkTreshold integer value, minimum number of factor levels for factor variables to be considered for shrinkage ('re' basis).
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{gam_model}{the fitted cause-specific GAM model object from mgcv::gam}
#'   \item{times}{vector of unique event times in the training data for the event of interest}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_gam"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'
#' @importFrom mgcv gam cox.ph s
#' @importFrom stats as.formula predict
#' @export
CRModel_GAM <- function(data, expvars, timevar, eventvar, failcode = 1,
                       shrinkTreshold = 10, ntimes = 50, verbose = FALSE) {

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

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == failcode]
  if (length(event_times) == 0) {
    stop("No events of type ", failcode, " in training data. Cannot fit competing risks model.")
  }

  # Store event time range for reference
  time_range <- c(0, max(event_times))

  # ============================================================================
  # Model Fitting - Cause-Specific GAM Models for ALL Competing Events
  # ============================================================================
  # Identify all unique event types (excluding censoring = 0)
  all_event_types <- sort(unique(XYTrain[[eventvar]][XYTrain[[eventvar]] != 0]))

  if (verbose) {
    cat("Fitting cause-specific GAM models for all event types:", paste(all_event_types, collapse = ", "), "\n")
  }

  # Store models for all event types
  gam_models_all_causes <- vector("list", length(all_event_types))
  names(gam_models_all_causes) <- as.character(all_event_types)

  # Fit a separate GAM model for each event type
  for (cause in all_event_types) {
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

    formGAM <- stats::as.formula(paste(timevar, "~", formGAM_rhs))
    if (verbose) print(formGAM)

    # Fit the GAM model with Cox PH family for cause-specific modeling
    gam_model_cause <- tryCatch(
      mgcv::gam(
        formGAM,
        family = mgcv::cox.ph(),
        data = XYTrain_cause,
        weights = XYTrain_cause$status_cs, # Use cause-specific event indicator as weights
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
  gam_model <- gam_models_all_causes[[as.character(failcode)]]

  if (is.null(gam_model)) {
    stop("Failed to fit GAM model for the event of interest (failcode = ", failcode, ")")
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific GAM model fitting complete.\n")

  result <- list(
    gam_model = gam_model,
    gam_models_all_causes = gam_models_all_causes,  # All cause-specific models for Aalen-Johansen
    all_event_types = all_event_types,  # Event type codes
    times = sort(unique(event_times)),
    varprof = varprof,
    model_type = "cr_gam",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
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
#' @param failcode integer, the code for the event of interest for CIF prediction.
#'   If NULL (default), uses the failcode from the model training.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @export
Predict_CRModel_GAM <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

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

  # Handle failcode parameter
  if (is.null(failcode)) {
    failcode <- modelout$failcode  # Use the failcode from training
  } else {
    if (!is.numeric(failcode) || length(failcode) != 1 || failcode < 1) {
      stop("'failcode' must be a positive integer")
    }
    if (!failcode %in% modelout$all_event_types) {
      stop("failcode ", failcode, " was not present in training data. Available event types: ",
           paste(modelout$all_event_types, collapse = ", "))
    }
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
  # Make Predictions using Aalen-Johansen Estimator
  # ============================================================================
  # Get survival predictions from ALL cause-specific GAM models
  cause_specific_survs <- vector("list", length(modelout$all_event_types))
  names(cause_specific_survs) <- as.character(modelout$all_event_types)

  # Predict survival for each cause
  for (cause in modelout$all_event_types) {
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

    # Extract survival probabilities
    surv_probs_cause <- sf$surv
    surv_times_cause <- sf$time

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

  # Get the time grid from the first model
  surv_times <- sf$time

  # Use Aalen-Johansen estimator to calculate proper CIF
  cif_matrix <- aalenJohansenCIF(
    cause_specific_survs = cause_specific_survs,
    times = surv_times,
    event_of_interest = failcode
  )

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(newtimes)) {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times, observations]
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