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
  # Model Fitting - Cause-Specific GAM Model
  # ============================================================================
  if (verbose) cat("Fitting cause-specific GAM model for event type", failcode, "...\n")

  # Create cause-specific data: event = 1 if failcode, 0 if censored, NA if other competing event
  XYTrain$status_cs <- ifelse(XYTrain[[eventvar]] == 0, 0,  # censored
                              ifelse(XYTrain[[eventvar]] == failcode, 1, NA))  # event of interest or competing

  # Remove competing events (treat as censored for this cause-specific model)
  XYTrain_cs <- XYTrain[!is.na(XYTrain$status_cs), , drop = FALSE]

  if (nrow(XYTrain_cs) < 10) {
    stop("Insufficient data after removing competing events. Need at least 10 observations.")
  }

  # Identify variable types for GAM formula construction
  numvars <- expvars[which(sapply(as.data.frame(XYTrain_cs[, expvars, drop=FALSE]), is.numeric))]
  fctvars <- expvars[which(sapply(as.data.frame(XYTrain_cs[, expvars, drop=FALSE]), function(x) {
    (is.factor(x) | is.character(x))
  }))]

  # Convert characters to factors
  for (i in fctvars) {
    if (is.character(XYTrain_cs[[i]])) {
      XYTrain_cs[[i]] <- as.factor(XYTrain_cs[[i]])
    }
  }

  # Separate variables based on number of levels for smoothing/shrinkage
  catvarstoshrink <- if (length(fctvars) > 0) {
    fctvars[which(sapply(as.data.frame(XYTrain_cs[, fctvars, drop=FALSE]), function(x) {
      length(levels(x)) > shrinkTreshold
    }))]
  } else {
    c()
  }
  catvarsnottoshrink <- if (length(fctvars) > 0) {
    fctvars[which(sapply(as.data.frame(XYTrain_cs[, fctvars, drop=FALSE]), function(x) {
      length(levels(x)) <= shrinkTreshold
    }))]
  } else {
    c()
  }
  numvarstosmooth <- if (length(numvars) > 0) {
    numvars[which(sapply(as.data.frame(XYTrain_cs[, numvars, drop=FALSE]), function(x) {
      length(unique(x[!is.na(x)])) > shrinkTreshold
    }))]
  } else {
    c()
  }
  numvarsnottosmooth <- if (length(numvars) > 0) {
    numvars[which(sapply(as.data.frame(XYTrain_cs[, numvars, drop=FALSE]), function(x) {
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
  gam_model <- mgcv::gam(
    formGAM,
    family = mgcv::cox.ph(),
    data = XYTrain_cs,
    weights = XYTrain_cs$status_cs, # Use cause-specific event indicator as weights
    select = TRUE # Enable shrinkage via double penalty approach
  )

  # Create baseline model for prediction using score2proba approach
  # Get linear predictors for training data
  train_linear_preds <- stats::predict(gam_model,
                                       newdata = XYTrain_cs,
                                       type = "link")

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

  # Store baseline model in the GAM model object
  gam_model$baseline_model <- baseline_info$model
  gam_model$baseline_sf <- baseline_info$sf

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cause-specific GAM model fitting complete.\n")

  result <- list(
    gam_model = gam_model,
    times = sort(unique(XYTrain_cs[[timevar]][XYTrain_cs$status_cs == 1])),
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
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats predict
#' @export
Predict_CRModel_GAM <- function(modelout, newdata, newtimes = NULL) {

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
  # Get linear predictors from the GAM model
  linear_preds <- stats::predict(modelout$gam_model,
                                 newdata = newdata_prepared,
                                 type = "link", # Get linear predictor
                                 se.fit = TRUE)

  # Use the stored baseline hazard from the training fit
  # and the new scores (linear_preds$fit) to get survival probabilities for the new data
  sf <- survival::survfit(
    modelout$gam_model$baseline_model, # Use the stored Cox model
    newdata = data.frame("score" = linear_preds$fit),
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