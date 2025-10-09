#' @title CRModel_FineGray
#'
#' @description Fit a Fine-Gray model for competing risks using penalized regression
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, ...)
#' @param event_codes character vector identifying the event code(s) to model.
#'   Fine-Gray supports exactly one competing event. If NULL (default), the
#'   first non-zero event code observed in the data is used.
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{fg_model}{the fitted Fine-Gray model object from fastcmprsk::fastCrrp}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "fine_gray"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{train_data}{processed training data used for fitting}
#'   \item{scaling}{list with meanTrain and sdTrain for variable scaling}
#'   \item{loadings}{SVD loadings for dimensionality reduction}
#'
#' @importFrom fastcmprsk fastCrrp Crisk
#' @importFrom stats AIC model.matrix quantile sd
#' @importFrom utils head tail
#' @export
CRModel_FineGray <- function(data, expvars, timevar, eventvar, event_codes = NULL,
                            ntimes = 50, verbose = FALSE) {

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
    event_codes <- available_events[1]
  }

  event_codes <- as.character(event_codes)

  if (length(event_codes) != 1) {
    stop("Fine-Gray model supports exactly one event code. Received ", length(event_codes), ".")
  }

  if (!event_codes %in% available_events) {
    stop("Requested event code ", event_codes, " not present in training data. Available codes: ",
         paste(available_events, collapse = ", "))
  }

  failcode <- suppressWarnings(as.numeric(event_codes))
  if (is.na(failcode)) {
    stop("Fine-Gray requires numeric event codes. Unable to coerce '", event_codes, "' to numeric.")
  }

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == failcode]
  if (length(event_times) == 0) {
    stop("No events of type ", event_codes, " in training data. Cannot fit competing risks model.")
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

  # SVD for dimensionality reduction (keep top 20 components or fewer)
  svdcovmat <- svd(covmat_scaled)
  n_components <- min(c(20, ncol(covmat_scaled)))
  Feat <- (covmat_scaled %*% svdcovmat$v)[, 1:n_components]

  # Fit model with default elastic net parameters (lambda=0.01, alpha=0.5)
  # For simplicity, we'll use fixed parameters instead of grid search
  fg_model <- tryCatch(
    fastcmprsk::fastCrrp(
      fastcmprsk::Crisk(XYTrain[[timevar]], XYTrain[[eventvar]]) ~ Feat,
      lambda = 0.01,
      alpha = 0.5,
      penalty = "ENET",
      standardize = FALSE,
      max.iter = 5000
    ),
    error = function(e) {
      stop("Failed to fit Fine-Gray model: ", e$message)
    }
  )

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


#' @title Predict_CRModel_FineGray
#'
#' @description Get predictions from a fitted Fine-Gray competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_FineGray'
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), generates 50 equally-spaced points from 0 to max observed time.
#'   Can be any positive values - interpolation handles all time points.
#' @param event_of_interest character or numeric scalar indicating the event code
#'   for which CIFs should be returned. If NULL (default), uses the event code
#'   stored in the fitted model. Fine-Gray models can only predict the event they
#'   were trained on.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated
#'     (always includes time 0)}
#'
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit model
#' model <- CRModel_FineGray(data, expvars = c("x1", "x2"),
#'                          timevar = "time", eventvar = "event")
#'
#' # Predict on test data
#' preds <- Predict_CRModel_FineGray(model, test_data)
#'
#' # Predict at specific times
#' preds_custom <- Predict_CRModel_FineGray(model, test_data,
#'                                         newtimes = c(30, 60, 90, 180, 365))
#' }
Predict_CRModel_FineGray <- function(modelout, newdata, newtimes = NULL, event_of_interest = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_finegray")) {
    stop("'modelout' must be output from CRModel_FineGray")
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
    stop("Fine-Gray models can only return CIFs for a single event of interest")
  }

  if (!identical(event_of_interest, modelout$event_codes)) {
    stop("Fine-Gray models can only predict for the event they were trained on (event code = ",
         modelout$event_codes, "). Requested event code: ", event_of_interest)
  }

  failcode <- suppressWarnings(as.numeric(event_of_interest))
  if (is.na(failcode)) {
    stop("Fine-Gray requires numeric event codes. Unable to coerce '", event_of_interest, "' to numeric.")
  }

  # Generate default times if not specified
  if (!is.null(newtimes)) {
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
  # Create model matrix and apply same preprocessing as training
  covmat <- stats::model.matrix(~ -1 + ., data = newdata_prepared)

  # Apply scaling
  covmat_scaled <- scale(covmat,
                         center = modelout$scaling$meanTrain,
                         scale = modelout$scaling$sdTrain)

  # Apply SVD transformation
  n_components <- min(c(20, ncol(covmat_scaled)))
  Feat <- (covmat_scaled %*% modelout$loadings)[, 1:n_components]

  # Get baseline cumulative hazard from the model
  # The model provides Breslow jumps at specific times
  baseline_times <- modelout$fg_model$breslowJump[, 1]
  baseline_haz <- modelout$fg_model$breslowJump[, 2]

  # Compute CIF for each observation
  n_obs <- nrow(Feat)
  cif_matrix <- matrix(NA, nrow = length(baseline_times) + 1, ncol = n_obs)

  # Time 0 has CIF = 0
  cif_matrix[1, ] <- 0

  # Compute CIF at each baseline time point
  for (i in 1:n_obs) {
    # Linear predictor
    lp <- sum(Feat[i, ] * unlist(modelout$fg_model$coef))

    # Cumulative hazard for this observation
    cumhaz <- exp(lp) * baseline_haz

    # Cumulative incidence function (subdistribution hazard)
    cif_vals <- 1 - exp(-cumsum(cumhaz))

    cif_matrix[-1, i] <- cif_vals
  }

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(newtimes)) {
    # Return predictions in native time grid: [times, observations]
    result_cifs <- cif_matrix  # cif_matrix is already [times+1, observations]
    result_times <- c(0, baseline_times)
  } else {
    # Interpolate to new time points
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))

    # Use the standard CIF interpolation utility function
    pred_cifs <- cifMatInterpolaltor(
      probsMat = t(cif_matrix),  # cifMatInterpolaltor expects [observations, times]
      times = c(0, baseline_times),
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
