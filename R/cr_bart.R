#' @title CRModel_BART
#'
#' @description Fit a BART model for competing risks outcomes using BART::crisk.bart.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2 where 0=censored, 1=event of interest, 2=competing event)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param K parameter 'K' for BART (number of time points for discrete approximation, default: 10)
#' @param ntree number of trees in BART model (default: 100)
#' @param ndpost number of posterior draws to save (default: 1000)
#' @param nskip number of MCMC iterations to discard as burn-in (default: 100)
#' @param keepevery interval at which to keep posterior draws (default: 10)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{bart_model}{the fitted BART model object from BART::crisk.bart}
#'   \item{times}{vector of unique event times in the training data}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_bart"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{x_train}{design matrix for the training data}
#'   \item{times_train}{time variable values from training data}
#'   \item{delta_train}{event variable values from training data}
#'
#' @importFrom BART crisk.bart
#' @importFrom stats model.matrix
#' @export
CRModel_BART <- function(data, expvars, timevar, eventvar, failcode = 1,
                        K = 10, ntree = 100, ndpost = 1000, nskip = 100,
                        keepevery = 10, verbose = FALSE) {

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
  # Fit BART model with retry logic for robustness
  failcount <- 0
  bart_model <- NULL

  while (is.null(bart_model) && failcount < 4) {
    failcount <- failcount + 1
    bart_model <- tryCatch(
      BART::crisk.bart(
        x.train = x_train,
        times = times_train,
        delta = delta_train,
        x.test = x_train,  # Include training data for predictions
        K = K,
        ntree = ntree,
        ndpost = ndpost,
        nskip = nskip,
        keepevery = keepevery,
        numcut = 2  # Small number of cuts for speed
      ),
      error = function(e) {
        if (verbose) message("BART fitting attempt ", failcount, " failed: ", e$message)
        NULL
      }
    )
  }

  if (is.null(bart_model)) {
    stop("Failed to fit BART model after ", failcount, " attempts")
  }

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] != 0, timevar]))

  # Get time range
  time_range <- range(data[data[[eventvar]] != 0, timevar])

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
    failcode = failcode,
    time_range = time_range,
    x_train = x_train,
    times_train = times_train,
    delta_train = delta_train
  )

  class(result) <- "ml4t2e_cr_bart"
  return(result)
}


#' @title Predict_CRModel_BART
#'
#' @description Get predictions from a CR BART model for a test dataset.
#'
#' @param modelout the output from 'CRModel_BART' (a list containing model and metadata)
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the model's native time points.
#'   Can be any positive values - interpolation handles all time points.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated
#'     (always includes time 0)}
#'
#' @importFrom BART crisk.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_BART <- function(modelout, newdata, newtimes = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_bart")) {
    stop("'modelout' must be output from CRModel_BART")
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
  # Create test model matrix
  x_test <- as.matrix(stats::model.matrix(~ -1 + ., data = newdata_prepared))

  # Prepare BART prediction structure
  # Note: BART crisk.pre.bart expects duplicated test data for some reason
  # (this seems to be a BART package requirement)
  test_duplicated <- rbind(x_test, x_test)

  pre <- BART::crisk.pre.bart(
    time = modelout$times_train,
    delta = modelout$delta_train,
    x.train = BART::bartModelMatrix(modelout$x_train),
    x.test = BART::bartModelMatrix(test_duplicated),
    x.train2 = BART::bartModelMatrix(modelout$x_train),
    x.test2 = BART::bartModelMatrix(test_duplicated),
    K = modelout$bart_model$K
  )

  # Predict using the fitted BART model
  pred <- predict(modelout$bart_model, pre$tx.test, pre$tx.test2)

  # Extract CIF predictions
  # pred$cif.test.mean is organized as [obs1_time1, obs1_time2, ..., obs1_timeK, obs2_time1, ...]
  # We need to reshape it to [observations, times]
  N <- nrow(x_test)  # Number of actual test observations
  K <- modelout$bart_model$K  # Number of time points
  
  # cif.test.mean has length N*K (but we duplicated test data, so it's 2N*K total)
  # Take only the first N observations (not the duplicated ones)
  cif_vector <- pred$cif.test.mean[1:(N*K)]
  
  # Reshape to matrix [observations, times]
  cif_matrix <- matrix(cif_vector, nrow = N, ncol = K, byrow = TRUE)

  # Add time 0 with CIF = 0
  cif_with_t0 <- cbind(0, cif_matrix)
  times_with_t0 <- c(0, modelout$bart_model$times)

  # ============================================================================
  # Apply Interpolation if needed
  # ============================================================================
  if (is.null(newtimes)) {
    # Use model's native time points
    result_cifs <- cif_with_t0  # Keep as [observations, times]
    result_times <- times_with_t0
  } else {
    # Interpolate to new time points
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))

    # Use cifMatInterpolaltor for interpolation
    result_cifs <- cifMatInterpolaltor(
      probsMat = cif_with_t0,  # [observations, times]
      times = times_with_t0,
      newtimes = newtimes
    )
    # cifMatInterpolaltor returns [newtimes, observations], transpose to [observations, newtimes]
    result_cifs <- t(result_cifs)
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
