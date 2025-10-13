#' @title SurvModel_BART
#'
#' @description Fit a BART model for survival outcomes using BART::surv.bart.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param K parameter 'K' for BART (number of time points for discrete approximation)
#' @param ntree number of trees in BART model
#'
#' @return a list containing the following objects:
#' model: fitted BART model object from surv.bart,
#' times: time points used in the BART model,
#' varprof: profile of explanatory variables,
#' expvars: character vector of explanatory variables used,
#' factor_levels: list containing factor levels for consistent prediction.
#'
#' @importFrom BART surv.bart
#' @importFrom stats model.matrix
#' @export
SurvModel_BART <- function(data,
                           expvars,
                           timevar,
                           eventvar,
                           K = 8,
                           ntree = 50) {

  if (missing(data)) stop("argument \"data\" is missing")
  if (missing(expvars)) stop("argument \"expvars\" is missing")
  if (missing(timevar)) stop("argument \"timevar\" is missing")
  if (missing(eventvar)) stop("argument \"eventvar\" is missing")

  # Variable profile
  varprof <- VariableProfile(data, expvars)

  # Store factor levels for prediction
  factor_levels <- lapply(data[, expvars, drop=FALSE], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })

  # Prepare data
  timevarvec <- data[[timevar]]
  eventvarvec <- as.integer(data[[eventvar]] == 1) # Ensure 0/1

  # Create model matrix
  x.train <- as.matrix(stats::model.matrix(~ -1 + ., data = data[, expvars, drop = FALSE]))

  # Fit BART model
  post <- surv.bart(
    x.train = x.train,
    times = timevarvec,
    delta = eventvarvec,
    x.test = x.train,
    K = K,
    ntree = ntree,
    ndpost = 200,
    nskip = 50,
    keepevery = 2L
  )

  # Get times from the fitted model
  times <- post$times

  # Return standardized output
  result <- list(
    model = post,
    times = times,
    varprof = varprof,
    expvars = expvars,
    factor_levels = factor_levels,
    x_train = x.train,
    times_train = timevarvec,
    delta_train = eventvarvec
  )
  class(result) <- c("ml4t2e_surv_bart", "list")
  return(result)
}


#' @title Predict_SurvModel_BART
#'
#' @description Make predictions using a fitted BART survival model.
#'
#' @param modelout the output from 'SurvModel_BART'
#' @param newdata the data for which the predictions are to be calculated
#' @param newtimes optional vector of new time points for interpolation. If NULL, uses model's native time points.
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @importFrom BART surv.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @export
Predict_SurvModel_BART <- function(modelout, newdata, newtimes = NULL) {

  if (missing(modelout)) stop("argument \"modelout\" is missing")
  if (missing(newdata)) stop("argument \"newdata\" is missing")

  # Prepare newdata: ensure factors have same levels as training data
  data_test <- newdata[, modelout$expvars, drop=FALSE]
  for (vari in modelout$expvars){
    if (!is.null(modelout$factor_levels[[vari]])) { # Check if var was factor in training
        train_levels <- modelout$factor_levels[[vari]]
        # Ensure the column exists in newdata before attempting to modify
        if (vari %in% colnames(data_test)) {
            # Convert to character first to handle potential new levels, then factor
            data_test[[vari]] <- factor(as.character(data_test[[vari]]), levels = train_levels)
        }
    }
  }

  # Create test model matrix
  x.test_orig <- stats::model.matrix(~ -1 + ., data = data_test)

  # Prepare BART prediction structure using stored training data
  pre <- surv.pre.bart(
    times = modelout$times_train,
    delta = modelout$delta_train,
    x.train = bartModelMatrix(modelout$x_train),
    x.test = bartModelMatrix(x.test_orig),
    K = modelout$model$K
  )

  # Predict using the fitted BART model
  pred <- predict(modelout$model, pre$tx.test)

  # Extract mean survival probabilities
  # CRITICAL FIX: pred$surv.test.mean is a VECTOR (not matrix) containing
  # survival probabilities for all patients concatenated: [patient1_times, patient2_times, ...]
  # We need to reshape it correctly with byrow=TRUE to get proper patient grouping
  PredMat <- matrix(pred$surv.test.mean, ncol = modelout$model$K, byrow = TRUE)

  # CRITICAL FIX: Sort times and reorder probabilities accordingly
  # BART times are not guaranteed to be sorted, which breaks monotonicity
  time_order <- order(modelout$times)
  sorted_times <- modelout$times[time_order]
  sorted_PredMat <- PredMat[, time_order, drop = FALSE]

  # Transpose to rows=times, cols=observations and add time 0 with probability 1
  Probs <- t(cbind(1, sorted_PredMat))
  Times <- c(0, sorted_times)

  # If newtimes specified, interpolate to those times
  if (!is.null(newtimes)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, newtimes = newtimes)
    Times <- newtimes
  }

  return(list(Probs = Probs, Times = Times))
}
