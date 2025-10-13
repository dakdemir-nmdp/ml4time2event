#' @title SurvModel_glmnet
#'
#' @description Fit a penalized Cox model for survival outcomes using glmnet.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param alpha elastic net mixing parameter (1=lasso, 0=ridge)
#' @param maxit maximum number of iterations for glmnet
#' @param nfolds number of folds for cross-validation in cv.glmnet
#'
#' @return a list containing the following objects:
#' model: cv.glmnet fit object,
#' times: unique event times from the training data,
#' varprof: profile of explanatory variables,
#' expvars: character vector of explanatory variables used,
#' factor_levels: list containing factor levels for consistent prediction.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom survival Surv survfit
#' @importFrom stats model.matrix
#' @export
SurvModel_glmnet <- function(data,
                             expvars,
                             timevar,
                             eventvar,
                             alpha = 0.5,
                             maxit = 5000,
                             nfolds = 30) {

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

  # Create model matrix
  TrainMat <- stats::model.matrix(~ . - 1, data = data[, expvars, drop = FALSE])

  # Create survival object
  yTrain <- survival::Surv(data[[timevar]], data[[eventvar]] == 1)

  # Store sample of training data for prediction (needed for survfit)
  sample_size <- min(500, nrow(data))
  sample_idx <- sample(seq_len(nrow(data)), sample_size)

  # Create model matrix for the sample
  train_sample <- data[sample_idx, expvars, drop = FALSE]
  TrainMat_sample <- stats::model.matrix(~ . - 1, data = train_sample)

  # Create survival object for the sample
  y_sample <- survival::Surv(data[[timevar]][sample_idx], data[[eventvar]][sample_idx] == 1)

  # Fit penalized Cox model
  cv.fit <- glmnet::cv.glmnet(
    x = TrainMat,
    y = yTrain,
    family = "cox",
    alpha = alpha,
    maxit = maxit,
    nfolds = nfolds
  )

  # Get unique event times
  sfitTrain <- survival::survfit(yTrain ~ 1)
  times <- sfitTrain$time

  # Return standardized output
  result <- list(
    model = cv.fit,
    times = times,
    varprof = varprof,
    expvars = expvars,
    factor_levels = factor_levels,
    train_sample = train_sample,
    train_matrix = TrainMat_sample,
    y_train = y_sample
  )
  class(result) <- c("ml4t2e_surv_glmnet", "list")
  return(result)
}

#' @title Predict_SurvModel_glmnet
#'
#' @description Get predictions from a glmnet survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_glmnet'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @importFrom survival survfit
#' @importFrom stats model.matrix
#' @export
Predict_SurvModel_glmnet <- function(modelout, newdata, newtimes = NULL) {

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

  # Create combined data to ensure factor levels match exactly
  combined_data <- rbind(modelout$train_sample, data_test)
  combined_mat <- stats::model.matrix(~ . - 1, data = combined_data)

  # Extract test portion
  n_train <- nrow(modelout$train_sample)
  test_mat <- combined_mat[-seq_len(n_train), , drop = FALSE]

  # Get training matrix column names for consistency
  train_cols <- colnames(modelout$train_matrix)
  
  # Ensure test matrix has the same columns as training matrix
  missing_cols <- setdiff(train_cols, colnames(test_mat))
  for (col in missing_cols) {
    test_mat <- cbind(test_mat, 0)
    colnames(test_mat)[ncol(test_mat)] <- col
  }
  
  # Reorder columns to match training matrix exactly
  test_mat <- test_mat[, train_cols, drop = FALSE]

  # Get predictions using survfit
  sfit <- survival::survfit(
    modelout$model,
    s = "lambda.min",
    x = modelout$train_matrix,
    y = modelout$y_train,
    newx = test_mat
  )

  # Extract survival probabilities and times
  surv_probs <- sfit$surv
  times <- sfit$time

  # Ensure time 0 with probability 1 is included
  if (sum(times == 0) == 0) {
    times <- c(0, times)
    surv_probs <- rbind(rep(1, ncol(surv_probs)), surv_probs)
  }

  Probs <- surv_probs
  Times <- times

  # If newtimes specified, interpolate to those times
  if (!is.null(newtimes)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, newtimes = newtimes)
    Times <- newtimes
  }

  return(list(
    Probs = Probs,
    Times = Times
  ))
}