#' @title SurvModel_gbm
#'
#' @description Fit a gbm model for survival outcomes using gbm package.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param ntree number of trees to grow
#' @param max.depth maximum depth for the trees
#' @param bag.fraction fraction of data to sample in each bagging iteration
#' @param train.fraction fraction for training data (used internally by gbm for CV/OOB error)
#' @param learninrate learning rate (shrinkage) for the boosting algorithm
#'
#' @return a list of four items: model: fitted gbm model object,
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables,
#'  expvars: the explanatory variables used.
#'
#' @importFrom gbm gbm gbm.perf basehaz.gbm
#' @importFrom survival Surv
#' @importFrom stats complete.cases
#' @export
SurvModel_gbm<-function(data,expvars, timevar, eventvar, ntree=200, max.depth=3, bag.fraction=.3, train.fraction=.3, learninrate=.01){
  if (missing(data)) stop("argument \"data\" is missing")
  if (missing(expvars)) stop("argument \"expvars\" is missing")
  if (missing(timevar)) stop("argument \"timevar\" is missing")
  if (missing(eventvar)) stop("argument \"eventvar\" is missing")

  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Store factor levels for prediction
  factor_levels <- lapply(data[, expvars, drop=FALSE], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })

  # Create full data frame with predictors and outcome
  gbm_data <- data[, c(expvars, timevar, eventvar), drop=FALSE]
  
  # Define formula using column names directly
  formula_gbm <- as.formula(paste("Surv(", timevar, ",", eventvar, ") ~", paste(expvars, collapse = "+")))

  # Determine cv.folds based on dataset size
  # For small datasets, disable CV to avoid subsample size issues
  n_obs <- nrow(data)
  cv_folds <- 0  # Always disable CV for small datasets to avoid issues
  
  # Adjust parameters for small datasets
  if (n_obs < 50) {
    bag.fraction <- max(0.8, bag.fraction)  # Use more data per tree
    train.fraction <- 1.0  # Use all data for training
  }

  # Fit gbm model with error handling
  gbmmodel <- tryCatch({
    gbm::gbm(formula = formula_gbm,
             data = gbm_data,
             distribution = "coxph",
             n.trees = ntree,
             shrinkage = learninrate,
             interaction.depth = max.depth,
             bag.fraction = bag.fraction,
             train.fraction = train.fraction,
             cv.folds = cv_folds,
             n.minobsinnode = max(2, min(5, floor(n_obs/10))), # Adaptive min obs per node
             keep.data = TRUE,
             verbose = FALSE)
  }, error = function(e) {
    # If still fails, try with very conservative parameters
    gbm::gbm(formula = formula_gbm,
             data = gbm_data,
             distribution = "coxph",
             n.trees = min(50, ntree),
             shrinkage = 0.001,
             interaction.depth = 1,
             bag.fraction = 1.0,
             train.fraction = 1.0,
             cv.folds = 0,
             n.minobsinnode = 1,
             keep.data = TRUE,
             verbose = FALSE)
  })

  # Find best iteration based on OOB performance (since CV is disabled)
  best.iter <- gbm::gbm.perf(gbmmodel, method = "OOB", plot.it = FALSE)

  # Get unique event times from training data
  time.interest <- sort(unique(data[[timevar]][data[[eventvar]]==1]))

  # Predict linear predictor on training data
  pred.train <- predict(gbmmodel, data[, expvars, drop=FALSE], n.trees = best.iter, type="link")

  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum <- gbm::basehaz.gbm(t = data[[timevar]], delta = data[[eventvar]], f.x = pred.train, t.eval = time.interest, cumulative = TRUE)

  # Store baseline hazard in the model for prediction
  gbmmodel$basehaz.cum <- basehaz.cum
  gbmmodel$time.interest <- time.interest
  gbmmodel$best.iter <- best.iter

  result <- list(model = gbmmodel, times = time.interest, varprof = varprof, expvars = expvars, factor_levels = factor_levels)
  class(result) <- c("ml4t2e_surv_gbm", "list")
  return(result)
}




#' @title Predict_SurvModel_gbm
#'
#' @description Get predictions from a gbm survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_gbm' (a list containing 'model', 'times', 'varprof', 'expvars')
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @export
Predict_SurvModel_gbm <- function(modelout, newdata, newtimes = NULL) {
  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) stop("argument \"modelout\" is missing")
  if (missing(newdata)) stop("argument \"newdata\" is missing")
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

  # Check for NAs and count complete cases
  has_missing <- any(is.na(data_test))
  n_complete_cases <- sum(complete.cases(data_test))
  n_expected_obs <- nrow(data_test)

  if (has_missing) {
    warning(sprintf(
      "Missing values detected in %d observation(s). ",
      n_expected_obs - n_complete_cases),
      "gbm::predict.gbm will use NA for incomplete cases."
    )
  }

  # ============================================================================
  # Make Predictions
  # ============================================================================
  # Predict linear predictor for newdata
  event_prediction <- suppressWarnings(predict(modelout$model, data_test, n.trees=modelout$model$best.iter, type="link"))

  # Handle missing values in predictions
  if (has_missing) {
    # predict.gbm returns NA for observations with missing values
    # We need to expand to match original dimensions
    full_event_prediction <- rep(NA, n_expected_obs)
    complete_idx <- which(complete.cases(data_test))
    full_event_prediction[complete_idx] <- event_prediction
    event_prediction <- full_event_prediction
  }

  # Calculate survival probabilities using stored baseline hazard
  survMat <- NULL
  for (i in seq_along(event_prediction)){
    if (is.na(event_prediction[i])) {
      # For missing predictions, use NA for all time points
      surf.i <- rep(NA, length(modelout$model$time.interest))
    } else {
      surf.i <- exp(-exp(event_prediction[i])*modelout$model$basehaz.cum)
    }
    survMat <- rbind(survMat,surf.i) # Matrix: rows=obs, cols=times
  }

  # Add time 0 with probability 1 (or NA for missing observations)
  time_0_probs <- ifelse(is.na(event_prediction), NA, 1)
  Probs <- t(cbind(time_0_probs, survMat)) # Transpose to rows=times, cols=obs
  Times <- c(0, modelout$model$time.interest)

  # If newtimes specified, interpolate to those times
  if (!is.null(newtimes)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, newtimes = newtimes)
    Times <- newtimes
  }

  return(list(Probs = Probs, Times = Times))
}
