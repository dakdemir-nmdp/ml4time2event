#'
#'
#'
#'
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
  # For small datasets, disable CV and use all trees
  n_obs <- nrow(data)
  cv_folds <- 0  # Disable CV for all cases in this simplified version

  # Calculate safe parameters to satisfy GBM constraint:
  # nTrain * bag.fraction > 2 * n.minobsinnode + 1
  # Start with conservative defaults
  safe_minobsinnode <- 1
  safe_bag_fraction <- 1.0
  safe_train_fraction <- 1.0

  # For larger datasets, we can use more aggressive parameters
  if (n_obs >= 100) {
    # Calculate max safe n.minobsinnode given bag.fraction
    # We want: n_obs * bag.fraction > 2 * n.minobsinnode + 1
    # Solve for n.minobsinnode: n.minobsinnode < (n_obs * bag.fraction - 1) / 2
    safe_bag_fraction <- min(bag.fraction, 0.8)
    max_safe_minobs <- floor((n_obs * safe_bag_fraction - 1) / 2)
    safe_minobsinnode <- max(1, min(5, max_safe_minobs))
  }

  # Fit gbm model with error handling
  message("Fitting Survival GBM model...")
  gbmmodel <- tryCatch({
    suppressMessages(gbm::gbm(formula = formula_gbm,
             data = gbm_data,
             distribution = "coxph",
             n.trees = ntree,
             shrinkage = learninrate,
             interaction.depth = max.depth,
             bag.fraction = safe_bag_fraction,
             train.fraction = safe_train_fraction,
             cv.folds = cv_folds,
             n.minobsinnode = safe_minobsinnode,
             keep.data = TRUE,
             verbose = FALSE))
  }, error = function(e) {
    # If still fails, try with minimal parameters
    warning("GBM failed with error: ", e$message, ". Trying minimal parameters...")
    suppressMessages(gbm::gbm(formula = formula_gbm,
             data = gbm_data,
             distribution = "coxph",
             n.trees = min(50, ntree),
             shrinkage = 0.01,
             interaction.depth = 1,
             bag.fraction = 1.0,
             train.fraction = 1.0,
             cv.folds = 0,
             n.minobsinnode = 1,
             keep.data = TRUE,
             verbose = FALSE))
  })

  # For small datasets, just use all trees to avoid OOB/CV issues
  if (n_obs < 200) {  # Use all trees for datasets smaller than 200
    best.iter <- gbmmodel$n.trees
  } else {
    # Find best iteration using available method
    best.iter <- tryCatch({
      if (gbmmodel$cv.folds > 0) {
        gbm::gbm.perf(gbmmodel, method = "cv", plot.it = FALSE)
      } else {
        gbm::gbm.perf(gbmmodel, method = "OOB", plot.it = FALSE)
      }
    }, error = function(e) {
      # Fallback: use all trees
      gbmmodel$n.trees
    })
  }

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




#'
#'
#'
#'
Predict_SurvModel_gbm <- function(modelout, newdata, new_times = NULL) {
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

  # If new_times specified, interpolate to those times
  if (!is.null(new_times)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, new_times = new_times)
    Times <- new_times
  }

  return(list(Probs = Probs, Times = Times))
}
