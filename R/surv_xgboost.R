#' @title xgb.train.surv (Internal Helper)
#' @description Internal function to train xgboost for survival and estimate baseline hazard.
#' Adapts xgboost training for Cox PH objective and calculates an optimized baseline hazard.
#' @param params list of xgboost parameters (must include objective='survival:cox', eval_metric='cox-nloglik').
#' @param data training data matrix (features only).
#' @param label training labels (negative time for censored, positive time for event).
#' @param weight optional observation weights.
#' @param nrounds number of boosting rounds.
#' @param watchlist list for monitoring evaluation metrics during training.
#' @param verbose verbosity level.
#' @param print_every_n print evaluation metric every n rounds.
#' @param early_stopping_rounds rounds to wait for improvement before stopping.
#' @param save_period frequency to save model during training (0 = never).
#' @param save_name filename for saved models.
#' @param xgb_model existing xgboost model to continue training from.
#' @param callbacks list of callback functions for training.
#' @param ... other arguments passed to xgboost::xgb.train.
#' @return An xgb.Booster object with an added 'baseline_hazard' component.
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom survival coxph Surv basehaz
#' @importFrom stats optim
#' @importFrom pec pec crps
#' @noRd
xgb.train.surv <- function(params = list(), data, label, weight = NULL, nrounds,
                           watchlist = list(), verbose = 1, print_every_n = 1L,
                           early_stopping_rounds = NULL, save_period = NULL,
                           save_name = "xgboost_surv.model", xgb_model = NULL, callbacks = list(), ...) {

  # Ensure correct objective and eval_metric
  if (length(params) > 0) {
    if (!is.null(params$objective) && params$objective != "survival:cox") stop("params objective must be set to survival:cox")
    if (!is.null(params$eval_metric) && params$eval_metric != "cox-nloglik") stop("params eval_metric must be set to cox-nloglik")
  }
  # Set defaults if not provided
  if (is.null(params$objective)) params$objective <- "survival:cox"
  if (is.null(params$eval_metric)) params$eval_metric <- "cox-nloglik"


  if(is.null(weight)) weight <- rep(1, nrow(data))

  # Create DMatrix
  data_DMatrix <- xgboost::xgb.DMatrix(data = data, label = label, weight = weight)

  # Train xgboost model
  xgboost_model <- xgboost::xgb.train(
    params = params, data = data_DMatrix, nrounds = nrounds, watchlist = watchlist, verbose = verbose,
    print_every_n = print_every_n, early_stopping_rounds = early_stopping_rounds, save_period = save_period,
    save_name = save_name, xgb_model = xgb_model, callbacks = callbacks, ...
  )


  # --- Proper Baseline Hazard Estimation (Breslow-like approach) ---
  # Create data frame with time, status, and linear predictors
  # XGBoost label convention: positive time = event, negative time = censored
  data_data.frame <- data.frame(data, time = abs(label), status = ifelse(label > 0, 1, 0))
  
  # Get linear predictors for all training observations
  lp_train <- xgboost:::predict.xgb.Booster(xgboost_model, data)  
  data_data.frame$lp <- lp_train
  data_data.frame$exp_lp <- exp(lp_train)

  # Get unique event times and estimate baseline hazard using Breslow method
  unique_event_times <- sort(unique(data_data.frame$time[data_data.frame$status == 1]))
  
  if (length(unique_event_times) == 0) {
    # No events - create minimal baseline hazard
    baseline_hazard <- data.frame(time = c(0, max(data_data.frame$time, na.rm = TRUE)), 
                                 hazard = c(0, 0.01))
  } else {
    # Breslow baseline hazard estimation
    baseline_hazard <- data.frame(time = numeric(0), hazard = numeric(0))
    
    for(t in unique_event_times) {
      # Number of events at time t
      d_t <- sum(data_data.frame$time == t & data_data.frame$status == 1)
      
      # Risk set at time t (all individuals still at risk)
      risk_set <- data_data.frame$time >= t
      
      # Sum of exp(lp) for individuals in risk set
      sum_exp_lp <- sum(data_data.frame$exp_lp[risk_set])
      
      # Breslow baseline hazard increment: d_t / sum(exp(lp))
      h0_t <- if(sum_exp_lp > 0) d_t / sum_exp_lp else 0.01
      
      baseline_hazard <- rbind(baseline_hazard, data.frame(time = t, hazard = h0_t))
    }
    
    # Ensure baseline hazard starts at time 0 with hazard 0
    if(nrow(baseline_hazard) > 0 && baseline_hazard[1, "time"] != 0) {
      baseline_hazard <- rbind(data.frame(time = 0, hazard = 0), baseline_hazard)
    }
  }

  # Store the optimized baseline hazard in the model object
  xgboost_model$baseline_hazard <- baseline_hazard
  class(xgboost_model) <- c("xgb.Booster.surv", class(xgboost_model)) # Add specific class
  return(xgboost_model)
}

#' @title predict.xgb.Booster.surv (Internal Helper)
#' @description Prediction method for survival xgboost models trained with xgb.train.surv.
#' @param object A model object of class 'xgb.Booster.surv'.
#' @param newdata New data matrix for prediction.
#' @param type Type of prediction: "risk" (linear predictor/log hazard ratio) or "surv" (survival probability).
#' @param times Optional numeric vector of times at which to predict survival probabilities. If NULL, uses times from baseline hazard.
#' @return Predicted risk scores or survival probabilities matrix (rows=observations, cols=times).
#' @importFrom stats approxfun
#' @noRd
predict.xgb.Booster.surv <- function(object, newdata, type = "risk", times = NULL) {
  # Predict the linear predictor (log hazard ratio) using the base xgboost model
  lp <- xgboost:::predict.xgb.Booster(object, newdata) # Ensure calling the base predict method

  if (type == "risk") {
    return(lp) # Return linear predictor
  } else if (type == "surv") {
    baseline_hazard <- object$baseline_hazard
    if (is.null(baseline_hazard)) {
        stop("Baseline hazard not found in the model object. Cannot predict survival probabilities.")
    }

    # Determine evaluation times
    if (!is.null(times)) {
      eval_times <- sort(unique(c(0, times))) # Ensure 0 is included and times are sorted
    } else {
      eval_times <- baseline_hazard[, "time"] # Use times from baseline hazard
    }

    # Calculate baseline cumulative hazard H0(t) at evaluation times
    # Sum baseline hazard increments up to each time point
    H0_t <- numeric(length(eval_times))
    for(i in seq_along(eval_times)) {
      # Sum baseline hazard increments up to time eval_times[i]
      relevant_times <- baseline_hazard$time[baseline_hazard$time <= eval_times[i]]
      if(length(relevant_times) > 0) {
        H0_t[i] <- sum(baseline_hazard$hazard[baseline_hazard$time %in% relevant_times])
      } else {
        H0_t[i] <- 0
      }
    }

    # Calculate survival probability: S(t|X) = exp(-H0(t) * exp(LP))
    # LP is the linear predictor from xgboost
    expLP <- exp(lp)
    surv <- exp(-outer(expLP, H0_t, "*")) # outer product: rows=observations, cols=times

    colnames(surv) <- eval_times
    return(surv)
  } else {
    stop('type must be one of "risk", "surv"')
  }
}





#' @title SurvModel_xgboost
#'
#' @description Fit an xgboost model for survival outcomes.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param eta learning rate (default: 0.01)
#' @param max_depth maximum depth of trees (default: 5)
#' @param nrounds number of boosting rounds (default: 100)
#' @param ... additional parameters passed to xgboost
#'
#' @return a list of three items: model: fitted xgboost model object (class xgb.Booster.surv),
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables.
#'
#' @importFrom xgboost xgb.DMatrix
#' @importFrom stats model.matrix
#' @export
SurvModel_xgboost<-function(data, expvars, timevar, eventvar, eta = 0.01, max_depth = 5, nrounds = 100, ...){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Prepare data matrix (handle factors)
  X<-stats::model.matrix(~-1+., data[,c(expvars), drop=FALSE])
  # Clean column names for xgboost if necessary
  # colnames(X)<-make.names(colnames(X), unique = TRUE) # Optional: if names are problematic

  # Prepare labels: negative time for censored, positive time for event
  ytimes = data[[timevar]]
  yevents<-as.numeric(data[[eventvar]] == 1) # Ensure 0/1
  yTrain<-ytimes
  yTrain[yevents==0]<-(-yTrain[yevents==0])

  # Define xgboost parameters with user-provided values
  params <- list(objective='survival:cox',
                 eval_metric = "cox-nloglik",
                 eta = eta,
                 max_depth = max_depth,
                 ...)

  # Train the survival xgboost model using the internal helper function
  bst <- xgb.train.surv(data=X, label=yTrain, params = params, nrounds = nrounds, verbose = 0)

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] == 1, timevar]))

  # Store feature names in the model for prediction consistency
  bst$feature_names <- colnames(X)

  return(list(model = bst, times = times, varprof = varprof, expvars = expvars))
}







#' @title Predict_SurvModel_xgboost
#'
#' @description Get predictions from an xgboost survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_xgboost'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @importFrom stats model.matrix
#' @importFrom xgboost xgb.DMatrix
#' @export
Predict_SurvModel_xgboost<-function(modelout, newdata){
  # Get the feature names from the stored model
  feature_names <- modelout$model$feature_names

  # Apply the same model.matrix transformation to newdata as was done during training
  # This ensures factors are properly converted to dummy variables
  X_new <- stats::model.matrix(~-1+., newdata[, modelout$expvars, drop=FALSE])

  # Ensure X_new has the same columns as training (in case of missing factor levels)
  missing_cols <- setdiff(feature_names, colnames(X_new))
  for(col in missing_cols){
    X_new[[col]] <- 0 # Add missing columns with 0
  }
  # Ensure correct column order
  X_new <- X_new[, feature_names, drop = FALSE]

  # Create DMatrix for prediction
  dtest <- xgboost::xgb.DMatrix(data=X_new)

  # Get the actual model object - handle both 'model' and 'learner' fields
  actual_model <- if (!is.null(modelout$model)) modelout$model else modelout$learner
  
  if (is.null(actual_model)) {
    stop("No model found in modelout$model or modelout$learner")
  }
  
  # Get evaluation times from the baseline hazard stored in the model
  times <- sort(unique(c(0, actual_model$baseline_hazard$time)))

  # Predict survival probabilities using the specific predict method
  preds <- predict.xgb.Booster.surv(actual_model, dtest, type = "surv", times=times)

  # Transpose predictions to match standard output: rows=times, cols=observations
  return(list(Probs=t(preds), Times=times))
}
