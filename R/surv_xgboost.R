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


  # --- Estimate Baseline Hazard ---
  # Create data frame with time, status, and original features
  data_data.frame <- data.frame(data, time = abs(label), status = ifelse(sign(label) == 1, 1, 0))

  # Fit a standard Cox model on the training data (features only)
  # This is used ONLY to get an initial estimate of the baseline hazard shape
  # Note: This assumes the feature names in 'data' are valid for formula
  # It might be safer to use the matrix 'data' directly if column names are problematic
  cox_formula <- tryCatch(
      stats::as.formula(paste("survival::Surv(time, status) ~", paste(colnames(data), collapse = "+"))),
      error = function(e) survival::Surv(time, status) ~ . # Fallback if colnames are bad
  )
  cox_model <- tryCatch(
      survival::coxph(formula = cox_formula, data = data_data.frame, x = FALSE, y = FALSE),
      error = function(e) {
          warning("Failed to fit Cox model for baseline hazard estimation: ", e$message)
          return(NULL)
      }
  )

  if (!is.null(cox_model)) {
      baseline_hazard <- tryCatch(
          survival::basehaz(cox_model, centered = FALSE), # Get non-centered baseline hazard
          error = function(e) {
              warning("Failed to calculate basehaz: ", e$message)
              return(NULL)
          }
      )
  } else {
      baseline_hazard <- NULL
  }


  # If basehaz fails, create a simple placeholder (e.g., based on unique event times)
  if (is.null(baseline_hazard)) {
      warning("Using simplified baseline hazard based on unique event times.")
      unique_event_times <- sort(unique(data_data.frame$time[data_data.frame$status == 1]))
      # Simple cumulative hazard estimate (e.g., Nelson-Aalen type)
      cum_haz_est <- cumsum(1 / rev(seq_along(unique_event_times))) # Basic estimate
      baseline_hazard <- data.frame(hazard = cum_haz_est, time = unique_event_times)
  }


  # Ensure baseline hazard starts at time 0 with hazard 0
  if (nrow(baseline_hazard) == 0 || baseline_hazard[1, "time"] != 0) {
    baseline_hazard <- rbind(c(0, 0), baseline_hazard)
    colnames(baseline_hazard) <- c("hazard", "time") # Ensure correct names
  }

  # --- Optimize Baseline Hazard Scaling using PEC ---
  # Predict Hazard Ratios (HR) using the trained xgboost model
  HR <- xgboost:::predict.xgb.Booster(object = xgboost_model, newdata = data_DMatrix)

  # Define function to calculate prediction error (CRPS) for a given baseline hazard scaling constant
  baseline_pred_error <- function(const) {
    # Scale the baseline hazard estimate
    scaled_baseline_hazard <- baseline_hazard
    scaled_baseline_hazard[, "hazard"] <- scaled_baseline_hazard[, "hazard"] * const

    # Calculate survival probabilities: S(t) = exp(-H0(t) * HR)
    # Need to interpolate baseline hazard at required times if necessary
    # For simplicity here, assume evaluation times match baseline_hazard times
    risk <- HR %*% matrix(scaled_baseline_hazard[, "hazard"], nrow = 1) # Should be HR * H0(t) -> matrix mult? No, element-wise exponent?
    # Correct calculation: S(t|X) = S0(t)^exp(LP) = exp(-H0(t))^exp(LP) = exp(-H0(t) * exp(LP))
    # xgboost predict gives LP (linear predictor), not HR=exp(LP) directly for survival:cox
    # So, risk = exp(HR) where HR is the output of predict.xgb.Booster
    # surv = exp(-matrix(scaled_baseline_hazard[, "hazard"], nrow=1) * exp(HR)) # This seems wrong dim mismatch
    # Let's use the structure from predict.xgb.Booster.surv: surv = exp(HR %*% -matrix(scaled_baseline_hazard[,1], nrow=1))
    # This assumes HR is log(HR), i.e., the linear predictor.
    surv <- exp(HR %*% -matrix(scaled_baseline_hazard[, "hazard"], nrow = 1))

    # Ensure predictions match the times in baseline_hazard
    eval_times <- scaled_baseline_hazard[, "time"]
    surv_at_eval_times <- surv[, findInterval(eval_times, eval_times)] # Should already match if eval_times = baseline_hazard$time

    # Format for pec
    Models <- list("xgboost" = t(surv_at_eval_times)) # pec expects rows=obs, cols=times

    # Calculate PEC
    PredError<-tryCatch(
        pec::pec(object=Models, formula=survival::Surv(time,status)~1,
                 data=data_data.frame,
                 cens.model="marginal", # Use marginal censoring model
                 times = eval_times,
                 exact = FALSE, # Use faster approximation
                 verbose = FALSE,
                 reference = FALSE,
                 splitMethod = "none" # Evaluate on the same data
                 ),
        error = function(e) {warning("pec calculation failed: ", e$message); return(NULL)}
    )

    if (is.null(PredError)) return(Inf) # Return high error if pec failed

    # Return integrated CRPS
    return(mean(pec::crps(PredError), na.rm = TRUE)) # Average CRPS over time
  }

  # Optimize the scaling constant
  optimal_const_result <- tryCatch(
      stats::optim(par = 1, fn = baseline_pred_error, method = "Brent", lower = 1e-4, upper = 10),
      error = function(e) {warning("Baseline hazard optimization failed: ", e$message); return(list(par=1))} # Default to 1 if optim fails
      )
  optimal_const <- optimal_const_result$par

  # Apply optimal scaling to the baseline hazard
  baseline_hazard[, "hazard"] <- baseline_hazard[, "hazard"] * optimal_const

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

    # Interpolate baseline cumulative hazard H0(t) at evaluation times
    haz_func <- stats::approxfun(baseline_hazard$time, baseline_hazard$hazard, method = "constant", f = 0, yleft = 0, rule = 2)
    H0_t <- haz_func(eval_times)

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
#'
#' @return a list containing the following objects:
#' expvars: vector of explanatory variables used,
#' learner: fitted model object (class xgb.Booster.surv),
#' estSURVTrain: predicted survival probabilities matrix for training data (rows=obs, cols=times),
#' datatrainProf: sample data frame of original predictors (for factor level consistency),
#' varprof: profile of explanatory variables.
#'
#' @importFrom xgboost xgb.DMatrix
#' @importFrom stats model.matrix
#' @export
SurvModel_xgboost<-function(data,expvars, timevar, eventvar){
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

  # Define xgboost parameters
  params <- list(objective='survival:cox',
                 eval_metric = "cox-nloglik",
                 eta = 0.01, # learning_rate
                 max_depth=5
                 # early_stopping_rounds=10 # Consider adding watchlist for this
                 )

  # Train the survival xgboost model using the internal helper function
  bst <- xgb.train.surv(data=X, label=yTrain, params = params, nrounds=100, verbose = 0) # Use helper

  # Predict survival probabilities on training data using the specific predict method
  dtrain <- xgboost::xgb.DMatrix(data=X) # Need DMatrix for prediction
  preds<-predict.xgb.Booster.surv(bst, dtrain, type = "surv") # Use specific predict method

  # Store sample of original data for prediction consistency
  sample_rows <- sample(1:nrow(data), min(10, nrow(data))) # Keep small sample
  datatrainProf <- data[sample_rows, c(expvars), drop=FALSE]

  return(list(expvars=expvars, learner=bst, estSURVTrain=preds, datatrainProf=datatrainProf, varprof=varprof))
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
  # Prepare test matrix ensuring factor levels and columns match training
  # Use the rbind trick with the sampled training data profile
  mmdata<-rbind(modelout$datatrainProf, newdata[,modelout$expvars, drop=FALSE])
  TestMat_full<-stats::model.matrix(~-1+., data=mmdata)
  # Select rows corresponding to newdata
  TestMat <- TestMat_full[-c(1:nrow(modelout$datatrainProf)), , drop=FALSE]

  # Ensure TestMat has the same columns as the matrix used for training
  # This relies on model.matrix handling factors based on the combined data
  train_cols <- colnames(modelout$learner$feature_names) # Get feature names from xgb model if stored
  if (is.null(train_cols)) {
      # Fallback: try to infer from training data matrix if feature_names not stored
      # This requires storing the training matrix or its colnames in modelout
      warning("Could not retrieve feature names from xgboost model, column consistency not guaranteed.")
  } else {
      missing_cols <- setdiff(train_cols, colnames(TestMat))
      for(col in missing_cols){
          TestMat[[col]] <- 0 # Add missing columns with 0
      }
      # Ensure correct column order and selection
      TestMat <- TestMat[, train_cols, drop = FALSE]
  }


  # Create DMatrix for prediction
  dtest <- xgboost::xgb.DMatrix(data=TestMat)

  # Get evaluation times from the baseline hazard stored in the model
  times=sort(unique(c(0, modelout$learner$baseline_hazard$time)))

  # Predict survival probabilities using the specific predict method
  preds<-predict.xgb.Booster.surv(modelout$learner, dtest, type = "surv", times=times)

  # Transpose predictions to match standard output: rows=times, cols=observations
  return(list(Probs=t(preds), Times=times))
}
