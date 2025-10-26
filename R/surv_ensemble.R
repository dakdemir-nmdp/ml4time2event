#' @title RunSurvModels
#'
#' @description Utility function to run several survival models at the same time.
#' Fits two Random Forest models (one with all vars, one with top N vars) automatically.
#'
#' @param datatrain data frame with training data (explanatory and outcome variables)
#' @param ExpVars character vector of names of all explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param models character vector of additional models to fit, chosen from:
#' "glmnet" (penalized Cox with elastic net), "coxph" (Cox with backward selection),
#' "rulefit", "xgboost", "gam", "gbm", "ExpSurvReg", "WeibSurvReg", "bart".
#' @param ntreeRF number of trees for Random Forest models.
#' @param nvars number of top variables (based on RF importance) to use for the second RF model and potentially others.
#' @param run_rf logical; if `FALSE`, skips fitting the baseline Random Forest models.
#' @param ... Additional arguments passed to individual model fitting functions (currently not explicitly handled, consider adding specific args).
#' @return A list containing:
#'   - input: list with ExpVars, ExpVars2 (top N vars), timevar, eventvar.
#'   - RF_Model: Output from SurvModel_RF using all ExpVars.
#'   - RF_Model2: Output from SurvModel_RF using top ExpVars2.
#'   - model_status: Named logical vector indicating which models succeeded (TRUE) or failed (FALSE).
#'   - ... other fitted model outputs named according to the 'models' input (e.g., glmnet_Model, CPH_Model).
#'
#' @export
RunSurvModels<-function(datatrain, ExpVars, timevar, eventvar,
                        models=c("glmnet","coxph","rulefit","xgboost","gam","gbm","ExpSurvReg","WeibSurvReg","bart","deepsurv"),
                        ntreeRF=300, nvars=20, run_rf = TRUE, ...){

  if (missing(datatrain) || is.null(datatrain) || !is.data.frame(datatrain)) {
    stop("'datatrain' must be a data frame")
  }
  if (missing(ExpVars) || length(ExpVars) == 0) {
    stop("'ExpVars' must be a non-empty character vector")
  }
  if (!all(ExpVars %in% colnames(datatrain))) {
    stop("All ExpVars must be present in datatrain")
  }
  if (missing(timevar) || !timevar %in% colnames(datatrain)) {
    stop("'timevar' must be a column in datatrain")
  }
  if (missing(eventvar) || !eventvar %in% colnames(datatrain)) {
    stop("'eventvar' must be a column in datatrain")
  }
  if (!is.null(models) && !is.character(models)) {
    stop("'models' must be a character vector")
  }

  # Initialize model status tracker
  model_status <- c()

  # Ensure event var is 0/1
  datatrain[[eventvar]] <- as.numeric(datatrain[[eventvar]] == 1)

  # Convert character columns to factors (needed by some models like RF)
  datatrainFact<-datatrain
  for (i in ExpVars){ # Iterate through predictors only
    if (is.character(datatrainFact[[i]])){
      datatrainFact[[i]]<-as.factor(datatrainFact[[i]])
    }
  }

  ExpVars2 <- ExpVars
  RF_Model <- NULL
  RF_Model2 <- NULL

  if (isTRUE(run_rf)) {
    # --- Fit Base Random Forest Models ---
    # Assuming SurvModel_RF is loaded/available
    RF_Model<-tryCatch(
        SurvModel_RF(data=datatrainFact, expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
        error = function(e) {
          message("Failed fitting RF_Model (all vars): ", e$message)
          return(NULL)
        }
    )
    model_status["RF_Model"] <- !is.null(RF_Model)

    if (!is.null(RF_Model)) {
        got_importance <- FALSE

        if (!got_importance && !is.null(RF_Model$model) && !is.null(RF_Model$model$importance) && !all(is.na(RF_Model$model$importance))) {
            importance_scores <- RF_Model$model$importance
            if (!is.null(importance_scores) && (is.matrix(importance_scores) || is.vector(importance_scores))) {
                if (is.matrix(importance_scores)) {
                    ExpVars2 <- names(sort(importance_scores[,1], decreasing=TRUE))[seq_len(min(length(ExpVars), nvars))]
                } else {
                    ExpVars2 <- names(sort(importance_scores, decreasing=TRUE))[seq_len(min(length(ExpVars), nvars))]
                }
                got_importance <- TRUE
            }
        }

        if (!got_importance && requireNamespace("randomForestSRC", quietly = TRUE) && !is.null(RF_Model$model)) {
            tryCatch({
                imp_obj <- randomForestSRC::vimp(RF_Model$model)
                if (!is.null(imp_obj) && !is.null(imp_obj$importance)) {
                    importance_scores <- imp_obj$importance
                    if (is.matrix(importance_scores)) {
                        ExpVars2 <- names(sort(importance_scores[,1], decreasing=TRUE))[seq_len(min(length(ExpVars), nvars))]
                    } else {
                        ExpVars2 <- names(sort(importance_scores, decreasing=TRUE))[seq_len(min(length(ExpVars), nvars))]
                    }
                    got_importance <- TRUE
                }
            }, error = function(e) {
                # vimp function failed, continue to next method
            })
        }

        if (!got_importance && !is.null(RF_Model$model$var.used) && length(RF_Model$model$var.used) > 0) {
            var_freq <- table(RF_Model$model$var.used)
            var_names <- names(var_freq)
            if (length(var_names) > 0) {
                ExpVars2 <- names(sort(var_freq, decreasing=TRUE))[seq_len(min(length(var_names), nvars))]
                got_importance <- TRUE
            }
        }

        if (!got_importance) {
            warning("Could not get importance from RF_Model, using all variables for RF_Model2.")
            ExpVars2 <- ExpVars
        }
    } else {
        warning("RF_Model is NULL, using all variables for RF_Model2.")
        ExpVars2 <- ExpVars
    }

    RF_Model2<-tryCatch(
        SurvModel_RF(data=datatrainFact, expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
        error = function(e) {
          message("Failed fitting RF_Model2 (top vars): ", e$message)
          return(NULL)
        }
    )
    model_status["RF_Model2"] <- !is.null(RF_Model2)
  } else {
    message("Skipping baseline Random Forest models (run_rf = FALSE).")
    model_status["RF_Model"] <- NA
    model_status["RF_Model2"] <- NA
  }


  # --- Fit Additional Requested Models ---
  # Use ExpVars2 (top variables) for these models as in original logic
  ExpVarsForOthers <- ExpVars2

  # Initialize list to store model outputs
  input2<-vector(mode="list")
  if (isTRUE(run_rf)) {
    input2<-c(input2,list(RF_Model=RF_Model, RF_Model2=RF_Model2)) # Add base RF models first
  }

  # Helper function for fitting with error catching
  fit_model <- function(model_name, fit_func, ...) {
      cat("Fitting", model_name, "...\n")
      model_out <- tryCatch(
          fit_func(...),
          error = function(e) {
              message("Failed fitting ", model_name, ": ", e$message)
              return(NULL)
          }
      )
      return(model_out)
  }

  if ("bart"%in% models){
    # Assuming SurvModel_BART is loaded/available
    bartout<-fit_model("bart", SurvModel_BART, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, ntree=50) # Reduced ntree for speed?
    input2<-c(input2,list(bart_Model=bartout))
    model_status["bart_Model"] <- !is.null(bartout)
  }

  if ("deepsurv"%in% models){
    # Assuming SurvModel_DeepSurv is loaded/available
    deepsurv_Model<-fit_model("deepsurv", SurvModel_DeepSurv, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, size=5, decay=0.01, maxit=500)
    input2<-c(input2,list(deepsurv_Model=deepsurv_Model))
    model_status["deepsurv_Model"] <- !is.null(deepsurv_Model)
  }

  if ("ExpSurvReg" %in% models) {
    # Assuming SurvModel_SurvReg is loaded/available
    # Original code had a loop reducing variables if fit failed, keeping similar logic
    ExpVarsW <- ExpVarsForOthers
    survregexp_Model <- NULL
    while (length(ExpVarsW) > 0 && is.null(survregexp_Model)) {
      survregexp_Model <- fit_model("ExpSurvReg", SurvModel_SurvReg, data=datatrainFact, expvars=ExpVarsW, timevar=timevar, eventvar=eventvar, dist = "exponential")
      if (is.null(survregexp_Model) && length(ExpVarsW) > 1) {
          warning("ExpSurvReg failed with ", length(ExpVarsW), " vars, removing last one and retrying.")
          ExpVarsW <- ExpVarsW[-length(ExpVarsW)] # Remove last variable and retry
      } else if (is.null(survregexp_Model)) {
          break # Stop if only 1 var left and it failed
      }
    }
    input2<-c(input2,list(survregexp_Model=survregexp_Model))
    model_status["survregexp_Model"] <- !is.null(survregexp_Model)
  }

  if ("WeibSurvReg" %in% models) {
    # Assuming SurvModel_SurvReg is loaded/available
    ExpVarsW <- ExpVarsForOthers
    survregweib_Model <- NULL
    while (length(ExpVarsW) > 0 && is.null(survregweib_Model)) {
      survregweib_Model <- fit_model("WeibSurvReg", SurvModel_SurvReg, data=datatrainFact, expvars=ExpVarsW, timevar=timevar, eventvar=eventvar, dist = "weibull")
       if (is.null(survregweib_Model) && length(ExpVarsW) > 1) {
          warning("WeibSurvReg failed with ", length(ExpVarsW), " vars, removing last one and retrying.")
          ExpVarsW <- ExpVarsW[-length(ExpVarsW)]
      } else if (is.null(survregweib_Model)) {
          break
      }
    }
     input2<-c(input2,list(survregweib_Model=survregweib_Model))
     model_status["survregweib_Model"] <- !is.null(survregweib_Model)
  }


  if ("glmnet" %in% models){
    # Use penalized Cox (elastic net) from SurvModel_Cox
    glmnet_Model <- fit_model("glmnet",
                              SurvModel_Cox,
                              data = datatrainFact,
                              expvars = ExpVarsForOthers,
                              timevar = timevar,
                              eventvar = eventvar,
                              varsel = "penalized",
                              alpha = 0.5,  # Elastic net (can be passed via ...)
                              nfolds = 10,
                              verbose = FALSE)
    input2 <- c(input2, list(glmnet_Model = glmnet_Model))
    model_status["glmnet_Model"] <- !is.null(glmnet_Model)
  }
  if ("coxph" %in% models){
    # Standard Cox with backward selection
    CPH_Model <- fit_model("coxph",
                          SurvModel_Cox,
                          data = datatrainFact,
                          expvars = ExpVarsForOthers,
                          timevar = timevar,
                          eventvar = eventvar,
                          varsel = "backward",
                          penalty = "AIC",
                          verbose = FALSE)
    input2 <- c(input2, list(CPH_Model = CPH_Model))
    model_status["CPH_Model"] <- !is.null(CPH_Model)
  }
  if ("rulefit"%in% models){
    # Assuming SurvModel_rulefit is loaded/available
    # Rulefit often benefits from more variables, maybe use ExpVars instead of ExpVarsForOthers? Sticking to original logic for now.
    RuleFit_Model<-fit_model("rulefit", SurvModel_rulefit, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, ntree=ntreeRF, nsample=min(c(500,ceiling(.3*nrow(datatrainFact)))))
    input2<-c(input2,list(RuleFit_Model=RuleFit_Model))
    model_status["RuleFit_Model"] <- !is.null(RuleFit_Model)
  }
  if ("xgboost"%in% models){
    # Assuming SurvModel_xgboost is loaded/available
    xgboost_Model<-fit_model("xgboost", SurvModel_xgboost, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(xgboost_Model=xgboost_Model))
    model_status["xgboost_Model"] <- !is.null(xgboost_Model)
  }
  if ("gam"%in% models){
     # Assuming SurvModel_GAM is loaded/available
     # Original code had fallback to ExpVars2 if ExpVars failed, simplifying here
    gam_Model<-fit_model("gam", SurvModel_GAM, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(gam_Model=gam_Model))
    model_status["gam_Model"] <- !is.null(gam_Model)
  }
  if ("gbm"%in% models){
    # Assuming SurvModel_gbm is loaded/available
    gbm_Model<-fit_model("gbm", SurvModel_gbm, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(gbm_Model=gbm_Model))
    model_status["gbm_Model"] <- !is.null(gbm_Model)
  }

  # --- Return Results ---
  input_params<-list(ExpVars=ExpVars, ExpVars2=ExpVars2, timevar=timevar, eventvar=eventvar)

  # Print summary of model fitting
  valid_status <- !is.na(model_status)
  n_total <- sum(valid_status)
  n_success <- if (n_total > 0) sum(model_status[valid_status]) else 0
  message(sprintf("Model fitting complete: %d/%d models succeeded", n_success, n_total))
  failed_models <- names(model_status)[valid_status & !model_status]
  if (length(failed_models) > 0) {
    message("Failed models: ", paste(failed_models, collapse=", "))
  }

  # Combine input parameters, model status, and fitted models
  result <- c(list(input=input_params, model_status=model_status), input2)

  # Convert to S3 class
  class(result) <- c("SurvEnsemble", "list")

  return(result)
}


#' @title ComputeSuperLearnerWeights
#' @description Compute and store super learner weights for a fitted ensemble
#' @param ensemble_models Fitted ensemble models (output from RunSurvModels)
#' @param training_data Training data used to fit the models
#' @param eval_times Time points for weight optimization
#' @param loss_type Loss function type ("mse" or "loglik")
#' @return Updated ensemble object with super learner weights
#' @export
ComputeSuperLearnerWeights <- function(ensemble_models, training_data, eval_times = NULL, loss_type = "mse") {
  if (!inherits(ensemble_models, "SurvEnsemble")) {
    stop("ensemble_models must be output from RunSurvModels")
  }
  
  if (is.null(eval_times)) {
    # Use quantiles of observed times as default
    times_col <- ensemble_models$input$timevar
    if (!is.null(training_data[[times_col]])) {
      eval_times <- quantile(training_data[[times_col]], probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    } else {
      stop("eval_times must be provided or training_data must contain the time variable")
    }
  }
  
  # Get training predictions from all successful models
  message("Computing super learner weights from training data...")
  training_preds <- tryCatch({
    PredictSurvModels(
      models = ensemble_models,
      newdata = training_data,
      new_times = eval_times,
      ensemble_method = "average"  # Just to get individual predictions
    )$ModelPredictions
  }, error = function(e) {
    stop("Failed to compute training predictions: ", e$message)
  })
  
  # Build observed survival matrix
  observed_survival <- tryCatch({
    buildObservedSurvivalMatrix(
      data = training_data,
      timevar = ensemble_models$input$timevar,
      eventvar = ensemble_models$input$eventvar,
      eval_times = eval_times
    )
  }, error = function(e) {
    stop("Failed to build observed survival matrix: ", e$message)
  })
  
  # Optimize weights
  optimal_weights <- tryCatch({
    optimizeSuperLearnerWeights(
      predictions_list = training_preds,
      actual_surv = observed_survival,
      loss_type = loss_type
    )
  }, error = function(e) {
    stop("Super learner weight optimization failed: ", e$message)
  })
  
  message("Super learner weights computed: ", paste(names(optimal_weights), "=", 
                                                   round(optimal_weights, 3), collapse = ", "))
  
  # Store weights in ensemble object
  ensemble_models$super_learner_weights <- optimal_weights
  ensemble_models$super_learner_eval_times <- eval_times
  ensemble_models$super_learner_loss_type <- loss_type
  
  return(ensemble_models)
}


#' @title PredictSurvModels
#'
#' @description Utility function to get predictions from several survival models
#' at the same time and calculate an ensemble prediction.
#'
#' @param models List containing the fitted models (output from 'RunSurvModels').
#' @param newdata Data frame with new data for prediction.
#' @param new_times Numeric vector of times for which predictions are required.
#' @param models_to_use Character vector of model names to use for ensemble prediction.
#'   If NULL (default), uses all successfully fitted models. Model names should match
#'   those in the models list (e.g., "RF_Model", "glmnet_Model", "CPH_Model").
#' @param ensemble_method Character specifying ensemble method: "average" (default),
#'   "weighted", or "super_learner".
#' @param model_weights Named numeric vector of weights for weighted ensemble method.
#'   Must sum to 1 (or will be normalized). Only used when ensemble_method="weighted" or
#'   when ensemble_method="super_learner" with pre-computed weights.
#' @param super_learner_training_data Data frame with training data for super learner
#'   weight optimization. Only required when ensemble_method="super_learner", model_weights=NULL,
#'   and weights are not pre-computed using ComputeSuperLearnerWeights().
#' @param super_learner_timevar Name of time variable in super_learner_training_data.
#' @param super_learner_eventvar Name of event variable in super_learner_training_data.
#' @return A list containing:
#'   - ModelPredictions: A list of individual model predictions (interpolated probability matrices).
#'   - NewProbs: Ensembled prediction matrix (averaged on cumulative hazard scale).
#'   - models_used: Character vector of model names actually used in the ensemble.
#'   - ensemble_method: The ensemble method used.
#'
#' @export
PredictSurvModels<-function(models, newdata, new_times, models_to_use=NULL,
                             ensemble_method="average", model_weights=NULL, 
                             super_learner_training_data=NULL, 
                             super_learner_timevar=NULL, 
                             super_learner_eventvar=NULL){

  if (missing(models) || is.null(models) || !is.list(models)) {
    stop("'models' must be a list returned by RunSurvModels")
  }
  if (missing(newdata)) {
    stop("'newdata' must be provided")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame")
  }
  if (missing(new_times)) {
    stop("'new_times' must be provided")
  }
  if (!is.numeric(new_times)) {
    stop("'new_times' must be numeric")
  }

  # Ensure new_times includes 0 and is sorted
  new_times <- sort(unique(c(0, new_times)))

  # --- Model Selection and Validation ---
  # Get model status if available (backward compatibility: if not present, assume all non-NULL models succeeded)
  if (!is.null(models$model_status)) {
    model_status <- models$model_status
  } else {
    # Create status based on which models are non-NULL
    model_names <- setdiff(names(models), c("input", "model_status"))
    model_status <- sapply(model_names, function(m) !is.null(models[[m]]))
    names(model_status) <- model_names
  }

  # Get successful models
  valid_success_idx <- which(!is.na(model_status) & model_status)
  successful_models <- names(model_status)[valid_success_idx]

  # Validate and filter models_to_use
  if (!is.null(models_to_use)) {
    if (!is.character(models_to_use) || length(models_to_use) == 0) {
      stop("'models_to_use' must be a non-empty character vector")
    }

    # Check if requested models exist
    missing_models <- setdiff(models_to_use, names(models))
    if (length(missing_models) > 0) {
      unknown_models <- setdiff(missing_models, names(model_status))
      if (length(unknown_models) > 0) {
        stop("The following models were requested but not found: ",
             paste(unknown_models, collapse=", "))
      }
    }

    # Check if requested models succeeded
    failed_requested <- setdiff(models_to_use, successful_models)
    if (length(failed_requested) > 0) {
      warning("The following requested models failed during training and will be excluded: ",
              paste(failed_requested, collapse=", "))
      models_to_use <- intersect(models_to_use, successful_models)
    }

    if (length(models_to_use) == 0) {
      stop("No successfully fitted models available among those requested")
    }

    active_models <- models_to_use
  } else {
    # Use all successful models
    active_models <- successful_models
  }

  if (length(active_models) == 0) {
    stop("No successfully fitted models available for prediction")
  }

  message(sprintf("Using %d models for ensemble prediction: %s",
                  length(active_models), paste(active_models, collapse=", ")))

  # Prepare newdata (e.g., handle factors - basic handling here)
  newdataFactor<-newdata
  # A more robust factor handling would use levels from training data (e.g., stored in varprof)
  for (i in colnames(newdataFactor)){
    if (is.character(newdataFactor[[i]])){
      newdataFactor[[i]]<-as.factor(newdataFactor[[i]])
    }
  }

  # --- Get Predictions from Individual Models ---
  ModelPredictionsList <- list() # To store interpolated predictions

  # Helper function for prediction and interpolation
  predict_and_interp <- function(model_name, predict_func, model_obj) {
      if (!is.null(model_obj)) {
          cat("Predicting", model_name, "...\n")
          pred_out <- tryCatch(
              predict_func(model_obj, newdata=newdataFactor),
              error = function(e) {
                  warning("Prediction failed for ", model_name, ": ", e$message)
                  return(NULL)
              }
          )
          if (!is.null(pred_out) && !is.null(pred_out$Probs) && !is.null(pred_out$Times)) {
              # Assuming survprobMatInterpolator is loaded/available
              interp_probs <- tryCatch(
                  survprobMatInterpolator(probsMat=pred_out$Probs, times=pred_out$Times, new_times=new_times),
                   error = function(e) {
                      warning("Interpolation failed for ", model_name, ": ", e$message)
                      return(NULL)
                  }
              )
              return(interp_probs)
          }
      }
      return(NULL)
  }

  # Predict and interpolate for each active model
  # Assuming Predict_* functions are loaded/available
  if ("RuleFit_Model" %in% active_models) {
    ModelPredictionsList[["RuleFit_Model"]] <- predict_and_interp("RuleFit", Predict_SurvModel_rulefit, models$RuleFit_Model)
  }
  if ("RF_Model" %in% active_models) {
    ModelPredictionsList[["RF_Model"]] <- predict_and_interp("RF", Predict_SurvModel_RF, models$RF_Model)
  }
  if ("RF_Model2" %in% active_models) {
    ModelPredictionsList[["RF_Model2"]] <- predict_and_interp("RF2", Predict_SurvModel_RF, models$RF_Model2)
  }
  # Both glmnet and coxph now use Predict_SurvModel_Cox (unified Cox interface)
  if ("glmnet_Model" %in% active_models) {
    ModelPredictionsList[["glmnet_Model"]] <- predict_and_interp("glmnet", Predict_SurvModel_Cox, models$glmnet_Model)
  }
  if ("CPH_Model" %in% active_models) {
    ModelPredictionsList[["CPH_Model"]] <- predict_and_interp("CoxPH", Predict_SurvModel_Cox, models$CPH_Model)
  }
  if ("bart_Model" %in% active_models) {
    ModelPredictionsList[["bart_Model"]] <- predict_and_interp("BART", Predict_SurvModel_BART, models$bart_Model)
  }
  if ("deepsurv_Model" %in% active_models) {
    ModelPredictionsList[["deepsurv_Model"]] <- predict_and_interp("DeepSurv", Predict_SurvModel_DeepSurv, models$deepsurv_Model)
  }
  if ("gam_Model" %in% active_models) {
    ModelPredictionsList[["gam_Model"]] <- predict_and_interp("GAM", Predict_SurvModel_GAM, models$gam_Model)
  }
  if ("gbm_Model" %in% active_models) {
    ModelPredictionsList[["gbm_Model"]] <- predict_and_interp("GBM", Predict_SurvModel_gbm, models$gbm_Model)
  }
  if ("survregexp_Model" %in% active_models) {
    ModelPredictionsList[["survregexp_Model"]] <- predict_and_interp("ExpSurvReg", Predict_SurvModel_SurvReg, models$survregexp_Model)
  }
  if ("survregweib_Model" %in% active_models) {
    ModelPredictionsList[["survregweib_Model"]] <- predict_and_interp("WeibSurvReg", Predict_SurvModel_SurvReg, models$survregweib_Model)
  }
  if ("xgboost_Model" %in% active_models) {
    ModelPredictionsList[["xgboost_Model"]] <- predict_and_interp("XGBoost", Predict_SurvModel_xgboost, models$xgboost_Model)
  }

  # Filter out NULL predictions (from failed models/predictions/interpolations)
  ValidPredictions <- Filter(Negate(is.null), ModelPredictionsList)

  if (length(ValidPredictions) == 0) {
      warning("No valid predictions obtained from any model.")
      return(list(ModelPredictions=ValidPredictions, NewProbs=NULL, models_used=character(0),
                  ensemble_method=ensemble_method))
  }

  # --- Ensemble Predictions ---
  # Apply the selected ensemble method
  NewProbs <- tryCatch(
      {
        if (ensemble_method == "average") {
          survprobMatListAveraging(ValidPredictions)
        } else if (ensemble_method == "weighted") {
          if (is.null(model_weights)) {
            stop("model_weights must be provided for weighted ensemble method")
          }
          survprobMatWeightedAveraging(ValidPredictions, model_weights)
        } else if (ensemble_method == "super_learner") {
          # Super learner implementation with proper stacking
          if (is.null(model_weights)) {
            # Check if weights are already stored in the ensemble
            if (!is.null(models$super_learner_weights)) {
              stored_weights <- models$super_learner_weights
              # Filter weights for available models
              available_weights <- stored_weights[intersect(names(stored_weights), names(ValidPredictions))]
              if (length(available_weights) > 0) {
                # Normalize weights for available models
                available_weights <- available_weights / sum(available_weights)
                message("Using pre-computed super learner weights: ", paste(names(available_weights), "=", 
                                                                           round(available_weights, 3), collapse = ", "))
                survprobMatWeightedAveraging(ValidPredictions, available_weights)
              } else {
                message("No pre-computed weights available for current models, falling back to equal-weight averaging")
                survprobMatListAveraging(ValidPredictions)
              }
            } else if (!is.null(super_learner_training_data) && 
                       !is.null(super_learner_timevar) && 
                       !is.null(super_learner_eventvar)) {
              # Compute weights on-the-fly (backward compatibility)
              message("Computing super learner weights using provided training data...")
              warning("Consider using ComputeSuperLearnerWeights() during training for better performance")
              
              # Get training predictions from same models
              training_preds <- tryCatch({
                PredictSurvModels(
                  models = models,
                  newdata = super_learner_training_data,
                  new_times = new_times,
                  models_to_use = names(ValidPredictions),
                  ensemble_method = "average"  # Just to get individual predictions
                )$ModelPredictions
              }, error = function(e) {
                warning("Failed to compute training predictions for super learner: ", e$message)
                return(NULL)
              })
              
              if (!is.null(training_preds)) {
                # Build observed survival matrix
                observed_survival <- tryCatch({
                  buildObservedSurvivalMatrix(
                    data = super_learner_training_data,
                    timevar = super_learner_timevar,
                    eventvar = super_learner_eventvar,
                    eval_times = new_times
                  )
                }, error = function(e) {
                  warning("Failed to build observed survival matrix: ", e$message)
                  return(NULL)
                })
                
                if (!is.null(observed_survival)) {
                  # Optimize weights
                  optimal_weights <- tryCatch({
                    optimizeSuperLearnerWeights(
                      predictions_list = training_preds,
                      actual_surv = observed_survival,
                      loss_type = "mse"
                    )
                  }, error = function(e) {
                    warning("Super learner weight optimization failed: ", e$message)
                    return(NULL)
                  })
                  
                  if (!is.null(optimal_weights)) {
                    message("Super learner weights: ", paste(names(optimal_weights), "=", 
                                                            round(optimal_weights, 3), collapse = ", "))
                    survprobMatWeightedAveraging(ValidPredictions, optimal_weights)
                  } else {
                    message("Falling back to equal-weight averaging")
                    survprobMatListAveraging(ValidPredictions)
                  }
                } else {
                  message("Falling back to equal-weight averaging")
                  survprobMatListAveraging(ValidPredictions)
                }
              } else {
                message("Falling back to equal-weight averaging")
                survprobMatListAveraging(ValidPredictions)
              }
            } else {
              message("Super learner requires either pre-computed weights, model_weights, or training data parameters")
              message("Use ComputeSuperLearnerWeights() to pre-compute weights during training")
              message("Falling back to equal-weight averaging")
              survprobMatListAveraging(ValidPredictions)
            }
          } else {
            # Use provided super learner weights
            message("Using provided super learner weights: ", paste(names(model_weights), "=", 
                                                                   round(model_weights, 3), collapse = ", "))
            survprobMatWeightedAveraging(ValidPredictions, model_weights)
          }
        } else {
          stop("Unknown ensemble_method: ", ensemble_method)
        }
      },
      error = function(e) {
          warning("Ensemble averaging failed: ", e$message)
          return(NULL)
      }
  )

  models_used <- names(ValidPredictions)

  return(list(
    ModelPredictions = ValidPredictions,
    NewProbs = NewProbs,
    NewTimes = new_times,
    models_used = models_used,
    ensemble_method = ensemble_method
  ))
}
