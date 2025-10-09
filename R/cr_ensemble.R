#' @title RunCRModels
#' @description  utility function to run several CR models at the same time
#' @param datatrain data frame with explanatory and outcome variables
#' @param ExpVars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded as 0,1,2),
#' 1 is the event of interest
#' @param models  a vector of models to be fitted to be chosen among
#' "FG", "rulefit", "bart", "cox", "xgboost", "gam", "survreg". Two random forest models are
#' automatically fitted and dont need to be listed here
#' @param ntreeRF number of trees for Random Forest models
#' @param varsel  logical indicating whether variable selection to be
#' applied before fitting models "FG", "rulefit", "bart", "cox", "xgboost", "gam", "survreg"
#' @return a list of two items. First is a list containing ExpVars: all of the explanatory vars,
#'  ExpVars2: a subset of the explanatory variables selected by random forest,
#'  timevar: time variable,
#'  eventvar: event variable,
#'  model_status: Named logical vector indicating which models succeeded (TRUE) or failed (FALSE).
#'  The second is also a list that contain the individual model outputs.
#' @export
RunCRModels<-function(datatrain, ExpVars, timevar, eventvar, models=c("FG", "rulefit", "bart", "cox"), ntreeRF=300, varsel=FALSE){

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

  datatrainFact<-datatrain

  for (i in seq_len(ncol(datatrainFact))){
    if (is.character(datatrainFact[,i])){
      datatrainFact[,i]<-as.factor(datatrainFact[,i])
    }
  }
  # Assuming CRModel_RF is loaded/available
  RF_Model<-tryCatch(
    CRModel_RF(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
    error = function(e) {
      message("Failed fitting RF_Model (all vars): ", e$message)
      return(NULL)
    }
  )
  model_status["RF_Model"] <- !is.null(RF_Model)

  # Get variable importance - use the first event's importance
  if (!is.null(RF_Model) && !is.null(RF_Model$rf_model$importance)) {
    imp_scores <- RF_Model$rf_model$importance[,1]
    # Ensure we have at least some variables selected
    if (all(is.na(imp_scores)) || all(imp_scores <= 0, na.rm = TRUE)) {
      # If all importances are NA or <= 0, use all variables
      ExpVars2 <- ExpVars
    } else {
      # Sort by importance and take top variables (at least 1, at most half)
      # Replace 1:min(...) with seq_len(min(...)) to handle edge cases
      n_vars <- min(length(ExpVars), max(1, length(ExpVars) %/% 2))
      ExpVars2 <- names(sort(imp_scores, decreasing = TRUE, na.last = TRUE))[seq_len(n_vars)]
    }
  } else {
    warning("Could not get importance from RF_Model, using all variables for RF_Model2.")
    ExpVars2 <- ExpVars # Fallback to all variables
  }

  RF_Model2<-tryCatch(
    CRModel_RF(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
    error = function(e) {
      message("Failed fitting RF_Model2 (top vars): ", e$message)
      return(NULL)
    }
  )
  model_status["RF_Model2"] <- !is.null(RF_Model2)

 if (!varsel){
  if ("FG" %in% models){
    # Assuming CRModel_FineGray is loaded/available
    FG_Model<-tryCatch(CRModel_FineGray(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      message("Failed fitting FG: ", e$message)
      return(NULL)
    })
    model_status["FG_Model"] <- !is.null(FG_Model)
  }
  if ("bart" %in% models){
    # Assuming CRModel_BART is loaded/available
    BART_Model<-tryCatch(CRModel_BART(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, ntree=ntreeRF), error=function(e){
      message("Failed fitting BART: ", e$message)
      return(NULL)
    })
    model_status["BART_Model"] <- !is.null(BART_Model)
  }
  if ("cox" %in% models){
    # Assuming CRModel_Cox is loaded/available
    Cox_Model<-tryCatch(CRModel_Cox(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      message("Failed fitting COX: ", e$message)
      return(NULL)
    })
    model_status["Cox_Model"] <- !is.null(Cox_Model)
  }

  if ("rulefit" %in% models){
    # Assuming CRModel_rulefit is loaded/available
    rulefit_Model<-tryCatch(CRModel_rulefit(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
       message("Failed fitting Rulefit: ", e$message)
      return(NULL)
    })
    model_status["rulefit_Model"] <- !is.null(rulefit_Model)
  }
  if ("xgboost" %in% models){
    # Assuming CRModel_xgboost is loaded/available
  xgboost_Model<-tryCatch(CRModel_xgboost(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, event_codes=1, nrounds=100), error=function(e){
       message("Failed fitting XGBoost: ", e$message)
      return(NULL)
    })
    model_status["xgboost_Model"] <- !is.null(xgboost_Model)
  }
  if ("gam" %in% models){
    # Assuming CRModel_GAM is loaded/available
    gam_Model<-tryCatch(CRModel_GAM(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, event_of_interest=1), error=function(e){
       message("Failed fitting GAM: ", e$message)
      return(NULL)
    })
    model_status["gam_Model"] <- !is.null(gam_Model)
  }
  if ("survreg" %in% models){
    # Assuming CRModel_SurvReg is loaded/available
    survreg_Model<-tryCatch(CRModel_SurvReg(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, dist="exponential"), error=function(e){
       message("Failed fitting SurvReg: ", e$message)
      return(NULL)
    })
    model_status["survreg_Model"] <- !is.null(survreg_Model)
  }
 } else {
   if ("FG" %in% models){
     FG_Model<-tryCatch(CRModel_FineGray(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       message("Failed fitting FG: ", e$message)
       return(NULL)
     })
     model_status["FG_Model"] <- !is.null(FG_Model)
   }
   if ("bart" %in% models){
     BART_Model<-tryCatch(CRModel_BART(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, ntree=ntreeRF), error=function(e){
       message("Failed fitting BART: ", e$message)
       return(NULL)
     })
     model_status["BART_Model"] <- !is.null(BART_Model)
   }
   if ("cox" %in% models){
     Cox_Model<-tryCatch(CRModel_Cox(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       message("Failed fitting COX: ", e$message)
       return(NULL)
     })
     model_status["Cox_Model"] <- !is.null(Cox_Model)
   }

   if ("rulefit" %in% models){
     rulefit_Model<-tryCatch(CRModel_rulefit(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       message("Failed fitting Rulefit: ", e$message)
       return(NULL)
     })
     model_status["rulefit_Model"] <- !is.null(rulefit_Model)
   }
   if ("xgboost" %in% models){
  xgboost_Model<-tryCatch(CRModel_xgboost(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, event_codes=1, nrounds=100), error=function(e){
       message("Failed fitting XGBoost: ", e$message)
       return(NULL)
     })
     model_status["xgboost_Model"] <- !is.null(xgboost_Model)
   }
   if ("gam" %in% models){
     gam_Model<-tryCatch(CRModel_GAM(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, event_of_interest=1), error=function(e){
       message("Failed fitting GAM: ", e$message)
       return(NULL)
     })
     model_status["gam_Model"] <- !is.null(gam_Model)
   }
   if ("survreg" %in% models){
     survreg_Model<-tryCatch(CRModel_SurvReg(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, event_of_interest=1, dist="exponential"), error=function(e){
       message("Failed fitting SurvReg: ", e$message)
       return(NULL)
     })
     model_status["survreg_Model"] <- !is.null(survreg_Model)
   }
 }



  input<-list(ExpVars=ExpVars,ExpVars2=ExpVars2, timevar=timevar, eventvar=eventvar)
  input2<-vector(mode="list")
  input2<-c(input2,list(RF_Model=RF_Model,RF_Model2=RF_Model2))

  if ("FG" %in% models){
    input2<-c(input2,list(FG_Model=FG_Model))
  }

  if ("bart" %in% models){
    input2<-c(input2,list(BART_Model=BART_Model))
  }

  if ("cox" %in% models){
    input2<-c(input2,list(Cox_Model=Cox_Model))
  }
  if ("rulefit" %in% models){
    input2<-c(input2,list(rulefit_Model=rulefit_Model))
  }
  if ("xgboost" %in% models){
    input2<-c(input2,list(xgboost_Model=xgboost_Model))
  }
  if ("gam" %in% models){
    input2<-c(input2,list(gam_Model=gam_Model))
  }
  if ("survreg" %in% models){
    input2<-c(input2,list(survreg_Model=survreg_Model))
  }

  # Print summary of model fitting
  n_success <- sum(model_status)
  n_total <- length(model_status)
  message(sprintf("Model fitting complete: %d/%d models succeeded", n_success, n_total))
  if (n_success < n_total) {
    failed_models <- names(model_status)[!model_status]
    message("Failed models: ", paste(failed_models, collapse=", "))
  }

  # Combine and convert to S3 class
  result <- c(list(input=input, model_status=model_status), input2)
  class(result) <- c("CREnsemble", "list")

  return(result)
}


#' @title ComputeCRSuperLearnerWeights
#' @description Compute and store super learner weights for a fitted CR ensemble
#' @param ensemble_models Fitted CR ensemble models (output from RunCRModels)
#' @param training_data Training data used to fit the models
#' @param eval_times Time points for weight optimization
#' @param loss_type Loss function type ("mse" or "loglik")
#' @return Updated ensemble object with super learner weights
#' @export
ComputeCRSuperLearnerWeights <- function(ensemble_models, training_data, eval_times = NULL, loss_type = "mse") {
  if (!inherits(ensemble_models, "CREnsemble")) {
    stop("ensemble_models must be output from RunCRModels")
  }
  
  if (is.null(eval_times)) {
    # Use quantiles of observed times as default
    times_col <- ensemble_models$input[[ensemble_models$input$timevar]]
    if (!is.null(training_data[[times_col]])) {
      eval_times <- quantile(training_data[[times_col]], probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    } else {
      stop("eval_times must be provided or training_data must contain the time variable")
    }
  }
  
  # Get training predictions from all successful models
  message("Computing CR super learner weights from training data...")
  training_preds <- tryCatch({
    PredictCRModels(
      models = ensemble_models,
      newdata = training_data,
      newtimes = eval_times,
      ensemble_method = "average"  # Just to get individual predictions
    )$ModelPredictions
  }, error = function(e) {
    stop("Failed to compute training predictions: ", e$message)
  })
  
  # Build observed CIF matrix
  observed_cif <- tryCatch({
    buildObservedCIFMatrix(
      data = training_data,
      timevar = ensemble_models$input$timevar,
      eventvar = ensemble_models$input$eventvar,
      cause_of_interest = 1,  # Assuming event 1 is of interest
      eval_times = eval_times
    )
  }, error = function(e) {
    stop("Failed to build observed CIF matrix: ", e$message)
  })
  
  # Optimize weights
  optimal_weights <- tryCatch({
    optimizeSuperLearnerWeights(
      predictions_list = training_preds,
      actual_surv = observed_cif,
      loss_type = loss_type
    )
  }, error = function(e) {
    stop("Super learner weight optimization failed: ", e$message)
  })
  
  message("CR super learner weights computed: ", paste(names(optimal_weights), "=", 
                                                      round(optimal_weights, 3), collapse = ", "))
  
  # Store weights in ensemble object
  ensemble_models$super_learner_weights <- optimal_weights
  ensemble_models$super_learner_eval_times <- eval_times
  ensemble_models$super_learner_loss_type <- loss_type
  
  return(ensemble_models)
}


#' @title  PredictCRModels
#' @description  Utility function to get predictions from  several competing risks models
#' at the same time. Also adds an ensemble prediction for the probabilities.
#' @param models the output from 'RunCRModels'
#' @param newdata the data for which the predictions are to be calculated
#' @param newtimes the times for which the predictions obtained
#' @param models_to_use Character vector of model names to use for ensemble prediction.
#'   If NULL (default), uses all successfully fitted models. Model names should match
#'   those in the models list (e.g., "RF_Model", "FG_Model", "Cox_Model").
#' @param ensemble_method Character specifying ensemble method: "average" (default),
#'   "weighted", or "super_learner".
#' @param model_weights Named numeric vector of weights for weighted ensemble method.
#'   Must sum to 1 (or will be normalized). Only used when ensemble_method="weighted" or
#'   when ensemble_method="super_learner" with pre-computed weights.
#' @param super_learner_training_data Data frame with training data for super learner
#'   weight optimization. Required when ensemble_method="super_learner" and model_weights=NULL
#'   and weights are not pre-computed.
#' @param super_learner_timevar Name of time variable in super_learner_training_data.
#' @param super_learner_eventvar Name of event variable in super_learner_training_data.
#' @return a list containing the following objects.
#' ModelPredictions: a list of individual model predictions,
#' NewProbs: ensembled predictions for CIFs,
#' models_used: Character vector of model names actually used in the ensemble,
#' ensemble_method: The ensemble method used.
#'
#' @export
PredictCRModels<-function(models, newdata, newtimes, models_to_use=NULL,
                           ensemble_method="average", model_weights=NULL,
                           super_learner_training_data=NULL, 
                           super_learner_timevar=NULL, 
                           super_learner_eventvar=NULL){

  if (missing(models) || is.null(models) || !is.list(models)) {
    stop("'models' must be a list returned by RunCRModels")
  }
  if (missing(newdata)) {
    stop("'newdata' must be provided")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame")
  }
  if (missing(newtimes)) {
    stop("'newtimes' must be provided")
  }
  if (!is.numeric(newtimes)) {
    stop("'newtimes' must be numeric")
  }

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
  successful_models <- names(model_status)[model_status]

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

  # --- Predictions ---
  if ("RF_Model" %in% active_models && !is.null(models$RF_Model)){
    # Assuming Predict_CRModel_RF is loaded/available
    Predict_RF<-tryCatch(
      Predict_CRModel_RF(models$RF_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for RF_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("RF_Model2" %in% active_models && !is.null(models$RF_Model2)){
    Predict_RF2<-tryCatch(
      Predict_CRModel_RF(models$RF_Model2, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for RF_Model2: ", e$message)
        return(NULL)
      }
    )
  }
  if ("FG_Model" %in% active_models && !is.null(models$FG_Model)){
    # Assuming Predict_CRModel_FineGray is loaded/available
    Predict_FG<-tryCatch(
      Predict_CRModel_FineGray(models$FG_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for FG_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("BART_Model" %in% active_models && !is.null(models$BART_Model)){
    # Assuming Predict_CRModel_BART is loaded/available
    Predict_BART<-tryCatch(
      Predict_CRModel_BART(models$BART_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for BART_Model: ", e$message)
        return(NULL)
      }
    )
  }

  if ("Cox_Model" %in% active_models && !is.null(models$Cox_Model)){
    # Assuming Predict_CRModel_Cox is loaded/available
    Predict_Cox<-tryCatch(
      Predict_CRModel_Cox(models$Cox_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for Cox_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("rulefit_Model" %in% active_models && !is.null(models$rulefit_Model)){
    # Assuming Predict_CRModel_rulefit is loaded/available
    Predict_rulefit<-tryCatch(
      Predict_CRModel_rulefit(models$rulefit_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for rulefit_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("xgboost_Model" %in% active_models && !is.null(models$xgboost_Model)){
    # Assuming Predict_CRModel_xgboost is loaded/available
    Predict_xgboost<-tryCatch(
      Predict_CRModel_xgboost(models$xgboost_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for xgboost_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("gam_Model" %in% active_models && !is.null(models$gam_Model)){
    # Assuming Predict_CRModel_GAM is loaded/available
    Predict_gam<-tryCatch(
      Predict_CRModel_GAM(models$gam_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for gam_Model: ", e$message)
        return(NULL)
      }
    )
  }
  if ("survreg_Model" %in% active_models && !is.null(models$survreg_Model)){
    # Assuming Predict_CRModel_SurvReg is loaded/available
    Predict_survreg<-tryCatch(
      Predict_CRModel_SurvReg(models$survreg_Model, newdata=newdata),
      error = function(e) {
        warning("Prediction failed for survreg_Model: ", e$message)
        return(NULL)
      }
    )
  }

  # Assuming cifMatInterpolaltor is loaded/available
  ModelPredictions<-list()
  if ("RF_Model" %in% active_models && !is.null(models$RF_Model) && exists("Predict_RF") && !is.null(Predict_RF)){
    newprobsRF<-cifMatInterpolaltor(probsMat=t(Predict_RF$CIFs),times=Predict_RF$Times, newtimes=newtimes)
    if (!is.null(newprobsRF)) {
      ModelPredictions[["RF_Model"]] <- newprobsRF
    }
  }
  if ("RF_Model2" %in% active_models && !is.null(models$RF_Model2) && exists("Predict_RF2") && !is.null(Predict_RF2)){
    newprobsRF2<-cifMatInterpolaltor(probsMat=t(Predict_RF2$CIFs),times=Predict_RF2$Times, newtimes=newtimes)
    if (!is.null(newprobsRF2)) {
      ModelPredictions[["RF_Model2"]] <- newprobsRF2
    }
  }
  if ("FG_Model" %in% active_models && !is.null(models$FG_Model) && exists("Predict_FG") && !is.null(Predict_FG)){
    newprobsFG<-cifMatInterpolaltor(probsMat=t(Predict_FG$CIFs),times=Predict_FG$Times, newtimes=newtimes)
    if (!is.null(newprobsFG)) {
      ModelPredictions[["FG_Model"]] <- newprobsFG
    }
  }
  if ("BART_Model" %in% active_models && !is.null(models$BART_Model) && exists("Predict_BART") && !is.null(Predict_BART)){
    newprobsBART<-cifMatInterpolaltor(probsMat=t(Predict_BART$CIFs),times=Predict_BART$Times, newtimes=newtimes)
    if (!is.null(newprobsBART)) {
      ModelPredictions[["BART_Model"]] <- newprobsBART
    }
  }
  if ("Cox_Model" %in% active_models && !is.null(models$Cox_Model) && exists("Predict_Cox") && !is.null(Predict_Cox)){
    newprobsCox<-cifMatInterpolaltor(probsMat=t(Predict_Cox$CIFs),times=Predict_Cox$Times, newtimes=newtimes)
    if (!is.null(newprobsCox)) {
      ModelPredictions[["Cox_Model"]] <- newprobsCox
    }
  }
  if ("rulefit_Model" %in% active_models && !is.null(models$rulefit_Model) && exists("Predict_rulefit") && !is.null(Predict_rulefit)){
    newprobsrulefit<-cifMatInterpolaltor(probsMat=t(Predict_rulefit$CIFs),times=Predict_rulefit$Times, newtimes=newtimes)
    if (!is.null(newprobsrulefit)) {
      ModelPredictions[["rulefit_Model"]] <- newprobsrulefit
    }
  }
  if ("xgboost_Model" %in% active_models && !is.null(models$xgboost_Model) && exists("Predict_xgboost") && !is.null(Predict_xgboost)){
    newprobsxgboost<-cifMatInterpolaltor(probsMat=t(Predict_xgboost$CIFs),times=Predict_xgboost$Times, newtimes=newtimes)
    if (!is.null(newprobsxgboost)) {
      ModelPredictions[["xgboost_Model"]] <- newprobsxgboost
    }
  }
  if ("gam_Model" %in% active_models && !is.null(models$gam_Model) && exists("Predict_gam") && !is.null(Predict_gam)){
    newprobsgam<-cifMatInterpolaltor(probsMat=t(Predict_gam$CIFs),times=Predict_gam$Times, newtimes=newtimes)
    if (!is.null(newprobsgam)) {
      ModelPredictions[["gam_Model"]] <- newprobsgam
    }
  }
  if ("survreg_Model" %in% active_models && !is.null(models$survreg_Model) && exists("Predict_survreg") && !is.null(Predict_survreg)){
    newprobsurvreg<-cifMatInterpolaltor(probsMat=t(Predict_survreg$CIFs),times=Predict_survreg$Times, newtimes=newtimes)
    if (!is.null(newprobsurvreg)) {
      ModelPredictions[["survreg_Model"]] <- newprobsurvreg
    }
  }

  # Remove any NULL predictions (safety)
  ModelPredictions <- Filter(Negate(is.null), ModelPredictions)

  # Check if we have any valid predictions
  if (length(ModelPredictions) == 0) {
    warning("No valid predictions obtained from any model.")
    return(list(ModelPredictions=ModelPredictions, NewProbs=NULL, models_used=character(0),
                ensemble_method=ensemble_method))
  }

  models_used <- names(ModelPredictions)

  # Apply the selected ensemble method
  NewProbs <- tryCatch(
    {
      if (ensemble_method == "average") {
        cifMatListAveraging(ModelPredictions, type = "prob")
      } else if (ensemble_method == "weighted") {
        if (is.null(model_weights)) {
          stop("model_weights must be provided for weighted ensemble method")
        }
        cifMatWeightedAveraging(ModelPredictions, model_weights, type = "prob")
      } else if (ensemble_method == "super_learner") {
        # Super learner implementation with proper stacking
        if (is.null(model_weights)) {
          # Check if weights are already stored in the ensemble
          if (!is.null(models$super_learner_weights)) {
            stored_weights <- models$super_learner_weights
            # Filter weights for available models
            available_weights <- stored_weights[intersect(names(stored_weights), names(ModelPredictions))]
            if (length(available_weights) > 0) {
              # Normalize weights for available models
              available_weights <- available_weights / sum(available_weights)
              message("Using pre-computed CR super learner weights: ", paste(names(available_weights), "=", 
                                                                           round(available_weights, 3), collapse = ", "))
              cifMatWeightedAveraging(ModelPredictions, available_weights, type = "prob")
            } else {
              message("No pre-computed weights available for current models, falling back to equal-weight averaging")
              cifMatListAveraging(ModelPredictions, type = "prob")
            }
          } else if (!is.null(super_learner_training_data) && 
                     !is.null(super_learner_timevar) && 
                     !is.null(super_learner_eventvar)) {
            # Compute weights on-the-fly (backward compatibility)
            message("Computing CR super learner weights using provided training data...")
            warning("Consider using ComputeCRSuperLearnerWeights() during training for better performance")
            
            # Get training predictions from same models
            training_preds <- tryCatch({
              PredictCRModels(
                models = models,
                newdata = super_learner_training_data,
                newtimes = newtimes,
                models_to_use = names(ModelPredictions),
                ensemble_method = "average"  # Just to get individual predictions
              )$ModelPredictions
            }, error = function(e) {
              warning("Failed to compute training predictions for CR super learner: ", e$message)
              return(NULL)
            })
            
            if (!is.null(training_preds)) {
              # Build observed CIF matrix
              observed_cif <- tryCatch({
                buildObservedCIFMatrix(
                  data = super_learner_training_data,
                  timevar = super_learner_timevar,
                  eventvar = super_learner_eventvar,
                  cause_of_interest = 1,  # Assuming event 1 is of interest
                  eval_times = newtimes
                )
              }, error = function(e) {
                warning("Failed to build observed CIF matrix: ", e$message)
                return(NULL)
              })
              
              if (!is.null(observed_cif)) {
                # Optimize weights
                optimal_weights <- tryCatch({
                  optimizeSuperLearnerWeights(
                    predictions_list = training_preds,
                    actual_surv = observed_cif,
                    loss_type = "mse"
                  )
                }, error = function(e) {
                  warning("CR super learner weight optimization failed: ", e$message)
                  return(NULL)
                })
                
                if (!is.null(optimal_weights)) {
                  message("CR super learner weights: ", paste(names(optimal_weights), "=", 
                                                             round(optimal_weights, 3), collapse = ", "))
                  cifMatWeightedAveraging(ModelPredictions, optimal_weights, type = "prob")
                } else {
                  message("Falling back to equal-weight averaging")
                  cifMatListAveraging(ModelPredictions, type = "prob")
                }
              } else {
                message("Falling back to equal-weight averaging")
                cifMatListAveraging(ModelPredictions, type = "prob")
              }
            } else {
              message("Falling back to equal-weight averaging")
              cifMatListAveraging(ModelPredictions, type = "prob")
            }
          } else {
            message("CR super learner requires either pre-computed weights, model_weights, or training data parameters")
            message("Use ComputeCRSuperLearnerWeights() to pre-compute weights during training")
            message("Falling back to equal-weight averaging")
            cifMatListAveraging(ModelPredictions, type = "prob")
          }
        } else {
          # Use provided super learner weights
          message("Using provided CR super learner weights: ", paste(names(model_weights), "=", 
                                                                     round(model_weights, 3), collapse = ", "))
          cifMatWeightedAveraging(ModelPredictions, model_weights, type = "prob")
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

  if (!is.null(NewProbs)) {
    NewProbs <- rbind(NewProbs) # Ensure it's a matrix
  }

  return(list(ModelPredictions=ModelPredictions, NewProbs=NewProbs, models_used=models_used,
              ensemble_method=ensemble_method))
}

cifMatListAveraging <- function(list_mats, type = "prob", na.rm = TRUE) {
  if (length(list_mats) == 0) {
    return(NULL)
  }

  # Check for consistent dimensions
  dims <- lapply(list_mats, dim)
  if (length(unique(dims)) > 1) {
    stop("All matrices in the list must have the same dimensions")
  }

  # Handle NA values
  if (na.rm) {
    averaged_mat <- Reduce("+", lapply(list_mats, function(mat) {
      mat[is.na(mat)] <- 0
      mat
    })) / length(list_mats)
  } else {
    averaged_mat <- Reduce("+", list_mats) / length(list_mats)
  }

  return(averaged_mat)
}
