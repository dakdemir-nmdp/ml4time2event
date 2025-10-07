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
#' @param ... Additional arguments passed to individual model fitting functions (currently not explicitly handled, consider adding specific args).
#' @return A list containing:
#'   - input: list with ExpVars, ExpVars2 (top N vars), timevar, eventvar.
#'   - RF_Model: Output from SurvModel_RF using all ExpVars.
#'   - RF_Model2: Output from SurvModel_RF using top ExpVars2.
#'   - ... other fitted model outputs named according to the 'models' input (e.g., glmnet_Model, CPH_Model).
#'
#' @export
RunSurvModels<-function(datatrain, ExpVars, timevar, eventvar, models=c("glmnet","coxph","rulefit","xgboost","gam","gbm","ExpSurvReg","WeibSurvReg","bart","deepsurv"), ntreeRF=300, nvars=20, ...){

  # Ensure event var is 0/1
  datatrain[[eventvar]] <- as.numeric(datatrain[[eventvar]] == 1)

  # Convert character columns to factors (needed by some models like RF)
  datatrainFact<-datatrain
  for (i in ExpVars){ # Iterate through predictors only
    if (is.character(datatrainFact[[i]])){
      datatrainFact[[i]]<-as.factor(datatrainFact[[i]])
    }
  }

  # --- Fit Base Random Forest Models ---
  # Assuming SurvModel_RF is loaded/available
  RF_Model<-tryCatch(
      SurvModel_RF(data=datatrainFact, expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
      error = function(e) {print("Failed fitting RF_Model (all vars)"); warning(e$message); return(NULL)}
  )

  # Select top variables based on importance from the first RF model
  if (!is.null(RF_Model) && !is.null(RF_Model$hd.obj$importance)) {
      importance_scores <- RF_Model$hd.obj$importance
      # Handle different importance structures (matrix or vector)
      if (is.matrix(importance_scores)) {
          # Use first column if matrix (assuming it's the relevant score)
          ExpVars2<-names(sort(importance_scores[,1], decreasing=TRUE))[1:min(length(ExpVars), nvars)]
      } else {
          ExpVars2<-names(sort(importance_scores, decreasing=TRUE))[1:min(length(ExpVars), nvars)]
      }
  } else {
      warning("Could not get importance from RF_Model, using all variables for RF_Model2.")
      ExpVars2 <- ExpVars # Fallback to all variables
  }

  RF_Model2<-tryCatch(
      SurvModel_RF(data=datatrainFact, expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF),
      error = function(e) {print("Failed fitting RF_Model2 (top vars)"); warning(e$message); return(NULL)}
  )


  # --- Fit Additional Requested Models ---
  # Use ExpVars2 (top variables) for these models as in original logic
  ExpVarsForOthers <- ExpVars2

  # Initialize list to store model outputs
  input2<-vector(mode="list")
  input2<-c(input2,list(RF_Model=RF_Model, RF_Model2=RF_Model2)) # Add base RF models first

  # Helper function for fitting with error catching
  fit_model <- function(model_name, fit_func, ...) {
      cat("Fitting", model_name, "...\n")
      model_out <- tryCatch(
          fit_func(...),
          error = function(e) {
              print(paste("Failed fitting", model_name, ":", e$message))
              warning(e$message)
              return(NULL)
          }
      )
      return(model_out)
  }

  if ("bart"%in% models){
    # Assuming SurvModel_BART is loaded/available
    bartout<-fit_model("bart", SurvModel_BART, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, ntree=50) # Reduced ntree for speed?
    input2<-c(input2,list(bart_Model=bartout))
  }

  if ("deepsurv"%in% models){
    # Assuming SurvModel_DeepSurv is loaded/available
    deepsurv_Model<-fit_model("deepsurv", SurvModel_DeepSurv, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, size=5, decay=0.01, maxit=500)
    input2<-c(input2,list(deepsurv_Model=deepsurv_Model))
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
  }
  if ("rulefit"%in% models){
    # Assuming SurvModel_rulefit is loaded/available
    # Rulefit often benefits from more variables, maybe use ExpVars instead of ExpVarsForOthers? Sticking to original logic for now.
    RuleFit_Model<-fit_model("rulefit", SurvModel_rulefit, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar, ntree=ntreeRF, nsample=min(c(500,ceiling(.3*nrow(datatrainFact)))))
    input2<-c(input2,list(RuleFit_Model=RuleFit_Model))
  }
  if ("xgboost"%in% models){
    # Assuming SurvModel_xgboost is loaded/available
    xgboost_Model<-fit_model("xgboost", SurvModel_xgboost, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(xgboost_Model=xgboost_Model))
  }
  if ("gam"%in% models){
     # Assuming SurvModel_GAM is loaded/available
     # Original code had fallback to ExpVars2 if ExpVars failed, simplifying here
    gam_Model<-fit_model("gam", SurvModel_GAM, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(gam_Model=gam_Model))
  }
  if ("gbm"%in% models){
    # Assuming SurvModel_gbm is loaded/available
    gbm_Model<-fit_model("gbm", SurvModel_gbm, data=datatrainFact, expvars=ExpVarsForOthers, timevar=timevar, eventvar=eventvar)
    input2<-c(input2,list(gbm_Model=gbm_Model))
  }

  # --- Return Results ---
  input_params<-list(ExpVars=ExpVars, ExpVars2=ExpVars2, timevar=timevar, eventvar=eventvar)
  # Combine input parameters and fitted models
  return(c(list(input=input_params), input2))
}



#' @title PredictSurvModels
#'
#' @description Utility function to get predictions from several survival models
#' at the same time and calculate an ensemble prediction.
#'
#' @param models List containing the fitted models (output from 'RunSurvModels').
#' @param newdata Data frame with new data for prediction.
#' @param newtimes Numeric vector of times for which predictions are required.
#' @return A list containing:
#'   - ModelPredictions: A list of individual model predictions (interpolated probability matrices).
#'   - NewProbs: Ensembled prediction matrix (averaged on cumulative hazard scale).
#'
#' @export
PredictSurvModels<-function(models, newdata, newtimes){

  # Ensure newtimes includes 0 and is sorted
  newtimes <- sort(unique(c(0, newtimes)))

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
                  survprobMatInterpolator(probsMat=pred_out$Probs, times=pred_out$Times, newtimes=newtimes),
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

  # Predict and interpolate for each model present in the input list
  # Assuming Predict_* functions are loaded/available
  ModelPredictionsList$newprobsrulefit <- predict_and_interp("RuleFit", Predict_SurvModel_rulefit, models$RuleFit_Model)
  ModelPredictionsList$newprobsRF <- predict_and_interp("RF", Predict_SurvModel_RF, models$RF_Model)
  ModelPredictionsList$newprobsRF2 <- predict_and_interp("RF2", Predict_SurvModel_RF, models$RF_Model2)
  # Both glmnet and coxph now use Predict_SurvModel_Cox (unified Cox interface)
  ModelPredictionsList$newprobsglmnet <- predict_and_interp("glmnet", Predict_SurvModel_Cox, models$glmnet_Model)
  ModelPredictionsList$newprobscph <- predict_and_interp("CoxPH", Predict_SurvModel_Cox, models$CPH_Model)
  ModelPredictionsList$newprobsbart_Model <- predict_and_interp("BART", Predict_SurvModel_BART, models$bart_Model)
  ModelPredictionsList$newprobsdeepsurv_Model <- predict_and_interp("DeepSurv", Predict_SurvModel_DeepSurv, models$deepsurv_Model)
  ModelPredictionsList$newprobsgam_Model <- predict_and_interp("GAM", Predict_SurvModel_GAM, models$gam_Model)
  ModelPredictionsList$newprobsgbm_Model <- predict_and_interp("GBM", Predict_SurvModel_gbm, models$gbm_Model)
  ModelPredictionsList$newprobssurvregexp_Model <- predict_and_interp("ExpSurvReg", Predict_SurvModel_SurvReg, models$survregexp_Model)
  ModelPredictionsList$newprobssurvregweib_Model <- predict_and_interp("WeibSurvReg", Predict_SurvModel_SurvReg, models$survregweib_Model)
  ModelPredictionsList$newprobsxgboost <- predict_and_interp("XGBoost", Predict_SurvModel_xgboost, models$xgboost_Model)

  # Filter out NULL predictions (from failed models/predictions/interpolations)
  ValidPredictions <- Filter(Negate(is.null), ModelPredictionsList)

  if (length(ValidPredictions) == 0) {
      warning("No valid predictions obtained from any model.")
      return(list(ModelPredictions=ValidPredictions, NewProbs=NULL))
  }

  # --- Ensemble Predictions ---
  # Average valid predictions on cumulative hazard scale
  # Assuming survprobMatListAveraging is loaded/available
  NewProbs<-tryCatch(
      survprobMatListAveraging(ValidPredictions),
      error = function(e) {
          warning("Ensemble averaging failed: ", e$message)
          return(NULL)
      }
  )

  return(list(ModelPredictions=ValidPredictions, NewProbs=NewProbs))
}
