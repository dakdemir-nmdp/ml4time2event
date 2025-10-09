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
#'  eventvar: event variable. The second is also a list that contain the individual model outputs.
#' @export
RunCRModels<-function(datatrain, ExpVars, timevar, eventvar, models=c("FG", "rulefit", "bart", "cox"), ntreeRF=300, varsel=FALSE){

  datatrainFact<-datatrain

  for (i in 1:ncol(datatrainFact)){
    if (is.character(datatrainFact[,i])){
      datatrainFact[,i]<-as.factor(datatrainFact[,i])
    }
  }
  # Assuming CRModel_RF is loaded/available
  RF_Model<-CRModel_RF(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)
  # Get variable importance - use the first event's importance
  imp_scores <- RF_Model$rf_model$importance[,1]
  # Ensure we have at least some variables selected
  if (all(is.na(imp_scores)) || all(imp_scores <= 0, na.rm = TRUE)) {
    # If all importances are NA or <= 0, use all variables
    ExpVars2 <- ExpVars
  } else {
    # Sort by importance and take top variables (at least 1, at most half)
    n_vars <- min(length(ExpVars), max(1, length(ExpVars) %/% 2))
    ExpVars2 <- names(sort(imp_scores, decreasing = TRUE, na.last = TRUE))[seq_len(n_vars)]
  }
  RF_Model2<-CRModel_RF(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)

 if (!varsel){
  if ("FG" %in% models){
    # Assuming CRModel_FineGray is loaded/available
    FG_Model<-tryCatch(CRModel_FineGray(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting FG")
      return()
    })
  }
  if ("bart" %in% models){
    # Assuming CRModel_BART is loaded/available
    BART_Model<-tryCatch(CRModel_BART(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, ntree=ntreeRF), error=function(e){
      print("Failed fitting BART")
      return()
    })
  }
  if ("cox" %in% models){
    # Assuming CRModel_Cox is loaded/available
    Cox_Model<-tryCatch(CRModel_Cox(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting COX")
      return()
    })
  }

  if ("rulefit" %in% models){
    # Assuming CRModel_rulefit is loaded/available
    rulefit_Model<-tryCatch(CRModel_rulefit(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
       print("Failed fitting Rulefit")
      return()
    })
  }
  if ("xgboost" %in% models){
    # Assuming CRModel_xgboost is loaded/available
  xgboost_Model<-tryCatch(CRModel_xgboost(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, event_codes=1, nrounds=100), error=function(e){
       print("Failed fitting XGBoost")
      return()
    })
  }
  if ("gam" %in% models){
    # Assuming CRModel_GAM is loaded/available
    gam_Model<-tryCatch(CRModel_GAM(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, failcode=1), error=function(e){
       print("Failed fitting GAM")
      return()
    })
  }
  if ("survreg" %in% models){
    # Assuming CRModel_SurvReg is loaded/available
    survreg_Model<-tryCatch(CRModel_SurvReg(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, failcode=1, dist="exponential"), error=function(e){
       print("Failed fitting SurvReg")
      return()
    })
  }
 } else {
   if ("FG" %in% models){
     FG_Model<-tryCatch(CRModel_FineGray(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       print("Failed fitting FG")
       return()
     })
   }
   if ("bart" %in% models){
     BART_Model<-tryCatch(CRModel_BART(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, ntree=ntreeRF), error=function(e){
       print("Failed fitting BART")
       return()
     })
   }
   if ("cox" %in% models){
     Cox_Model<-tryCatch(CRModel_Cox(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       print("Failed fitting COX")
       return()
     })
   }

   if ("rulefit" %in% models){
     rulefit_Model<-tryCatch(CRModel_rulefit(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
       print("Failed fitting Rulefit")
       return()
     })
   }
   if ("xgboost" %in% models){
  xgboost_Model<-tryCatch(CRModel_xgboost(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, event_codes=1, nrounds=100), error=function(e){
       print("Failed fitting XGBoost")
       return()
     })
   }
   if ("gam" %in% models){
     gam_Model<-tryCatch(CRModel_GAM(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, failcode=1), error=function(e){
       print("Failed fitting GAM")
       return()
     })
   }
   if ("survreg" %in% models){
     survreg_Model<-tryCatch(CRModel_SurvReg(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, failcode=1, dist="exponential"), error=function(e){
       print("Failed fitting SurvReg")
       return()
     })
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

  return(c(list(input=input),input2))
}






#' @title  PredictCRModels
#' @description  Utility function to get predictions from  several competing risks models
#' at the same time. Also adds an ensemble prediction for the probabilities.
#' @param models the output from 'RunCRModels'
#' @param newdata the data for which the predictions are to be calculated
#' @param newtimes the times for which the predictions obtained
#' @return a list containing the following objects.
#' ModelPredictions: a list of individual model predictions,
#' NewProbs: ensembled predictions for CIFs
#'
#' @export
PredictCRModels<-function(models, newdata, newtimes){
  # for (i in 1:ncol(newdata)){
  #   if (is.factor(newdata[,i])){
  #     newdata[,i]<-as.character(newdata[,i])
  #   }
  # }

  if (!is.null(models$RF_Model)){
    # Assuming Predict_CRModel_RF is loaded/available
    Predict_RF<-Predict_CRModel_RF(models$RF_Model, newdata=newdata)
  }
  if (!is.null(models$RF_Model2)){
    Predict_RF2<-Predict_CRModel_RF(models$RF_Model2, newdata=newdata)
  }
  if (!is.null(models$FG_Model)){
    # Assuming Predict_CRModel_FineGray is loaded/available
    Predict_FG<-Predict_CRModel_FineGray(models$FG_Model, newdata=newdata)
  }
  if (!is.null(models$BART_Model)){
    # Assuming Predict_CRModel_BART is loaded/available
    Predict_BART<-Predict_CRModel_BART(models$BART_Model, newdata=newdata)
  }

  if (!is.null(models$Cox_Model)){
    # Assuming Predict_CRModel_Cox is loaded/available
    Predict_Cox<-Predict_CRModel_Cox(models$Cox_Model, newdata=newdata)
  }
  if (!is.null(models$rulefit_Model)){
    # Assuming Predict_CRModel_rulefit is loaded/available
    Predict_rulefit<-Predict_CRModel_rulefit(models$rulefit_Model, newdata=newdata)
  }
  if (!is.null(models$xgboost_Model)){
    # Assuming Predict_CRModel_xgboost is loaded/available
    Predict_xgboost<-Predict_CRModel_xgboost(models$xgboost_Model, newdata=newdata)
  }
  if (!is.null(models$gam_Model)){
    # Assuming Predict_CRModel_GAM is loaded/available
    Predict_gam<-Predict_CRModel_GAM(models$gam_Model, newdata=newdata)
  }
  if (!is.null(models$survreg_Model)){
    # Assuming Predict_CRModel_SurvReg is loaded/available
    Predict_survreg<-Predict_CRModel_SurvReg(models$survreg_Model, newdata=newdata)
  }

  # Assuming cifMatInterpolaltor is loaded/available
  if (!is.null(models$RF_Model)){
    newprobsRF<-cifMatInterpolaltor(probsMat=t(Predict_RF$CIFs),times=Predict_RF$Times, newtimes=newtimes)
  }
  if (!is.null(models$RF_Model2)){
    newprobsRF2<-cifMatInterpolaltor(probsMat=t(Predict_RF2$CIFs),times=Predict_RF2$Times, newtimes=newtimes)
  }
  if (!is.null(models$FG_Model)){
    newprobsFG<-cifMatInterpolaltor(probsMat=t(Predict_FG$CIFs),times=Predict_FG$Times, newtimes=newtimes)
  }
  if (!is.null(models$BART_Model)){
    newprobsBART<-cifMatInterpolaltor(probsMat=t(Predict_BART$CIFs),times=Predict_BART$Times, newtimes=newtimes)
  }
  if (!is.null(models$Cox_Model)){
    newprobsCox<-cifMatInterpolaltor(probsMat=t(Predict_Cox$CIFs),times=Predict_Cox$Times, newtimes=newtimes)
  }
  if (!is.null(models$rulefit_Model)){
    newprobsrulefit<-cifMatInterpolaltor(probsMat=t(Predict_rulefit$CIFs),times=Predict_rulefit$Times, newtimes=newtimes)
  }
  if (!is.null(models$xgboost_Model)){
    newprobsxgboost<-cifMatInterpolaltor(probsMat=t(Predict_xgboost$CIFs),times=Predict_xgboost$Times, newtimes=newtimes)
  }
  if (!is.null(models$gam_Model)){
    newprobsgam<-cifMatInterpolaltor(probsMat=t(Predict_gam$CIFs),times=Predict_gam$Times, newtimes=newtimes)
  }
  if (!is.null(models$survreg_Model)){
    newprobsurvreg<-cifMatInterpolaltor(probsMat=t(Predict_survreg$CIFs),times=Predict_survreg$Times, newtimes=newtimes)
  }

  ModelPredictions<-list()
  if (!is.null(models$RF_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsRF=newprobsRF))
  }
  if (!is.null(models$RF_Model2)){
    ModelPredictions<-c(ModelPredictions,list(newprobsRF2=newprobsRF2))
  }
  if (!is.null(models$FG_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsFG=newprobsFG))
  }
  if (!is.null(models$BART_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsBART=newprobsBART))
  }
  if (!is.null(models$Cox_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsCox=newprobsCox))
  }
  if (!is.null(models$rulefit_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsrulefit=newprobsrulefit))
  }
  if (!is.null(models$xgboost_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsxgboost=newprobsxgboost))
  }
  if (!is.null(models$gam_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsgam=newprobsgam))
  }
  if (!is.null(models$survreg_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsurvreg=newprobsurvreg))
  }

  # Assuming cifMatListAveraging is loaded/available
  NewProbs<-cifMatListAveraging(ModelPredictions, type = "prob")
  NewProbs<-rbind(NewProbs) # Ensure it's a matrix
  return(list(ModelPredictions=ModelPredictions,NewProbs=NewProbs))
}
