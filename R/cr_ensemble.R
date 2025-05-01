#' @title RunCRModels
#' @description  utility function to run several CR models at the same time
#' @param datatrain data frame with explanatory and outcome variables
#' @param ExpVars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded as 0,1,2),
#' 1 is the event of interest
#' @param models  a vector of models to be fitted to be chosen among
#' "FG", "rulefit", "bart", "cox". Two random forest models are
#' automatically fitted and dont need to be listed here
#' @param ntreeRF number of trees for Random Forest models
#' @param varsel  logical indicating whether variable selection to be
#' applied before fitting models "FG", "rulefit", "bart", "cox"
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
  ExpVars2<-names(sort(RF_Model$hd.obj$importance[,1], decreasing=TRUE))[1:min(length(ExpVars),20)]
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
    # Assuming Predict_CRModel_FG is loaded/available
    Predict_FG<-Predict_CRModel_FG(models$FG_Model, newdata=newdata)
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

  # Assuming cifMatInterpolaltor is loaded/available
  if (!is.null(models$RF_Model)){
    newprobsRF<-cifMatInterpolaltor(probsMat=Predict_RF$CIFs,times=Predict_RF$Times, newtimes=newtimes)
  }
  if (!is.null(models$RF_Model2)){
    newprobsRF2<-cifMatInterpolaltor(probsMat=Predict_RF2$CIFs,times=Predict_RF2$Times, newtimes=newtimes)
  }
  if (!is.null(models$FG_Model)){
    newprobsFG<-cifMatInterpolaltor(probsMat=Predict_FG$CIFs,times=Predict_FG$Times, newtimes=newtimes)
  }
  if (!is.null(models$BART_Model)){
    newprobsBART<-cifMatInterpolaltor(probsMat=Predict_BART$CIFs,times=Predict_BART$Times, newtimes=newtimes)
  }
  if (!is.null(models$Cox_Model)){
    newprobsCox<-cifMatInterpolaltor(probsMat=Predict_Cox$CIFs,times=Predict_Cox$Times, newtimes=newtimes)
  }
  if (!is.null(models$rulefit_Model)){
    newprobsrulefit<-cifMatInterpolaltor(probsMat=Predict_rulefit$CIFs,times=Predict_rulefit$Times, newtimes=newtimes)
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

  # Assuming cifMatListAveraging is loaded/available
  NewProbs<-cifMatListAveraging(ModelPredictions)
  NewProbs<-rbind(NewProbs) # Ensure it's a matrix
  return(list(ModelPredictions=ModelPredictions,NewProbs=NewProbs))
}
