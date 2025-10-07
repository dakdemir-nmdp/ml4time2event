#' @title CRModel_Cox
#'
#' @description Fit a Cox model for CR outcomes using forward selection
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2)
#'
#' @return a list object with the following. cph_model: trained Cox PH model object,
#' times: unique times in the training data,
#' varprof: profile of explanatory variables.
#'
#' @importFrom survival coxph Surv survfit
#' @importFrom stats step as.formula
#' @export
CRModel_Cox<-function(data,expvars, timevar, eventvar){
  # Assuming VariableProfile is defined elsewhere or loaded via namespace
  # varprof<-VariableProfile(data, expvars)
  varprof <- list() # Placeholder
  form<-stats::as.formula(paste(paste("survival::Surv(",timevar, ",", eventvar,") ~1+", collapse = ""),paste(expvars,collapse="+"),sep=""))
  data[,eventvar]<-as.factor(as.numeric(data[,eventvar]))
  XYTrain<-data[,c(timevar,eventvar,expvars)]
  XYTrain$id<-as.character(1:nrow(XYTrain)) # Add 'id' column back
  cph_model <- survival::coxph(form, data = XYTrain,x=TRUE,y=TRUE, id=id) # Use id=id
  stepMod<-stats::step(cph_model, data = XYTrain) # Explicitly pass data to step
  # cph_modelTrainPredict<- survival::survfit(stepMod, newdata=XYTrain,times=sort(unique(XYTrain[,timevar]))) # Removed as not returned
  return(list(cph_model=stepMod,times=sort(unique(XYTrain[,timevar])),  varprof=varprof))
}

#' @title Predict_CRModel_Cox
#'
#' @description Get predictions from a Cox CR model for a test dataset
#'
#' @param modelout the output from 'CRModel_Cox'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. cph_modelTestPredict: the output of 'survival::survfit',
#'  CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @importFrom survival survfit
#' @export
Predict_CRModel_Cox<-function(modelout, newdata){
  cph_modelTestPredict<- survival::survfit(modelout$cph_model, newdata=newdata,times=modelout$times)
  return(list(cph_modelTestPredict=cph_modelTestPredict, CIFs=cbind(0,t(cph_modelTestPredict$pstate[,,2])),Times=c(0,cph_modelTestPredict$time)))
}
