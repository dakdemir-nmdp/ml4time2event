#' @title SurvModel_Cox
#'
#' @description Fit a Cox model for survival outcomes using forward selection via stats::step.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#'
#' @return a list of three items: cph_model: the final Cox PH model object after step selection,
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables.
#'
#' @importFrom survival coxph Surv survfit
#' @importFrom stats step as.formula
#' @export
SurvModel_Cox<-function(data,expvars, timevar, eventvar){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Define initial formula including all predictors
  form<-stats::as.formula(paste(paste("survival::Surv(",timevar, ",", eventvar,") ~1+", collapse = ""),paste(expvars,collapse="+"),sep=""))
  print(form) # Print the initial formula

  # Prepare data subset
  XYTrain<-data[,c(timevar,eventvar,expvars), drop=FALSE]

  # Fit initial Cox model
  cph_model <- survival::coxph(form, data = XYTrain, x=TRUE, y=TRUE) # Keep x=T, y=T for step

  # Perform forward selection using step (AIC-based)
  # Note: step by default performs backward elimination from the full model.
  # To perform forward selection, one would typically start with a null model
  # and specify scope. However, the original code calls step on the full model,
  # implying backward elimination was intended. We keep this behavior.
  stepMod<-stats::step(cph_model, direction = "backward", trace = 0) # Use direction="backward", suppress trace

  # Get unique event times from training data
  times<-sort(unique(XYTrain[XYTrain[[eventvar]] == 1, timevar])) # Get times where event occurred

  # cph_modelTrainPredict<- survival::survfit(stepMod, newdata=XYTrain,times=times) # Not returned

  return(list(cph_model=stepMod, times=times, varprof=varprof))
}

#' @title Predict_SurvModel_Cox
#'
#' @description Get predictions from a Cox survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_Cox' (list containing 'cph_model' and 'times')
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the times at which the probabilities are calculated (including 0),
#' cph_modelTestPredict: the raw output object from 'survival::survfit'.
#'
#' @importFrom survival survfit
#' @export
Predict_SurvModel_Cox<-function(modelout, newdata){
  # Predict survival curves using the selected model and training times
  cph_modelTestPredict<- survival::survfit(modelout$cph_model, newdata=newdata, times=modelout$times)

  # Extract probabilities and times
  Probs<-cph_modelTestPredict$surv # Matrix: rows=times, cols=observations
  Times<-cph_modelTestPredict$time

  # Add time 0 with probability 1 if missing
  if (sum(Times==0)==0){
    Times<-c(0,Times)
    Probs<-rbind(rep(1, ncol(Probs)),Probs)
  }

  return(list(Probs=Probs,
              Times=Times,
              cph_modelTestPredict=cph_modelTestPredict))
}
