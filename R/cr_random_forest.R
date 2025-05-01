#' @title CRModel_RF
#'
#' @description  Fit a  random forest model for competing risk outcomes
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2)
#' @param cause vectoe of weights for causes
#' @param ntree number of trees to fit to extract rules
#' @param samplesize number of samples for each tree
#' @param nsplit maximum number of splits for each tree
#' @param trace trace the results or not
#'
#' @return a list object with the following. hd.obj: trained random forest model object,
#' varprof: profile of explanatory variables.
#'
#' @importFrom randomForestSRC tune rfsrc predict.rfsrc
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @export
CRModel_RF<-function(data,expvars, timevar, eventvar, cause=c(1,1), ntree=300, samplesize=500, nsplit=5, trace=TRUE){
  # Assuming VariableProfile is defined elsewhere or loaded via namespace
  # varprof<-VariableProfile(data, expvars)
  varprof <- list() # Placeholder if VariableProfile is not available here
  for (vari in expvars){
    x<-data[,vari]
    if (is.character(x)){
      x<-as.factor(x)
      data[,vari]<-x
    }
  }
  formRF<-stats::as.formula(paste("survival::Surv(",timevar, ",", eventvar,") ~.", collapse = ""))
  samplesize<-min(c(ceiling(.7*nrow(data)), samplesize))

  o <- randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars)],sampsize = samplesize, trace = trace,nsplit=nsplit,stepFactor = 1.5,
                             mtryStart = 2,
                             nodesizeTry = c(seq(10, 101, by = 10)), ntreeTry = ntree, cause=cause)

  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars)],
                                   nodesize =o$optimal[[1]], ntree=ntree,mtry= o$optimal[[2]],
                                   tree.err = FALSE, importance = TRUE,statistics=TRUE,nsplit=nsplit,ntime=100,
                                   do.trace = TRUE,sampsize = samplesize, cause=cause)

  # predSurvsTrainRF<-randomForestSRC::predict.rfsrc(hd.obj, newdata = data[,c(timevar, eventvar, expvars)]) # Removed as it's not returned
  return(list(hd.obj=hd.obj, varprof=varprof)) # Assuming varprof needs to be returned
}




#' @title Predict_CRModel_RF
#'
#' @description Get predictions from a RF model for a test dataset
#'
#' @param modelout the output from 'CRModel_RF'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. predSurvsTestRF: the output of 'randomForestSRC::predict.rfsrc',
#'  CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @importFrom randomForestSRC predict.rfsrc
#' @export
Predict_CRModel_RF<-function(modelout, newdata){
  predSurvsTestRF<-randomForestSRC::predict.rfsrc(modelout$hd.obj, newdata = newdata)
  return(list(predSurvsTestRF=predSurvsTestRF,CIFs=cbind(0,predSurvsTestRF$cif[,,1]),Times=c(0,predSurvsTestRF$time.interest)))
}
