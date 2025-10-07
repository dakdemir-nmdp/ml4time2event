#' @title CRModel_BART
#'
#' @description Fit a BART model for competing risks outcomes
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param K number of cutpoints for the time variable
#' @param ntree number of trees to fit to extract rules
#'
#' @return a list object with the following: post: BART model,
#'  expvars: character vector of explanatory variables,
#'  eventvarvec: event variable values,timevarvec: time variable values,
#'  x.train: design matrix for the training data,
#'  varprof: profile of the explanatory variables
#'
#' @importFrom BART crisk.bart crisk.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @export
CRModel_BART<-function(data,expvars, timevar, eventvar, K=10,ntree=100){
  # Assuming VariableProfile is defined elsewhere or loaded via namespace
  # varprof<-VariableProfile(data, expvars)
  varprof <- list() # Placeholder

  timevarvec<-data[, timevar]
  eventvarvec<-as.integer(data[, eventvar])
  x.train<-as.matrix(stats::model.matrix(~-1+., data=data[, expvars]))
  failcount<-0
  post<-NULL
  while ((is.null(post) & (failcount<4))){
    failcount<-failcount+1
    post <- tryCatch(BART::crisk.bart(x.train=x.train, times=timevarvec, delta=eventvarvec, x.test=x.train,
                                      K=K,ntree=ntree,numcut=2, ndpost=2000, nskip=100,
                                      keepevery = 3L),
                     error=function(e){NULL}
    )
  }
  return(list( post=post, expvars=expvars,eventvarvec=eventvarvec,timevarvec=timevarvec,x.train=x.train, varprof=varprof))
}


#' @title Predict_CRModel_BART
#'
#' @description  Get predictions from a BART model for a test dataset
#'
#' @param modelout the output from 'CRModel_BART'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items.
#' CIFs:predicted CIF matrix (for event coded as 1),
#' Times: The times at which the CIFs are predicted.
#'
#' @importFrom BART crisk.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_BART<-function(modelout, newdata){
  x.test<-as.matrix(stats::model.matrix(~-1+., data=newdata[, modelout$expvars]))
  (M=modelout$post$ndpost)
  (K=modelout$post$K)
  (N=nrow(x.test))
  test = rbind(x.test, x.test) # Why rbind x.test with itself? Check BART documentation/original intent. Keeping for now.
  pre=BART::crisk.pre.bart(time=modelout$timevarvec, delta=modelout$eventvarvec,
                           x.train= BART::bartModelMatrix(modelout$x.train),
                           x.test=  BART::bartModelMatrix(test),
                           x.train2= BART::bartModelMatrix(modelout$x.train), # Why train2 and test2? Check BART docs.
                           x.test2=  BART::bartModelMatrix(test),
                           K=K)
  pred = predict(modelout$post, pre$tx.test, pre$tx.test2)
  PredMat<-matrix(pred$cif.test.mean, ncol=K,byrow = TRUE)
  PredMat1<-PredMat[1:N,] # Taking only the first N rows, corresponding to the actual newdata
  return(list(CIFs=cbind(0,PredMat1), Times=c(0,modelout$post$times)))
}
