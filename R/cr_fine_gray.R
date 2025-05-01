#' @title CRModel_FineGray
#'
#' @description Fit a Fine-Gray model for CR outcomes using penalized regression
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2)
#'
#' @return a list object with the following. allvals: AIC values for a range of shrinkage parameters,
#' fit: fit object from 'fastcmprsk::fastCrrp',
#' CIF.hat: estimated CIFs in training data,
#' expvars: character vector of explanatory variables,
#' timesfit: times for which CIFs are given,
#' varprof: profile of explanatory variables,
#' meanTrain: means of training features,
#' sdTrain: standard deviations of training features,
#' loadings: SVD loadings.
#'
#' @importFrom fastcmprsk fastCrrp Crisk
#' @importFrom stats AIC model.matrix quantile sd
#' @export
CRModel_FineGray<-function(data,expvars, timevar, eventvar){

  # Assuming VariableProfile is defined elsewhere or loaded via namespace
  # varprof<-VariableProfile(data, expvars)
  varprof <- list() # Placeholder

  timevarvec<-data[, timevar]
  eventvarvec<-as.integer(data[, eventvar])

  covmat<-stats::model.matrix(~-1+., data=data[, expvars])

  lambdavec<-seq(0.005,.99,length=20)
  alphavec<-seq(0.05,.95,length=20)

  allvals<-as.data.frame(expand.grid(lambdavec,alphavec))
  colnames(allvals)<-c("lambda","alpha")
  allvals$AIC<-apply(allvals,1,function(x){
    aic<- tryCatch(stats::AIC(fastcmprsk::fastCrrp(fastcmprsk::Crisk(timevarvec,eventvarvec)~covmat,lambda = x[1], alpha=x[2], penalty = "ENET")), error=function(e){Inf})
    return(aic)
  }
  )
  covmat<-scale(covmat, center=TRUE, scale=TRUE)
  meanTrain<-attr(covmat, "scaled:center")
  sdTrain<-attr(covmat, "scaled:scale")

  svdcovmat<-svd(covmat) # Changed from stats::svd
  Feat<-(covmat%*%svdcovmat$v)[,1:min(c(20,ncol(covmat)))]
  fit<-fastcmprsk::fastCrrp(fastcmprsk::Crisk(timevarvec,eventvarvec)~Feat,
                            lambda = allvals[which.min(allvals[,3]),1],
                            alpha=allvals[which.min(allvals[,3]),2],
                            penalty = "ENET",
                            standardize=FALSE,
                            max.iter = 5000)
  CIF.hat <- tcrossprod(exp(cbind(Feat)%*%c(unlist(fit$coef))),  fit$breslowJump[, 2]) #This is baseline cumulative hazard
  CIF.hat<-rbind(0,apply(CIF.hat, 1, cumsum))
  CIF.hat <- 1 - exp(-CIF.hat)
  return(list(allvals=allvals,fit=fit,CIF.hat=CIF.hat, expvars=expvars, timesfit= c(0,fit$breslowJump[,1]), varprof=varprof, meanTrain=meanTrain,sdTrain=sdTrain, loadings=svdcovmat$v))
}


#' @title Predict_CRModel_FG
#'
#' @description Get predictions from a Fine-Gray CR model for a test dataset
#'
#' @param modelout the output from 'CRModel_FineGray'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_FG<-function(modelout, newdata){
  covmat<-stats::model.matrix(~-1+., data=newdata[, modelout$expvars])
  Feat<-(scale(covmat, scale=modelout$sdTrain, center=modelout$meanTrain)%*%modelout$loadings)[,1:min(c(20,ncol(covmat)))]
  CIF.hat <- tcrossprod(exp(Feat%*%c(unlist(modelout$fit$coef))),  modelout$fit$breslowJump[, 2]) #This is baseline cumulative hazard
  CIF.hat<-apply(CIF.hat, 1, cumsum)
  CIF.hat <- rbind(0,1 - exp(-CIF.hat))
  return(list(CIFs=cbind(t(CIF.hat)),Times=c(modelout$timesfit)))
}
