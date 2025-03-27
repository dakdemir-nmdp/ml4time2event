## Module: CR Models Functions
## This file contains functions for competing risk models.
## Sections include:
##   - Random Forest Models (CRModel_RF, Predict_CRModel_RF)
##   - Cox Models (CRModel_Cox, Predict_CRModel_Cox)
##   - Fine-Gray Model (CRModel_FineGray, Predict_CRModel_FG)
##   - Rulefit Models (CRModel_rulefit, Predict_CRModel_rulefit)
##   - BART Models (CRModel_BART, Predict_CRModel_BART)
##   - Utility functions for ensemble prediction and interpolation.
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
#'
#'
#' @export
CRModel_RF<-function(data,expvars, timevar, eventvar, cause=c(1,1), ntree=300, samplesize=500, nsplit=5, trace=TRUE){
  varprof<-VariableProfile(data, expvars)
  for (vari in expvars){
    x<-data[,vari]
    if (is.character(x)){
      x<-as.factor(x)
      data[,vari]<-x
    }
  }
  formRF<-as.formula(paste("Surv(",timevar, ",", eventvar,") ~.", collapse = ""))
  samplesize<-min(c(ceiling(.7*nrow(data)), samplesize))

  o <- randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars)],sampsize = samplesize, trace = trace,nsplit=nsplit,stepFactor = 1.5,
                             mtryStart = 2,
                             nodesizeTry = c(seq(10, 101, by = 10)), ntreeTry = ntree, cause=cause)

  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars)],
                                   nodesize =o$optimal[[1]], ntree=ntree,mtry= o$optimal[[2]],
                                   tree.err = FALSE, importance = TRUE,statistics=TRUE,nsplit=nsplit,ntime=100,
                                   do.trace = TRUE,sampsize = samplesize, cause=cause)

  predSurvsTrainRF<-randomForestSRC::predict.rfsrc(hd.obj, newdata = data[,c(timevar, eventvar, expvars)])
  return(list(hd.obj=hd.obj, varprof=varprof))
}




#' @title Predict_CRModel_RF
#'
#' @description Get predictions from a RF model for a test dataset
#'
#' @param modelout the output from 'SurvModel_RF'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. predSurvsTestRF: the output of 'randomForestSRC::predict.rfsrc',
#'  CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @export
Predict_CRModel_RF<-function(modelout, newdata){
  predSurvsTestRF<-randomForestSRC::predict.rfsrc(modelout$hd.obj, newdata = newdata)
  return(list(predSurvsTestRF=predSurvsTestRF,CIFs=cbind(0,predSurvsTestRF$cif[,,1]),Times=c(0,predSurvsTestRF$time.interest)))
}



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
#' varprof: profile of explanatory variables.
#'
#' @export
CRModel_Cox<-function(data,expvars, timevar, eventvar){
  varprof<-VariableProfile(data, expvars)
  form<-as.formula(paste(paste("survival::Surv(",timevar, ",", eventvar,") ~1+", collapse = ""),paste(expvars,collapse="+"),sep=""))
  data[,eventvar]<-as.factor(as.numeric(data[,eventvar]))
  XYTrain<-data[,c(timevar,eventvar,expvars)]
  XYTrain$idvar<-as.character(1:nrow(XYTrain))
  cph_model <- survival::coxph(form, data = XYTrain,x=TRUE,y=TRUE, id=idvar)
  stepMod<-stats::step(cph_model)
  cph_modelTrainPredict<- survival::survfit(stepMod, newdata=XYTrain,times=sort(unique(XYTrain[,timevar])))
  return(list(cph_model=stepMod,times=sort(unique(XYTrain[,timevar])),  varprof=varprof))
}

#' @title Predict_CRModel_Cox
#'
#' @description Get predictions from a Cox CR model for a test dataset
#'
#' @param modelout the output from 'SurvModel_RF'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. cph_modelTestPredict: the output of 'survival::survfit',
#'  CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @export
Predict_CRModel_Cox<-function(modelout, newdata){
  cph_modelTestPredict<- survival::survfit(modelout$cph_model, newdata=newdata,times=modelout$times)
  return(list(cph_modelTestPredict=cph_modelTestPredict, CIFs=cbind(0,t(cph_modelTestPredict$pstate[,,2])),Times=c(0,cph_modelTestPredict$time)))
}











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
#' varprof: profile of explanatory variables.
#'
#' @export
CRModel_FineGray<-function(data,expvars, timevar, eventvar){

  varprof<-VariableProfile(data, expvars)

  timevarvec<-data[, timevar]
  eventvarvec<-as.integer(data[, eventvar])

  covmat<-model.matrix(~-1+., data=data[, expvars])

  lambdavec<-seq(0.005,.99,length=20)
  alphavec<-seq(0.05,.95,length=20)

  allvals<-as.data.frame(expand.grid(lambdavec,alphavec))
  colnames(allvals)<-c("lambda","alpha")
  allvals$AIC<-apply(allvals,1,function(x){
    aic<- tryCatch(AIC(fastcmprsk::fastCrrp(fastcmprsk::Crisk(timevarvec,eventvarvec)~covmat,lambda = x[1], alpha=x[2], penalty = "ENET")), error=function(e){Inf})
    return(aic)
  }
  )
  covmat<-scale(covmat, center=TRUE, scale=TRUE)
  meanTrain<-attr(covmat, "scaled:center")
  sdTrain<-attr(covmat, "scaled:scale")

  svdcovmat<-svd(covmat)
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
#' @description Get predictions from a Cox CR model for a test dataset
#'
#' @param modelout the output from 'SurvModel_RF'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. CIFs=predicted CIF matrix (for event coded as 1),
#'  Times: The times at which the CIFs are predicted.
#'
#' @export
Predict_CRModel_FG<-function(modelout, newdata){
  covmat<-model.matrix(~-1+., data=newdata[, modelout$expvars])
  Feat<-(scale(covmat, scale=modelout$sdTrain, center=modelout$meanTrain)%*%modelout$loadings)[,1:min(c(20,ncol(covmat)))]
  CIF.hat <- tcrossprod(exp(Feat%*%c(unlist(modelout$fit$coef))),  modelout$fit$breslowJump[, 2]) #This is baseline cumulative hazard
  CIF.hat<-apply(CIF.hat, 1, cumsum)
  CIF.hat <- rbind(0,1 - exp(-CIF.hat))
  return(list(CIFs=cbind(t(CIF.hat)),Times=c(modelout$timesfit)))
}


#' @title CRModel_rulefit
#'
#' @description Fit a rulefit model for competing risks outcomes
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param ntree number of trees to fit to extract rules
#' @param nsample number of samples for each tree
#' @param keepvars these variables will be used in each bagging iteration
#' @param cuttimes these cut times are used for calculating pseudo observation
#'
#' @return a list object with the following: datasamp: a sample of data,
#' ruleslist: list of the fitted rules,
#' ctreelist: list of tree models,
#' yTrain: survival object for outcome variables,
#' CRrulefitModel: fitted Fine-Gray model object,
#' expvars: explanatory variables used in the model
#' dupcols: Duplicated columns in training data (duplicated rules)
#' varprof: profile of explanatory variables.
#'
#'
#' @export
CRModel_rulefit <- function(data,
                            expvars,
                            timevar,
                            eventvar,
                            ntree = 300,
                            nsample = 300,
                            keepvars = NULL,
                            cuttimes =NULL
) {
  if (is.null(cuttimes)){

    x<-data[,timevar]

    cuttimes<-quantile(x, c(.1,.25,.50,.70,.8))
    cuttimesReg<-quantile(x, c(.25,.50,.60, .70, .80))

  }
  varprof<-VariableProfile(data, expvars)

  formClass <- as.formula(paste("ClassVar ~.", collapse = ""))
  formReg <- as.formula(paste("RegVar ~.", collapse = ""))
  formSurv <- as.formula(paste("Surv(",timevar,",", eventvar,"==1) ~.", collapse = ""))


  ctreelist <- lapply(1:ntree, function(repi){
    sampcols <-
      union(keepvars, sample(expvars, min(length(expvars), sample(c(
        1:10
      ), 1))))
    selmodel <- sample(c(2,3), 1)
   selmodel<-2
    if (selmodel == 1) {

      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(1:nrow(data), nsample, replace = T)
      datasampl <- data[samprows, colnames(data) %in% usevars]

      datasampl$ClassVar <-
        as.character(datasampl[, colnames(datasampl) %in% timevar] < sample(cuttimes, 1) &
                       datasampl[, colnames(datasampl) %in% eventvar] == 1)
      datasampl <-
        datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
      rpcontrol <-
        rpart::rpart.control(
          minsplit = rpois(1,1)+1,
          minbucket = rpois(1,30)+1,
          cp = 0.01 * runif(1),
          maxcompete = rpois(1,30)+1,
          maxsurrogate = rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formClass,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else if (selmodel == 2) {
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(1:nrow(data), nsample, replace = T)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      pout<-pseudo::pseudoci(time=datasampl[,timevar],event=datasampl[,eventvar],tmax=sample(cuttimesReg,1))
      datasampl$RegVar <-c(unlist(pout$pseudo$cause1))
      datasampl$RegVar<-(datasampl$RegVar-mean(datasampl$RegVar))/sd(datasampl$RegVar)
      datasampl <-
        datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
      rpcontrol <-
        rpart::rpart.control(
          minsplit = rpois(1,1)+1,
          minbucket = rpois(1,30)+1,
          cp = 0.01 * runif(1),
          maxcompete = rpois(1,30)+1,
          maxsurrogate = rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formReg,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else if (selmodel == 3) {
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(1:nrow(data), nsample, replace = T)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      rpcontrol <-
        rpart::rpart.control(
          minsplit = rpois(1,1)+1,
          minbucket = rpois(1,30)+1,
          cp = 0.01 * runif(1),
          maxcompete = rpois(1,30)+1,
          maxsurrogate = rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formSurv,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    }
    return(rpartmodel)
  })
  ruleslist <-
    lapply(ctreelist, function(x) {
      listrules(partykit::as.party(x))
    })
  ruleslist <-
    lapply(ctreelist, function(x) {
      listrules(partykit::as.party(x))
    })
  RulesTrain <-
    lapply(ctreelist, function(x)
      predict(partykit::as.party(x), newdata = data, type = "node"))
  RulesTrainMatList <- lapply(1:length(RulesTrain), function(x) {
    Train <- RulesTrain[[x]]
    rules <- ruleslist[[x]]
    Train <- factor(as.character(Train), levels = names(rules))
    if (length(unique(names(rules))) > 1) {
      MM <- model.matrix(~ -1 + Train)
      colnames(MM) <-
        paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
      MM
    } else {
      NULL
    }
  })
  TrainMat <- Reduce("cbind", RulesTrainMatList)
  TrainMat<-cbind(model.matrix(~.-1, data=data[,expvars]), TrainMat)
  usecols<-colnames(TrainMat)
  nrowforsample<-min(c(nrow(data),500))
  CRrulefitModel<-CRModel_FineGray(data=cbind(TrainMat, data[, colnames(data)%in%c(eventvar,timevar)]),expvars=usecols, timevar=timevar , eventvar=eventvar)
  return(
    list(
      datasamp=data[1:nrowforsample,expvars],
      ruleslist = ruleslist,
      ctreelist = ctreelist,
      CRrulefitModel=CRrulefitModel,
      expvars=expvars,
      usecols=usecols,
      varprof=varprof
    )
  )
}


#' @title Predict_CRModel_rulefit
#'
#' @description Get predictions from a rulefit model for a test dataset
#'
#' @param modelout the output from 'CRModel_rulefit'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items. TestMat: model matrix for the test data,
#' PredRFOut: prediction object from Fine-Gray model,
#' CIFs:predicted CIF matrix (for event coded as 1),
#' Times: The times at which the CIFs are predicted.
#'
#' @export
Predict_CRModel_rulefit <- function(modelout, newdata) {
  RulesTest <-
    lapply(modelout$ctreelist, function(x)
      predict(partykit::as.party(x), newdata = newdata, type = "node"))
  RulesTestMatList <- lapply(1:length(RulesTest), function(x) {
    Test <- RulesTest[[x]]
    rules <- modelout$ruleslist[[x]]
    Test <- factor(as.character(Test), levels = names(rules))
    if (length(unique(names(rules))) > 1) {
      MM <- model.matrix(~ -1 + Test)
      colnames(MM) <-
        paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
      MM
    } else {
      NULL
    }
  })
  TestMat <- Reduce("cbind", RulesTestMatList)
  MM<-model.matrix(~-1+., data=rbind(modelout$datasamp, newdata[,modelout$expvars]))
  MM<-MM[-c(1:nrow(modelout$datasamp)),]
  TestMat<-as.data.frame(cbind(MM, TestMat))
  PredRFOut<-Predict_CRModel_FG(modelout$CRrulefitModel,TestMat)
  return(list(
    TestMat=TestMat,
    PredRFOut=PredRFOut,
    CIFs=PredRFOut$CIFs,
    Times=PredRFOut$Times
  )
  )
}


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
#'
#' @export
CRModel_BART<-function(data,expvars, timevar, eventvar, K=10,ntree=100){
  varprof<-VariableProfile(data, expvars)
  timevarvec<-data[, timevar]
  eventvarvec<-as.integer(data[, eventvar])
  x.train<-as.matrix(model.matrix(~-1+., data=data[, expvars]))
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
#' @param modelout the output from 'Predict_CRModel_BART'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list with the following items.
#' CIFs:predicted CIF matrix (for event coded as 1),
#' Times: The times at which the CIFs are predicted.
#'
Predict_CRModel_BART<-function(modelout, newdata){
  x.test<-as.matrix(model.matrix(~-1+., data=newdata[, modelout$expvars]))
  (M=modelout$post$ndpost)
  (K=modelout$post$K)
  (N=nrow(x.test))
  test = rbind(x.test, x.test)
  pre=BART::crisk.pre.bart(time=modelout$timevarvec, delta=modelout$eventvarvec,
                           x.train= BART::bartModelMatrix(modelout$x.train),
                           x.test=  BART::bartModelMatrix(test),
                           x.train2= BART::bartModelMatrix(modelout$x.train),
                           x.test2=  BART::bartModelMatrix(test),
                           K=K)
  pred = predict(modelout$post, pre$tx.test, pre$tx.test2)
  PredMat<-matrix(pred$cif.test.mean, ncol=K,byrow = TRUE)
  PredMat1<-PredMat[1:N,]
  return(list(CIFs=cbind(0,PredMat1), Times=c(0,modelout$post$times)))
}





#' @title RunCRModels
#' @description  utility function to run several CR models at the same time
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded as 0,1,2),
#' 1 is the event of interest
#' @param models  a vector of models to be fitted to be chosen among
#' "FG", "rulefit", "bart", "cox". Two random forest models are
#' automatically fitted and dont need to be listed here
#' @param varsel  logical indicaiong whether variable selection to be
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
  RF_Model<-CRModel_RF(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)
  ExpVars2<-names(sort(RF_Model$hd.obj$importance[,1], decreasing=TRUE))[1:min(length(ExpVars),20)]
  RF_Model2<-CRModel_RF(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)

 if (!varsel){
  if ("FG" %in% models){
    FG_Model<-tryCatch(CRModel_FineGray(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting FG")
      return()
    })
  }
  if ("bart" %in% models){
    BART_Model<-tryCatch(CRModel_BART(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, ntree=ntreeRF), error=function(e){
      print("Failed fitting BART")
      return()
    })
  }
  if ("cox" %in% models){
    Cox_Model<-tryCatch(CRModel_Cox(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting COX")
      return()
    })
  }

  if ("rulefit" %in% models){
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






#' @title  PredictCRodels
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
    Predict_RF<-Predict_CRModel_RF(models$RF_Model, newdata=newdata)
  }
  if (!is.null(models$RF_Model2)){
    Predict_RF2<-Predict_CRModel_RF(models$RF_Model2, newdata=newdata)
  }
  if (!is.null(models$FG_Model)){
    Predict_FG<-Predict_CRModel_FG(models$FG_Model, newdata=newdata)
  }
  if (!is.null(models$BART_Model)){
    Predict_BART<-Predict_CRModel_BART(models$BART_Model, newdata=newdata)
  }

  if (!is.null(models$Cox_Model)){
    Predict_Cox<-Predict_CRModel_Cox(models$Cox_Model, newdata=newdata)
  }
  if (!is.null(models$rulefit_Model)){
    Predict_rulefit<-Predict_CRModel_rulefit(models$rulefit_Model, newdata=newdata)
  }
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
  NewProbs<-cifMatListAveraging(ModelPredictions)
  NewProbs<-rbind(NewProbs)
  return(list(ModelPredictions=ModelPredictions,NewProbs=NewProbs))
}






#' @title cifInterpolator
#'
#' @noRd
cifInterpolator<-function(x, probs, times){
  f<-approxfun(times, probs)
  sapply(x, function(xi)f(xi))
}

#' @title cifMatInterpolaltor
#'
#' @noRd
cifMatInterpolaltor<-function(probsMat, times,newtimes){
  interpolate1<-function(probs){
    cifInterpolator(newtimes,c(0,probs),c(0,times))
  }
  probsMat1<-apply(probsMat,1,interpolate1)
  probsMat2<-apply(probsMat1,2,function(x){replace(x, ((seq_along(x) <= which(x < cummax(x))[1])|is.na(x)), max(x))})
  probsMat2
}



#' @title cifMatListAveraging
#'
#' @noRd
cifMatListAveraging<-function(listprobsMat, type="CumHaz"){
  if (type=="CumHaz"){
  HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
  for (i in 1:length(listprobsMat)){
    HazzardArray[,,i]<--log(1-listprobsMat[[i]])
  }
  MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(x)))
  NewProbs<-1-exp(-MeanHazzard)
  }
  if (type=="prob"){
    ProbsArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
    for (i in 1:length(listprobsMat)){
      ProbsArray[,,i]<-listprobsMat[[i]]

    }
    NewProbs<-apply(ProbsArray, c(1,2),function(x)(mean(x)))
  }
  NewProbs
}



#' @title timedepConcordanceCR
#'
#' @description Get time dependent concordance for competing risk outcomes (only evemt1)
#' @param predCIF Predicted matrix of probabilities (CIFs)
#' @param predCIFtimes The times for which the probabilities are predicted
#' @param obstimes Observed times
#' @param obsevents Observed event indicator
#' @param TestMat test dataset
#'
#' @return output from the "pec::cindex' function.
#' @export
timedepConcordanceCR<-function(predCIF, predCIFtimes, obstimes, obsevents, TestMat){
  datforpec<-data.frame(time=obstimes,event=obsevents,TestMat)
  cindexTest<-pec::cindex(object=t(predCIF),
                          formula=as.formula("Hist(time, event)~."),
                          data=datforpec,eval.times=predCIFtimes)
  cindexTest
}
