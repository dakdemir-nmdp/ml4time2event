
## Module: Surv Models Functions
## This file contains functions for fitting and predicting survival models.
## Sections include:
##   - Internal utility functions for converting scores and probabilities (score2proba, proba2score)
##   - Integration utility (Integrator)
##   - Survival probability interpolation functions (survivalProbsInterpolator, survprobMatInterpolator, survprobMatListAveraging)
##   - Variable profiling (VariableProfile)
##   - Rule extraction from party objects (listrules)
##   - Survival model functions:
##         * Rulefit model for survival outcomes (SurvModel_rulefit, Predict_SurvModel_rulefit)
##         * GAM model for survival outcomes (SurvModel_GAM, Predict_SurvModel_GAM)
##         * Random Forest survival model (SurvModel_RF, Predict_SurvModel_RF)
##         * Cox model (SurvModel_Cox, Predict_SurvModel_Cox)
##         * Glmnet (SurvModel_glmnet, Predict_SurvModel_glmnet)
##         * GBM model (SurvModel_gbm, Predict_SurvModel_gbm)
##         * XGBoost model (SurvModel_xgboost, Predict_SurvModel_xgboost)
##         * DeepSurv model (SurvModel_deepsurv, Predict_SurvModel_deepsurv)
##         * Parametric survival regression (SurvModel_SurvReg, Predict_SurvModel_SurvReg)
##         * BART model (SurvModel_BART, Predict_SurvModel_BART)
##         - Additional utility functions for time-dependent concordance and Brier scores.
################################################################################
#' @title score2proba
#'
#' @description internal function: from linear score to survival probabilities
#'
#'
#' @noRd
score2proba <-
  function(datasurv, score, conf.int=0.95, which.est=c("point", "lower", "upper")) {
    which.est <- match.arg(which.est)
    pred <- rep(NA, length(score))
    names(pred) <- names(score)
    datacox<-cbind(datasurv, data.frame(score=score))
    colnames(datacox)[1:2]<-c("time","event")
    predm <- survival::coxph(survival::Surv(time, event) ~ score,data=datacox, init=1, control=survival::coxph.control(iter.max = 0))
    sf <- survival::survfit(predm, newdata=data.frame("score"=score), conf.int=conf.int)
    return(list(model=predm,sf=sf))
  }
## End of SurvModelsFunctions.R


################################################################################
#' @title integrate a curve over times
#'
#' @description internal function, function integrator
#' @noRd
Integrator<-function(times, scores, minmax=c(1,35), scale=FALSE){
  mask<-(times>=minmax[1]) & (times<=minmax[2])
  timesn<-times[mask]
  scoressn<-scores[mask]
  AUCsuperlearnMean = pracma::trapz(timesn,scoressn)
  if (scale){AUCsuperlearnMean<-AUCsuperlearnMean/(minmax[2]-minmax[1])}
  AUCsuperlearnMean
}

################################################################################
#' @title survivalProbsInterpolator
#'
#' @noRd
survivalProbsInterpolator<-function(x, probs, times){
  f<-approxfun(times, probs)
  sapply(x, function(xi)f(xi))
}

################################################################################
#' @title survprobMatInterpolator
#'
#' @noRd
survprobMatInterpolator<-function(probsMat, times,newtimes){
  interpolate1<-function(probs){
    y<-survivalProbsInterpolator(newtimes,c(probs),c(times))
    y
  }
  probsMat1<-apply(probsMat,1,interpolate1)
  probsMat2<-apply(probsMat1,2,function(x){replace(x, seq_along(x) <= which(x < cummax(x))[1], max(x))})
  probsMat2
}


################################################################################
#' @title survprobMatListAveraging
#'
#' @noRd
survprobMatListAveraging<-function(listprobsMat){
  HazzardArray<-array(dim=c(dim(listprobsMat[[1]]),length(listprobsMat)))
  for (i in 1:length(listprobsMat)){
    HazzardArray[,,i]<--log(listprobsMat[[i]]+1e-10)
  }
  MeanHazzard<-apply(HazzardArray, c(1,2),function(x)(mean(na.omit(x))))
  NewProbs<-exp(-MeanHazzard)
  NewProbs
}


################################################################################
#' @title proba2score
#'
#' @description internal function: from survival probabilities to  score
#'
#'
#' @noRd
proba2score <-
  function(Probs, Times, ll=0, ul=Inf, scale=FALSE) {
    ul<-min(max(Times),ul)
    ll<-max(min(Times),ll)
    Scores<-apply(Probs, 2,function(scores){Integrator(Times, scores, minmax=c(ll,ul), scale=scale)})
    Scores
  }


################################################################################
#' @title VariableProfile
#'
#' @noRd
VariableProfile<-function(data, expvars){
varprofile<-vector(mode="list", length=length(expvars))
names(varprofile)<-expvars
  for (vari in expvars){
    if (is.factor(data[,vari])){
    varprofile[[vari]]<-table(data[,vari])
    }
    if (is.numeric(data[,vari])){
      varprofile[[vari]]<-c(min(data[,vari]), max(data[,vari]))
    }
    if (is.character(data[,vari])){
      varprofile[[vari]]<-table(data[,vari])
    }
  }
varprofile
}

################################################################################
#' @title listrules
#'
#' @description internal function, extract rules from a party object (a tree)
#'
#'
#' @noRd
listrules<-function (x, i = NULL, ...)
{
  if (is.null(i))
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, listrules, x = x)
    names(ret) <- if (is.character(i))
      i
    else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data
  }
  else {
    fit <- NULL
    dat <- x$data
  }
  rule <- c()
  recFun <- function(node) {
    if (partykit::id_node(node) == i)
      return(NULL)
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index))
        index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) +
          1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels,
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    }
    else {
      if (is.null(index))
        index <- 1:length(kid)
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split),
                                                      Inf))
      sbreak <- breaks[index == whichkid, ]
      right <- partykit::right_split(split)
      srule <- c()
      if (is.finite(sbreak[1]))
        srule <- c(srule, paste(svar, ifelse(right, ">",
                                             ">="), sbreak[1]))
      if (is.finite(sbreak[2]))
        srule <- c(srule, paste(svar, ifelse(right, "<=",
                                             "<"), sbreak[2]))
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}




################################################################################
#' @title SurvModel_rulefit
#'
#' @description Fit a rulefit model for survival outcomes
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
#' TrainMat: fitted rules matrix,
#' ctreelist: list of tree models,
#' yTrain: survival object for outcome variables,
#' cv.fitRules: fitted glmnet survival models for a range of shrinkage parameters,
#' timesTrain: unique times in the outcome variable
#' expvars: explanatory variables used in the model
#'
#'
#'
#' @export
SurvModel_rulefit <- function(data,
                              expvars,
                              timevar,
                              eventvar,
                              ntree = 300,
                              nsample = 300,
                              keepvars = NULL,
                              cuttimes= NULL) {
  varprof<-VariableProfile(data, expvars)

  if (is.null(cuttimes)){
    x<-data[,timevar]
    cuttimes<-quantile(x, c(.1,.25,.50,.70,.90))
  }
  formRF1 <-
    as.formula(paste("survival::Surv(", timevar, ",", eventvar, ") ~.", collapse = ""))
  formRF2 <- as.formula(paste("ClassVar ~.", collapse = ""))

  ctreelist <- lapply(1:ntree, function(repi) {
    sampcols <-
      union(keepvars, sample(expvars, min(length(expvars), sample(c(
        1:10
      ), 1))))
    selmodel <- sample(c(1, 1, 1, 1, 1, 1, 2), 1)
    formRF <- list(formRF1, formRF2)[[selmodel]]
    if (selmodel == 1) {
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
      rpartmodel <- rpart::rpart(formRF,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else {
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
      rpartmodel <- rpart::rpart(formRF,
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

  yTrain = survival::Surv(data[, colnames(data) == timevar], as.numeric(data[, colnames(data) ==
                                                                     eventvar] == 1))
    cv.fitRules = glmnet::cv.glmnet(
      x = TrainMat,
      y = yTrain,
      # create survival object from the data
      alpha = .5,
      # lasso: alpha = 1; ridge: alpha=0
      family = "cox",
      # specify Cox PH model
      maxit = 1000
    )
    est.coefRules = coef(cv.fitRules, s = cv.fitRules$lambda.1se) # returns the p length coefficient vector
    # of the solution corresponding to lambda
    sfitTrainRules <-
      survival::survfit(
        cv.fitRules,
        s = cv.fitRules$lambda.1se,
        x = TrainMat,
        y = yTrain,
        newx = TrainMat
      )
    timesTrain <- sfitTrainRules$time
    nrowforsample<-min(c(nrow(data),500))
    return(
      list(
        datasamp=data[1:nrowforsample,expvars],
        varprof=varprof,
        ruleslist = ruleslist,
        TrainMat = TrainMat,
        ctreelist = ctreelist,
        yTrain = yTrain,
        cv.fitRules = cv.fitRules,
        timesTrain = timesTrain,
        expvars=expvars
        )
    )
}


################################################################################
#' @title Predict_SurvModel_rulefit
#'
#' @description Get predictions from a rulefit model for a test dataset
#'
#'
#' @param modelout the output from 'SurvModel_rulefit'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list object with three items. predSurvsTestRules: survival probability predictions,
#' timesTest: the times of the returned probabilities, sfitTestRules: prediction object
#'
#' @export
Predict_SurvModel_rulefit <- function(modelout, newdata) {

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
  TestMat<-cbind(MM, TestMat)
    sfitTestRules <-
      survival::survfit(
        modelout$cv.fitRules,
        s = modelout$cv.fitRules$lambda.min,
        x = modelout$TrainMat,
        y = modelout$yTrain,
        newx = TestMat
      )
    predSurvsTestRules <- sfitTestRules$surv
    timesTest <- sfitTestRules$time
    if (sum(timesTest==0)==0){
      timesTest<-c(0,timesTest)
      predSurvsTestRules<-rbind(rep(1, ncol(predSurvsTestRules)),predSurvsTestRules)
    }
    return(
      list(
        Probs = predSurvsTestRules,
        Times = timesTest,
        sfitTestRules=sfitTestRules
      )
    )
}



################################################################################
#' @title SurvModel_GAM
#'
#' @description Fit a GAM model for survival outcomes
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param shrinkTreshold integer value, minimum number of
#' factor levels for factor variables to be shrunk
#'
#' @return a list with two elements. b: the fitted GAM model,
#' model: predictions and related output
#'
#' @export
SurvModel_GAM <-
  function(data,
           expvars,
           timevar,
           eventvar,
           shrinkTreshold = 10) {
    varprof<-VariableProfile(data, expvars)
    data[, eventvar] <- as.numeric(data[, eventvar] == 1)
    numvars <-expvars[which(sapply(as.data.frame(data[, expvars]), is.numeric))]
    fctvars <-expvars[which(sapply(as.data.frame(data[, expvars]),function(x) {
        (is.factor(x) | is.character(x))}))]

    catvarstoshrink <-if (length(fctvars)>0){fctvars[which(sapply(as.data.frame(data[, fctvars]), function(x) {
        length(table(x)) > shrinkTreshold
      }))]} else {c()}
    catvarsnottoshrink <-if (length(fctvars)>0){fctvars[which(sapply(as.data.frame(data[, fctvars]), function(x) {
      length(table(x)) <= shrinkTreshold
    }))]} else {c()}
    numvarstosmooth <-if (length(numvars)>0){numvars[which(sapply(as.data.frame(data[, numvars]), function(x) {
      length(table(x)) > shrinkTreshold
    }))]} else {c()}
    numvarsnottosmooth <-if (length(numvars)>0){numvars[which(sapply(as.data.frame(data[, numvars]), function(x) {
      length(table(x)) <= shrinkTreshold
    }))]} else {c()}
    for (i in fctvars) {
      data[, i] <- as.factor(data[, i])
    }
    formGAM <- "~"
    for (vari in catvarstoshrink) {
      formGAM <-
        paste(formGAM, paste("s(", vari, ", bs='re')", sep = ""), sep = "+")
    }
    for (vari in numvarstosmooth) {
      formGAM <- paste(formGAM, paste("s(", vari, ")", sep = ""), sep = "+")
    }
    for (vari in c(numvarsnottosmooth, catvarsnottoshrink)) {
      formGAM <- paste(formGAM, vari, sep = "+")
    }
    formGAM <- paste(timevar, gsub("\\~\\+", "\\~", formGAM), sep = "")
    print(formGAM)

    b <- mgcv::gam(
      as.formula(formGAM),
      family = mgcv::cox.ph(),
      data = data,
      weights = data[, eventvar],
      select = TRUE
    )
    fvTrain <- predict(b,
                       newdata = data[, expvars],
                       type = "link",
                       se = TRUE)
    datasurv <- data[, c(timevar, eventvar)]
    times <- sort(unique(data[data[, eventvar] == 1, timevar]))
    preds <-
      score2proba(
        datasurv = datasurv,
        score = fvTrain$fit,
        conf.int = 0.95,
        which.est = "point"
      )
    return(list(
      b = b,
      model = preds,
      varprof=varprof
    ))
  }

################################################################################
#' @title Predict_SurvModel_GAM
#'
#' @description Get predictions from a GAM model for a test dataset
#'
#' @param modelout the output from 'SurvModel_GAM'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list object with three items. estSURVTest: survival probability predictions,
#' time.interest: the times of the returned probabilities
#' sf: prediction object
#'
#' @export
Predict_SurvModel_GAM <- function(modelout, newdata) {
  fvTest <- predict(modelout$b,
                    newdata = newdata,
                    type = "link",
                    se = TRUE)
  sf <-
    survival::survfit(
      modelout$model$model,
      newdata = data.frame("score" = fvTest$fit),
      conf.int = .95
    )
  estSURVTest = sf$surv
  time.interest = sf$time
  if (sum(time.interest==0)==0){
    time.interest<-c(0,time.interest)
    estSURVTest<-rbind(rep(1, ncol(estSURVTest)),estSURVTest)
  }
  return(list(
    Probs = estSURVTest,
    Times = time.interest,
    sf=sf
  ))
}








#' @title SurvModel_RF
#'
#' @description Fit a RF model for survival outcomes
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param ntree integer value,  number of trees to grow
#' @param samplesize integer value,  sample size for each grown tree
#' @param nsplit integer value, maximum number of splits for each tree
#' @param trace trace or not, logical
#'
#' @return the output of 'randomForestSRC::rfsrc'
#'
#' @export
SurvModel_RF<-function(data,expvars, timevar, eventvar, ntree=300, samplesize=500, nsplit=5, trace=TRUE){
  formRF<-as.formula(paste("Surv(",timevar, ",", eventvar,") ~.", collapse = ""))
  data[,eventvar]<-as.numeric(data[,eventvar]==1)
  varprof<-VariableProfile(data, expvars)
  for (vari in expvars){
    if (is.character(data[, vari])){
      data[, vari]<-as.factor(data[, vari])
    }
  }
  o <- randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars)],
                             splitrule="bs.gradient",samptype = "swor",sampsize = samplesize, trace = trace,nsplit=nsplit,stepFactor = 1.5,
                             mtryStart = 2,
                             nodesizeTry = c(seq(1, 101, by = 10)), ntreeTry = ntree)
  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars)],
                                   nodesize =o$optimal[[1]], ntree=ntree,mtry= o$optimal[[2]],
                                   tree.err = FALSE, importance = TRUE,statistics=TRUE,
                                   do.trace = trace, splitrule="bs.gradient",samptype = "swor",sampsize = samplesize, nsplit = nsplit)#, case.wt=pvec)
  return(list(hd.obj=hd.obj, varprof=varprof))
}

#' @title Predict_SurvModel_RF
#'
#' @description Get predictions from a RF model for a test dataset
#'
#' @param modelout the output from 'SurvModel_RF'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#'  Probs=predicted Survival probability matrix,
#'  Times: The times at which the CIFs are predicted,
#'  predSurvsTestRF:the output of 'randomForestSRC::predict.rfsrc'.
#'
#' @export
Predict_SurvModel_RF<-function(modelout, newdata){
  predSurvsTestRF<-randomForestSRC::predict.rfsrc(modelout$hd.obj, newdata = newdata)
  Probs<-cbind(1,predSurvsTestRF$survival)
  Times<-c(0,predSurvsTestRF$time.interest)
  return(list(
    Probs = t(Probs),
    Times = Times,
    predSurvsTestRF=predSurvsTestRF)
    )
}



#' @title SurvModel_Cox
#'
#' @description Fit a Cox model for survival outcomes using forward selection
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#'
#' @return a list of two items. cph_model: coxPH model object,
#'  times: unique times in the data
#' @export
SurvModel_Cox<-function(data,expvars, timevar, eventvar){
  form<-as.formula(paste(paste("survival::Surv(",timevar, ",", eventvar,") ~1+", collapse = ""),paste(expvars,collapse="+"),sep=""))
  print(form)
  varprof<-VariableProfile(data, expvars)

  data[,eventvar]<-as.numeric(data[,eventvar]==1)
  XYTrain<-data[,c(timevar,eventvar,expvars)]
  cph_model <- survival::coxph(form, data = XYTrain,x=TRUE,y=TRUE)
  stepMod<-stats::step(cph_model)
  times<-as.numeric(c(unlist(XYTrain[, timevar])))
  cph_modelTrainPredict<- survival::survfit(stepMod, newdata=XYTrain,times=sort(unique(times)))
  return(list(cph_model=stepMod,times=sort(unique(times)), varprof=varprof))
}

#' @title Predict_SurvModel_Cox
#'
#' @description Get predictions from a Cox model for a test dataset
#'
#' @param modelout the output from 'SurvModel_Cox'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @returns 'cph_modelTestPredict': predictions as the output of 'survival::survfit'
#'
#' @export
Predict_SurvModel_Cox<-function(modelout, newdata){
  cph_modelTestPredict<- survival::survfit(modelout$cph_model, newdata=newdata,times=modelout$times)
  Probs<-cph_modelTestPredict$surv
  Times<-cph_modelTestPredict$time
  return(list(Probs=rbind(1,Probs),
              Times=c(0,Times),
              cph_modelTestPredict=cph_modelTestPredict))
}





#' @title SurvModel_glmnet
#' @description Fit a Cox model for survival outcomes using glmnet
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param alpha   lasso: alpha = 1; ridge: alpha=0
#' @param maxit maximum number of iterations for glmnet
#'
#' @return a list containing the following objects. traindata: a sample of training data,
#' TrainMat: preprocessed version of training data,
#' yTrain: a sample for response variable,
#' expvars: a character vector of explanatory variables,
#' cv.fit: glmnet cv fit object,
#' est.coef: estimated model coefficients
#' @export
SurvModel_glmnet<-function(data,expvars, timevar, eventvar, alpha=.5, maxit=5000, nfolds=30){
  formRF<-as.formula(paste("survival::Surv(",timevar, ",", eventvar,") ~.", collapse = ""))
  varprof<-VariableProfile(data, expvars)

   data[,eventvar]<-as.numeric(data[,eventvar]==1)
  TrainMat<-model.matrix(~., data=data[,expvars])
  yTrain = survival::Surv(data[, colnames(data)==timevar], data[, colnames(data)==eventvar])
  cv.fit = glmnet::cv.glmnet(x = TrainMat,
                     y = yTrain, # create survival object from the data
                     alpha = alpha, # lasso: alpha = 1; ridge: alpha=0
                     family = "cox", # specify Cox PH model
                     maxit = maxit,
                     nfolds=nfolds
                     )
  est.coef = stats::coef(cv.fit, s = cv.fit$lambda.min) # returns the p length coefficient vector
  # of the solution corresponding to lambda
  sfitTrain<-survival::survfit(cv.fit, s =cv.fit$lambda.min, x = TrainMat, y=yTrain, newx=TrainMat)
  predSurvsTrain<-sfitTrain$surv
  lpTrain<-TrainMat%*%est.coef
  timesTrain<-sfitTrain$time
  sampledata<-sample(1:nrow(data),min(c(500,nrow(data))))
  out<-list(traindata=data[sampledata,expvars],TrainMat=TrainMat[sampledata,],yTrain=yTrain[sampledata,], expvars=expvars,cv.fit=cv.fit, est.coef=est.coef, varprof=varprof)
  return(out)
}



#' @title Predict_SurvModel_glmnet
#'
#' @description Get predictions from a glmnet model for a test dataset
#'
#' @param modelout the output from 'SurvModel_glmnet'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items: predSurvsTest: survival probability matrix,
#' timesTest: the unique times for which the probabilities are calculated
#'
#' @export
Predict_SurvModel_glmnet<-function(modelout, newdata){
  mmdata<-data.frame(rbind(modelout$traindata,newdata[,modelout$expvars]))
  TestMat<-model.matrix(~., data=mmdata)[-c(1:nrow(modelout$traindata)),]
  sfitTest<-survival::survfit(modelout$cv.fit, s=modelout$cv.fit$lambda.min, x = modelout$TrainMat, y=modelout$yTrain, newx=TestMat)
  predSurvsTest<-sfitTest$surv
  timesTest<-sfitTest$time
if (sum(timesTest==0)==0){
  timesTest<-c(0,timesTest)
  predSurvsTest<-rbind(1,predSurvsTest)
}
    return(list(Probs=predSurvsTest,Times=timesTest))
}



#' @title SurvModel_gbm
#'
#' @description Fit a gbm model for survival outcomes using glmnet
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param ntree number of trees to grow
#' @param max.depth maximum depth for the trees
#' @param bag.fraction fraction of data to sample in each bagging iterations
#' @param train.fraction fraction for training data
#' @param learninrate learning rate for the boosting algorithm
#'
#' @return a list which contains the following objects.
#' gbmmodel: model object,
#' best.iter: best iteration number
#' expvars: vector of the names of explanatory variables,
#' survMat: survival probabilities,
#' basehaz.cum: baseline cumulative hazard
#' time.interest: unique observed event times
#'
#'
#' @export
SurvModel_gbm<-function(data,expvars, timevar, eventvar, ntree=200, max.depth=3, bag.fraction=.4, train.fraction=.4, learninrate=.01){
  data[,eventvar]<-as.numeric(data[,eventvar]==1)
  numvars<-names(which(sapply(as.data.frame(data[,expvars]), is.numeric)))
  fctvars<-names(which(sapply(as.data.frame(data[,expvars]), function(x){(is.factor(x) | is.character(x))})))
  varprof<-VariableProfile(data, expvars)

   for (i in fctvars){
    data[,i]<-as.factor(data[,i])
  }
  gbmmodel <- gbm::gbm(survival::Surv(data[, timevar], data[, eventvar]) ~ .,
                   data=data[, colnames(data)%in%c(expvars)],
                   distribution="coxph",
                   n.trees=5,
                   shrinkage=learninrate,
                   interaction.depth=max.depth,
                   bag.fraction = bag.fraction,
                   train.fraction = train.fraction,
                   cv.folds = 2,
                   n.minobsinnode = 30,
                   keep.data = TRUE)
  best.iter <- gbm::gbm.perf(gbmmodel,method="cv")
  time.interest <- sort(unique(data[, timevar][data[, eventvar]==1]))
  pred.train <- predict(gbmmodel, data, n.trees = best.iter)
  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum <- gbm::basehaz.gbm(data[, timevar], data[, eventvar], pred.train, t.eval = time.interest, cumulative = TRUE)
  survMat<-NULL
  for (i in 1:length(pred.train)){
    surf.i <- exp(-exp(pred.train[i])*basehaz.cum)
    survMat<-rbind(survMat,surf.i)
  }
  return(list(gbmmodel=gbmmodel, best.iter=best.iter,expvars=expvars, survMat=survMat,basehaz.cum=basehaz.cum,time.interest=time.interest, varprof=varprof))

}




#' @title Predict_SurvModel_gbm
#'
#' @description Get predictions from a gbm model for a test dataset
#'
#' @param modelout the output from 'SurvModel_glmnet'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items: survMat: survival probability matrix,
#' time.interest: the unique times for which the probabilities are calculated
#'
#' @export
Predict_SurvModel_gbm<-function(modelout, newdata){
  data<-newdata[, colnames(newdata)%in%c(modelout$expvars)]
  for (i in 1:ncol(data)){
    if (!is.null(levels(data[,i]))){
      data[,i]<-as.character(data[,i])
    }
  }
  # # Prediction and evaluation metrics
  event_prediction <- suppressWarnings(predict(modelout$gbmmodel,data,n.trees=modelout$best.iter))
  #
  # surv <- Surv(data.test$time, data.test$event)
  survMat<-NULL
  for (i in 1:length(event_prediction)){
    surf.i <- exp(-exp(event_prediction[i])*modelout$basehaz.cum)
    survMat<-rbind(survMat,surf.i)
  }
  return(list(Probs=t(cbind(1,survMat)),Times=c(0,modelout$time.interest)))
}



#' @title xgb.train.surv
#' @description internal function
#' @noRd
xgb.train.surv <- function(params = list(), data, label, weight = NULL, nrounds,
                           watchlist = list(), verbose = 1, print_every_n = 1L,
                           early_stopping_rounds = NULL, save_period = NULL,
                           save_name = "xgboost_surv.model", xgb_model = NULL, callbacks = list(), ...) {

  if (length(params) > 0) {
    if (params$objective != "survival:cox") stop("params objective must be set to survival:cox")
    if (params$eval_metric != "cox-nloglik") stop("params eval_metric must be set to cox-nloglik")
  } else {
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik"
    )
  }

  if(is.null(weight)) weight <- rep(1, nrow(data))

  data_DMatrix <- xgboost::xgb.DMatrix(data = data, label = label, weight = weight)

  xgboost_model <- xgboost::xgb.train(
    params = params, data = data_DMatrix, nrounds = nrounds, watchlist = watchlist, verbose = verbose,
    print_every_n = print_every_n, early_stopping_rounds = early_stopping_rounds, save_period = save_period,
    save_name = "surv_xgboost.model", xgb_model = xgb_model, callbacks = callbacks, ...
  )


  # generate baseline hazard
  data_data.frame <- data.frame(data, time = abs(label), status = ifelse(sign(label) == 1, 1, 0))

  cox_model <- survival::coxph(formula = survival::Surv(time, status) ~ ., data = data_data.frame)
  baseline_hazard <- survival::basehaz(cox_model)

  if (baseline_hazard[1, 2] != 0) {
    baseline_hazard <- rbind(c(0, 0), baseline_hazard) # pec always requests time = 0 survival as well
  }

  HR <- xgboost:::predict.xgb.Booster(object = xgboost_model, newdata = data_DMatrix)
  baseline_pred <- function(const) {
    risk <- HR * const
    surv <- exp(risk %*% -matrix(baseline_hazard[, 1], nrow = 1))

    Models <- list(
      "xgboost" = surv
    )

    PredError<-pec::pec(object=Models,formula=Surv(time,status)~1,
                        data=data_data.frame,
                        cens.model="marginal",
                        times = baseline_hazard[, 2],
                        exact = F,
                        verbose = F,
                        reference = F,
                        splitMethod = "none"
    )

    return(pec::crps(PredError))
  }

  optimal_const <- optim(par = 1, fn = baseline_pred, method = "Brent", lower = 0, upper = 10)
  baseline_hazard[, 1] <- baseline_hazard[, 1] * optimal_const$par
  xgboost_model$baseline_hazard <- baseline_hazard
  class(xgboost_model) <- "xgb.Booster"
  return(xgboost_model)
}

#' @title xgb.train.surv
#' @description internal function
#' @noRd
predict.xgb.Booster.surv <- function(object, newdata, type = "risk", times = NULL) {
  risk <- xgboost:::predict.xgb.Booster(object, newdata)
  if (type == "risk") {
    return(risk)
  } else if (type == "surv") {
    if (!is.null(times)) {
      if (max(times) > max(object$baseline_hazard[, 2])) {
        object$baseline_hazard <- rbind(object$baseline_hazard, c(max(object$baseline_hazard[, 1]), max(times)))
      }
    } else {
      times <- object$baseline_hazard[, 2]
    }
    surv <- exp(risk %*% -matrix(object$baseline_hazard[, 1], nrow = 1))
    surv <- surv[, findInterval(times, object$baseline_hazard[, 2])]
    colnames(surv) <- times
    return(surv)
  } else {
    stop('type must be one of "risk", "surv"')
  }
}





#' @title SurvModel_xgboost
#'
#' @description Fit a xgboost model for survival outcomes using xgboost
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#'
#'
#' @return a list containing the following objects.
#' expvars: vector of explanatory variables,
#' learner: fitted model object,
#' estSURVTrain: predicted probabilities,
#' datatrainProf: sample data for internal use
#'
#' @export
SurvModel_xgboost<-function(data,expvars, timevar, eventvar){
  varprof<-VariableProfile(data, expvars)
  X<-model.matrix(~-1+., data[,c(expvars)])
  colnames(X)<-paste("v",1:ncol(X), sep="")
  ytimes = data[, c(timevar)]
  yevents<-as.numeric(data[, c(eventvar)])
  yTrain<-ytimes
  yTrain[yevents==0]<-(-yTrain[yevents==0])
  times<-sort(unique(data[data[,eventvar]==1, timevar]))

  params <- list(objective='survival:cox',
                 tree_method='hist',
                 eval_metric = "cox-nloglik",
                 learning_rate=0.01,
                 max_depth=5,
                 early_stopping_rounds=10)

  dtrain <- xgboost::xgb.DMatrix(data=X)
  bst=xgb.train.surv(data=X,label=yTrain,params = params, nrounds=100)
  preds<-predict.xgb.Booster.surv(bst, dtrain,type = "surv", times=sort(unique(c(0,times))))
  return(list(expvars=expvars, learner=bst,estSURVTrain=preds, datatrainProf=data[1:10,c(expvars)], varprof=varprof))
}







#' @title Predict_SurvModel_xgboost
#'
#' @description Get predictions from a xgboost model for a test dataset
#'
#' @param modelout the output from 'SurvModel_xgboost'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items: estSURVTest: survival probability matrix,
#' time.interest: the unique times for which the probabilities are calculated
#'
#' @export
Predict_SurvModel_xgboost<-function(modelout, newdata){
 for (vari in modelout$expvars){
   if (is.factor(modelout$datatrainProf[, vari])){
   newdata[,vari]<-factor(as.character(newdata[,vari]), levels=levels(modelout$datatrainProf[, vari]))
   }
  }
  X<-as.matrix(model.matrix(~-1+., newdata[,c(modelout$expvars)]))
  colnames(X)<-paste("v",1:ncol(X), sep="")
  dtest <- xgboost::xgb.DMatrix(data=X)
  times=as.numeric(colnames(modelout$estSURVTrain))
  preds<-predict.xgb.Booster.surv(modelout$learner, dtest,type = "surv", times=sort(unique(c(0,times))))
  return(list(Probs=t(preds),Times=sort(unique(c(0,times)))))
}





#' @title get_avesurv.deepsurv
#' internal function
#' @noRd
get_avesurv.deepsurv <- function(object, ...) {
  object <- get_indivsurv.deepsurv(object)
  surv <- as.vector(colMeans(object$surv))
  chaz <- -log(surv)
  time <- object$time
  out <- list(time = time, surv = surv, chaz=chaz)
  out$call <- match.call()
  class(out) <- "satsurv"
  out
}

#' @title get_indivsurv.deepsurv
#' @description internal function
#' @noRd
get_indivsurv.deepsurv <- function(object, newdata) {
  if (missing(newdata)) {
    newdata <- object$modelData
  }
  surv <- predict(object, newdata = newdata, type = "survival")
  times <- as.numeric(colnames(surv))
  #times <- times[1:(length(times)-1)]
  ss <- surv[, 1:length(times),drop=FALSE]
  out <- list(time = times, surv = ss, chaz = -log(ss))
  out$call <- match.call()
  class(out) <- "satsurv"
  out
}


#' @title predictSurvProb.deepsurv
#' @description internal function
#' @noRd
predictSurvProb.deepsurv <- function(object, newdata, times, ...){
  N <- nrow(newdata)
  sfit <- get_indivsurv.deepsurv(object, newdata = newdata)
  S <- sfit$surv
  Time <- sfit$time
  S <- matrix(S, ncol=length(Time))
  if(N == 1) S <- matrix(S, nrow = 1)
  p <-  cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
  if (nrow(p) != nrow(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


#' @title SurvModel_deepsurv
#'
#' @description Fit a deepsurv model for survival outcomes
#'  To use first do:
#'
#'  library(survivalmodels)
#'
#'  install_pycox(pip = TRUE, install_torch = TRUE)
#'
#'  install_keras(pip = TRUE, install_tensorflow = TRUE)
#'
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param num_nodes number of nodes in each of the hidden layers
#' @param dropout  dropout percentage for each hiden layer
#'
#' @return a list contating the following objects:
#' expvars: character vector of explanatory variable names
#' learner: fitted learner
#' times: observed event times
#'
#' @seealso [survivalmodels::deepsurv()] which this function wraps.
#' @export
SurvModel_deepsurv<-function(data,expvars, timevar, eventvar, num_nodes = c(100L,100L),dropout = 0.2, batch_size=128L,  epochs = 100L){
  #########################
  varprof<-VariableProfile(data, expvars)
  DatSurvXYTrain<-as.data.frame(model.matrix(~-1+.,model.frame(~-1+.,data[,c(timevar, eventvar, expvars)])))
  colnames(DatSurvXYTrain)[1:2]<-c("time","status")

  tmp<-capture.output(moddeepSurv<-survivalmodels::deepsurv(data = DatSurvXYTrain, frac = 0, activation = "relu",
                                                            num_nodes = num_nodes, dropout = dropout, early_stopping = FALSE, epochs = 100L,
                                                            batch_size = 256L, verbose = 0, batch_norm = TRUE,weight_decay=1e-4, lr=1e-5))
  return(list(expvars=expvars, learner=moddeepSurv, times=sort(unique(c(DatSurvXYTrain$time))), varprof=varprof))
}



#' @title Predict_SurvModel_deepsurv
#'
#' @description Get predictions from a deepsurv model for a test dataset
#'
#' @param modelout the output from 'SurvModel_glmnet'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following objects. estSURVTest: predicted survival probability matrix
#' time.interest: times at which the survival probabilities are provided.
#'
#' @export
Predict_SurvModel_deepsurv<-function(modelout, newdata){
  DatSurvXTest<-as.data.frame(model.matrix(~-1+.,model.frame(~-1+.,data=newdata[,c(modelout$expvars)])))
  preddeepSurvTest<-predictSurvProb.deepsurv(modelout$learner,newdata = DatSurvXTest,times=modelout$times)
  if (sum(modelout$times==0)==0){
    modelout$times<-c(0,modelout$times[1:(length(modelout$times)-1)])
    preddeepSurvTest<-preddeepSurvTest
  }
  return(list(Probs=t(preddeepSurvTest),Times=modelout$times))
}




#' @title savedeepsurv
#' @description Save a SurvModel_deepsurv model using ave_model_tf
#' @param modelout the output from 'SurvModel_deepsurv'
#' @param filename the name of the file to save the model to
#'
#' @export
savedeepsurv<-function(modelout, filename){
  deepsurvM$learner$model$save_net(filename)
  # save the net as pickle using reticulate
  reticulate::py_save_object(modelout$learner$model, filename)
  filenameRData<-paste(filename, ".RData", sep = "")
  save(modelout, file=filenameRData)
}

#' @title loaddeepsurv
#' @description Load a SurvModel_deepsurv model using load_model_tf
#' @param filename the name of the file to load the model from
#'
#' @return a SurvModel_deepsurv model
#' @export
loaddeepsurv<-function(filename){
  filenamerdata=paste(filename, ".RData", sep = "")
  load(filenamerdata)
  model <- reticulate::py_load_object(filename)
  modelout$learner$model<-model
  modelout$learner$model$load_net(filename)
  return(modelout)
}

#' @title SurvModel_SurvReg
#'
#' @description Fit parametric survival models
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param dist one of "weibull", "exponential", "gaussian", "logistic","lognormal" or "loglogistic"
#'
#' @return a list containing the following objects. survregOut: survreg model out,
#' predicttrain: predictions in the training data,
#' times.interest: unique times for which the probabilities are calculated,
#' expvars: explanatory variables used to build the model
#'
#' @export

SurvModel_SurvReg <- function(data, expvars, timevar, eventvar, dist = "exponential") {
  varprof <- VariableProfile(data, expvars)

  form <- as.formula(paste(paste("survival::Surv(", timevar, ",", eventvar, ") ~", collapse = ""), paste(expvars, collapse = "+"), sep = ""))
  data[, eventvar] <- as.numeric(data[, eventvar] == 1)
  dataYX <- data[, c(timevar, eventvar, expvars)]

  # Check if at least one variable is present
  if (length(expvars) == 0) {
    stop("No variables provided. At least one variable is required.")
  }

  # Forward selection with AIC
  selected_vars <- c()
  best_model <- NULL
  best_aic <- Inf

  while (length(selected_vars) < length(expvars)) {
    aic_values <- numeric()

    for (var in setdiff(expvars, selected_vars)) {
      formula <- as.formula(paste(paste("survival::Surv(", timevar, ",", eventvar, ") ~", collapse = ""), paste(c(selected_vars, var), collapse = "+"), sep = ""))
      model <- survreg(formula, data = dataYX, dist = dist, x = TRUE, y = TRUE)
      aic_values <- c(aic_values, AIC(model))
    }

    best_var <- setdiff(expvars, selected_vars)[which.min(aic_values)]
    best_aic_new <- aic_values[which.min(aic_values)]

    if (length(selected_vars) == 0 && best_aic_new >= best_aic) {
      best_var <- expvars[which.min(aic_values)]
    }

    if (best_aic_new < best_aic) {
      selected_vars <- c(selected_vars, best_var)
      best_model <- survreg(as.formula(paste(paste("survival::Surv(", timevar, ",", eventvar, ") ~", collapse = ""), paste(selected_vars, collapse = "+"), sep = "")), data = dataYX, dist = dist, x = TRUE, y = TRUE)
      best_aic <- best_aic_new
    } else {
      break
    }
  }

  if (is.null(best_model)) {
    stop("No valid model could be fit. Please check the data and variables.")
  }

  times.interest <- sort(unique(dataYX[, timevar]))

  predicttrain <- predict(best_model, newdata = dataYX, type = "quantile", p = seq(0, 1, length = 100))
  predicttrain <- apply(predicttrain, 1, function(x) { survivalProbsInterpolator(c(0, times.interest), 1 - seq(0, 1, length = 100), x) })

  return(list(survregOut = best_model, predicttrain = predicttrain, times.interest = c(0, times.interest), expvars = selected_vars, varprof = varprof))
}


#' @title Predict_SurvModel_SurvReg
#'
#' @description Make predictions using BART model, trained by the
#' SurvModel_SurvReg' function
#' @param modelout the output from 'SurvModel_SurvReg'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following objects. estSURVTest: predicted survival probability matrix,
#' time.interest: the times at which the probabilities are calculated
#'
#' @export
Predict_SurvModel_SurvReg <- function(modelout, newdata, times = NULL) {
  X <- newdata[, colnames(newdata) %in% modelout$expvars]

  predicttest <- predict(modelout$survregOut, newdata = X, type = "quantile", p = seq(0, 1, length = 100))
  if (is.null(times)) {
    times.interest <- modelout$times.interest
  } else {
    times.interest <- times
  }
  predicttest <- apply(predicttest, 1, function(x) {survivalProbsInterpolator(c(times.interest), 1 - seq(0, 1, length = 100), x) })
  return(list(Probs = predicttest, Times = c(times.interest)))
}



#' @title SurvModel_BART
#'
#' @description Fit a BART model
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param K  parameter 'K' for BART
#' @export ntree number of trees in BART model

SurvModel_BART<-function(data,expvars, timevar, eventvar, K=8, ntree=50){
  timevarvec<-data[, timevar]
  eventvarvec<-as.integer(data[, eventvar])
  varprof<-VariableProfile(data, expvars)
  x.train<-as.matrix(model.matrix(~-1+., data=data[, expvars]))
  failcount<-0
  post<-NULL
  while ((is.null(post) & (failcount<4))){
    failcount<-failcount+1
    post <- tryCatch(BART::surv.bart(x.train=x.train, times=timevarvec, delta=eventvarvec, x.test=x.train,
                                     K=K,ntree=ntree,ndpost=2000, nskip=500,
                                     keepevery = 2L),
                     error=function(e){NULL}
    )
  }
  return(list(post=post, expvars=expvars,eventvarvec=eventvarvec,timevarvec=timevarvec,x.train=x.train, varprof=varprof ))
}


#' @title Predict_SurvModel_BART
#'
#' @description Make predictions using BART model, trained by the 'SurvModel_BART' function
#' @param modelout the output from 'SurvModel_BART'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following objects. PredMat: predicted survival probability matrix,
#' times: the times at which the probabilities are calculated.
#'
#' @export
Predict_SurvModel_BART<-function(modelout, newdata, times=NULL){
  x.test<-as.matrix(model.matrix(~-1+., data=newdata[, modelout$expvars]))
  (M=modelout$post$ndpost)
  (K=modelout$post$K)
  (N=nrow(x.test))
  Xall<-dplyr::bind_rows(as.data.frame(x.test),as.data.frame(modelout$x.train))
  x.test<-Xall[1:N,]
  x.test<-x.test[,colnames(modelout$x.train)]
  x.test[is.na(x.test)]<-0
  pre=BART::surv.pre.bart(times=modelout$timevarvec, delta=modelout$eventvarvec,
                          x.train= BART::bartModelMatrix(modelout$x.train),
                          x.test=  BART::bartModelMatrix(x.test),
                          K=K)
  pred = predict(modelout$post, pre$tx.test)
  PredMat<-matrix(pred$surv.test.mean, ncol=K,byrow = TRUE)
  return(list(Probs=t(cbind(1,PredMat)), Times=c(0,modelout$post$times)))
}



####


#' @title timedepConcordance
#'
#' @description Get time dependent concordance
#' @param predsurv Predicted matrix of probabilities
#' @param predsurvtimes The times for which the probabilities are predicted
#' @param obstimes Observed times
#' @param obsevents Observed event indicator
#' @param ctimes The times at which the concordance index is calculated
#'
#' @return output from the "pec::cindex' function.
#' @export
timedepConcordance<-function(predsurv, predsurvtimes, obstimes, obsevents, ctimes=NULL){
  datforpec<-data.frame(time=obstimes,event=obsevents)
  if (is.null(ctimes)){
    ctimes<-predsurvtimes
  }
  cindexTest<-pec::cindex(object=t(predsurv),
                          formula=as.formula("Surv(time, event)~."),
                          data=datforpec,eval.times=ctimes)
  cindexTest
}


#' @title integratedC
#'
#' @description Integrated c-score
#'
#' @param times The times at which the concordance scores were calculated
#' @param scores Time dependent concordance scores
#' @param minmax The min max limits of the integral
#'
#' @return a scalar measure of accuracy
#' @export
integratedC<-function(times, scores, minmax=c(1,35)){
  mask<-(times>=minmax[1]) & (times<=minmax[2])
  timesn<-times[mask]
  scoressn<-scores[mask]
  AUCMean = pracma::trapz(timesn,scoressn)/(minmax[2]-minmax[1])
  AUCMean
}




#' @title BrierScore
#'
#' @description get time dependent brier score (PEC)
#'
#' @param predsurv Predicted matrix of probabilities
#' @param predsurvtimes The times for which the probabilities are predicted
#' @param obstimes Observed times
#' @param obsevents Observed event indicator
#'
#' @return the output from the 'pec::pec' function
#' @export
BrierScore<-function(predsurv, predsurvtimes, obstimes, obsevents){
  datforpec<-data.frame(time=obstimes,event=obsevents)
  BrierTest<-pec::pec(object=t(predsurv),
                      formula=Surv(time, event)~1,
                      data=datforpec, reference=FALSE)
  BrierTest
}



#' @title integratedBrier
#'
#' @description Integrated Brier score
#'
#' @param times The times at which the concordance scores were calculated
#' @param scores Time dependent Brier scores
#' @param minmax The min max limits of the integral
#'
#' @return a scalar measure of accuracy
#'
#' @export
integratedBrier<-function(times, scores, minmax=c(1,35)){
  mask<-(times>=minmax[1]) & (times<=minmax[2])
  timesn<-times[mask]
  scoressn<-scores[mask]
  AUCMean = pracma::trapz(timesn,scoressn)/(minmax[2]-minmax[1])
  AUCMean
}



##############


#' @title RunSurvModels
#'
#' @description utility function to run several survival models at the same time
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data
#' @param models  a vector of models to be fitted to be chosen among
#' "glmnet","coxph","rulefit","xgboost","deepsurv","gam","gbm",
#' "ExpSurvReg","WeibSurvReg","bart". Two random forest models are
#' automatically fitted and dont need to be listed here
#' @return a list of two items. First is a list containing ExpVars: all of the explanatory vars,
#'  ExpVars2: a subset of the explanatory variables selected by random forest,
#'  timevar: time variable,
#'  eventvar: event variable. The second is also a list that contain the individual model outputs.
#' @export
RunSurvModels<-function(datatrain, ExpVars, timevar, eventvar, models=c("glmnet","coxph","rulefit","xgboost","deepsurv","gam","gbm","ExpSurvReg","WeibSurvReg","bart"), ntreeRF=300,nvars=20, ...){
  datatrainFact<-datatrain
  for (i in 1:ncol(datatrainFact)){
    if (is.character(datatrainFact[,i])){
      datatrainFact[,i]<-as.factor(datatrainFact[,i])
    }
  }
  RF_Model<-SurvModel_RF(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)
  ExpVars2<-names(sort(RF_Model$hd.obj$importance, decreasing=TRUE))[1:min(length(ExpVars), nvars)]
  RF_Model2<-SurvModel_RF(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar, samplesize=min(c(500,ceiling(.3*nrow(datatrainFact)))), ntree=ntreeRF)

  if ("bart"%in% models){
    bartout<-tryCatch(SurvModel_BART(datatrain,ExpVars, timevar=timevar, eventvar=eventvar, ntree=50), error=function(e){
      print("Failed fitting bart")
      return()
    })
  }

  if ("ExpSurvReg" %in% models) {
    ExpVarsW <- ExpVars2
    while (length(ExpVarsW) > 0) {
      survregexp_Model <- tryCatch(
        SurvModel_SurvReg(datatrain, ExpVarsW, timevar = timevar, eventvar = eventvar, dist = "exponential"),
        error = function(e) {
          print("Failed fitting ExpSurv")
          return()
        }
      )
      if (!is.null(survregexp_Model)) {
        break
      }
      ExpVarsW <- ExpVarsW[-length(ExpVarsW)]
    }
  }

  if ("WeibSurvReg" %in% models) {
    ExpVarsW <- ExpVars2
    while (length(ExpVarsW) > 0) {
      survregweib_Model <- tryCatch(
        SurvModel_SurvReg(datatrain, ExpVarsW, timevar = timevar, eventvar = eventvar, dist = "weibull"),
        error = function(e) {
          print("Failed fitting WeibSurv")
          return()
        }
      )
      if (!is.null(survregweib_Model)) {
        break
      }
      ExpVarsW <- ExpVarsW[-length(ExpVarsW)]
    }
  }


  if ("glmnet" %in% models){
    glmnet_Model<-tryCatch(SurvModel_glmnet(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting glmnet")
      return()

    })
  }
  if ("coxph" %in% models){
    CPH_Model<-tryCatch(SurvModel_Cox(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Error in CPH")
      return()
      })
  }
  if ("rulefit"%in% models){

    RuleFit_Model<-tryCatch(SurvModel_rulefit(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar, ntree=ntreeRF,nsample=min(c(500,ceiling(.3*nrow(datatrainFact))))), error=function(e){
      print("Failed fitting Rulefit")
      return()
    })
  }
  if ("xgboost"%in% models){
    xgboost_Model<-tryCatch(SurvModel_xgboost(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting xgboost")
      return()
    })
  }
  if ("gam"%in% models){
    gam_Model<-tryCatch(tryCatch(SurvModel_GAM(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){SurvModel_GAM(data=datatrainFact,expvars=ExpVars2, timevar=timevar, eventvar=eventvar)}), error=function(e){
      print("Failed fitting GAM")
      return()
    })
  }
  if ("gbm"%in% models){
    gbm_Model<-tryCatch(SurvModel_gbm(data=datatrainFact,expvars=ExpVars, timevar=timevar, eventvar=eventvar), error=function(e){
      print("Failed fitting gbm")
      return()
    })
  }
  if ("deepsurv" %in% models) {
    ExpVarsW<-ExpVars2
    while (TRUE) {
      deepsurv_Model <- tryCatch(
        SurvModel_deepsurv(data = datatrainFact, expvars = ExpVarsW, timevar = timevar, eventvar = eventvar, ...),
        error = function(e) {
          print("Failed fitting DeepSurv")
          return()
        }
      )

      # Check if the model fitting was successful
      if (!inherits(deepsurv_Model, "try-error")) {
        pred <- tryCatch(
          Predict_SurvModel_deepsurv(deepsurv_Model, datatrainFact),
          error = function(e) {
            print("Failed making predictions")
            return()
          }
        )

        # Check if predictions were successful and if probs contain NA values
        if (!inherits(pred, "try-error") && !anyNA(pred$Probs)) {
          probs <- pred$Probs
          break  # Exit the loop since successful prediction was made without NA values
        }
      }

      if (length(ExpVarsW) > 1) {
        ExpVarsW <- ExpVarsW[-length(ExpVarsW)]  # Drop the last variable from ExpVars2
      } else {
        break  # Exit the loop if no more variables to drop
      }
    }
  }

  input<-list(ExpVars=ExpVars,ExpVars2=ExpVars2, timevar=timevar, eventvar=eventvar)
  input2<-vector(mode="list")
  input2<-c(input2,list(RF_Model=RF_Model,RF_Model2=RF_Model2))

  if ("glmnet" %in% models){
    input2<-c(input2,list(glmnet_Model=glmnet_Model))
  }
  if ("coxph" %in% models){
    input2<-c(input2,list(CPH_Model=CPH_Model))
  }
  if ("rulefit"%in% models){
    input2<-c(input2,list(RuleFit_Model=RuleFit_Model))
  }
  if ("xgboost"%in% models){
    input2<-c(input2,list(xgboost_Model=xgboost_Model))
  }
  if ("deepsurv"%in% models){
    input2<-c(input2,list(deepsurv_Model=deepsurv_Model))
  }
  if ("gam"%in% models){
    input2<-c(input2,list(gam_Model=gam_Model))
  }
  if ("gbm"%in% models){
    input2<-c(input2,list(gbm_Model=gbm_Model))
  }
  if ("bart"%in% models){
  input2<-c(input2,list(bart_Model=bartout))
  }
  if ("ExpSurvReg"%in% models){
    input2<-c(input2,list(survregexp_Model=survregexp_Model))
  }
  if ("WeibSurvReg"%in% models){
    input2<-c(input2,list(survregweib_Model=survregweib_Model))
  }
  return(c(list(input=input),input2))
}



#' @title PredictSurvModels
#'
#' @description Utility function to get predictions from  several models
#' at the same time. Also adds an ensemble prediction for the probabilities.
#' @param models the output from 'RunSurvModels'
#' @param newdata the data for which the predictions are to be calculated
#' @param newtimes the times for which the predictions obtained
#' @return a list containing the following objects.
#' ModelPredictions: a list of individual model predictions,
#' NewProbs: ensembled predictions
#'
#' @export
PredictSurvModels<-function(models, newdata, newtimes){
  newdataFactor<-newdata
  for (i in 1:ncol(newdata)){
    if (is.factor(newdata[,i])){
      newdata[,i]<-as.character(newdata[,i])
    }
  }
  if (!is.null(models$RuleFit_Model)){
    Predict_rulefit<-Predict_SurvModel_rulefit(models$RuleFit_Model, newdata=newdataFactor)
  }

  if (!is.null(models$RF_Model)){
    Predict_RF<-Predict_SurvModel_RF(models$RF_Model, newdata=newdataFactor)
  }
  if (!is.null(models$RF_Model2)){
    Predict_RF2<-Predict_SurvModel_RF(models$RF_Model2, newdata=newdataFactor)
  }
  if (!is.null(models$glmnet_Model)){
    Predict_glmnet<-Predict_SurvModel_glmnet(models$glmnet_Model, newdata=newdataFactor)
  }
  if (!is.null(models$CPH_Model)){
    Predict_CPH<-Predict_SurvModel_Cox(models$CPH_Model, newdata=newdataFactor)
  }
  if (!is.null(models$bart_Model)) {
    Predict_bart_Model<-Predict_SurvModel_BART(models$bart_Model, newdata=newdataFactor)
  }
  if (!is.null(models$deepsurv_Model)){
    Predict_deepsurv_Model<-Predict_SurvModel_deepsurv(models$deepsurv_Model, newdata=newdataFactor)
  }
  if (!is.null(models$gam_Model)){
    Predict_gam_Model<-Predict_SurvModel_GAM(models$gam_Model, newdata=newdataFactor)
  }
  if (!is.null(models$gbm_Model)){
    Predict_gbm_Model<-Predict_SurvModel_gbm(models$gbm_Model, newdata=newdataFactor)
  }
  if (!is.null(models$survregexp_Model)){
    Predict_survregexp_Model<-Predict_SurvModel_SurvReg(models$survregexp_Model, newdata=newdataFactor)
  }
  if (!is.null(models$survregweib_Model)){
    Predict_survregweib_Model<-Predict_SurvModel_SurvReg(models$survregweib_Model, newdata=newdataFactor)
    }
  if (!is.null(models$xgboost_Model)){
    Predict_xgboost_Model<-Predict_SurvModel_xgboost(models$xgboost_Model, newdata=newdataFactor)
    }
#########################
  if (!is.null(models$RuleFit_Model)){
    newprobsrulefit<-survprobMatInterpolator(probsMat=t(Predict_rulefit$Probs),times=Predict_rulefit$Times, newtimes=newtimes)
    }
  if (!is.null(models$RF_Model)){
    newprobsRF<-survprobMatInterpolator(probsMat=t(Predict_RF$Probs),times=Predict_RF$Times, newtimes=newtimes)
  }

  if (!is.null(models$RF_Model2)){
    newprobsRF2<-survprobMatInterpolator(probsMat=t(Predict_RF2$Probs),times=Predict_RF2$Times, newtimes=newtimes)
  }

  if (!is.null(models$glmnet_Model)){
    newprobsglmnet<-survprobMatInterpolator(probsMat=t(Predict_glmnet$Probs),times=Predict_glmnet$Times, newtimes=newtimes)
  }

  if (!is.null(models$CPH_Model)){
    newprobscph<-survprobMatInterpolator(probsMat=t(Predict_CPH$Probs),times=Predict_CPH$Times, newtimes=newtimes)
  }

   if (!is.null(models$xgboost_Model)){
     newprobsxgboost<-survprobMatInterpolator(probsMat=t(Predict_xgboost_Model$Probs),times=Predict_xgboost_Model$Times, newtimes=newtimes)
   }

  if (!is.null(models$deepsurv_Model)){

    newprobsdeepsurv_Model<-survprobMatInterpolator(probsMat=t(Predict_deepsurv_Model$Probs),times=Predict_deepsurv_Model$Times, newtimes=newtimes)
  }
  if (!is.null(models$gam_Model)){

    newprobsgam_Model<-survprobMatInterpolator(probsMat=t(Predict_gam_Model$Probs),times=Predict_gam_Model$Times, newtimes=newtimes)
    }

  if (!is.null(models$gbm_Model)){

    newprobsgbm_Model<-survprobMatInterpolator(probsMat=t(Predict_gbm_Model$Probs),times=Predict_gbm_Model$Times, newtimes=newtimes)
    }


  if (!is.null(models$survregexp_Model)){
    newprobssurvregexp<-survprobMatInterpolator(probsMat=t(Predict_survregexp_Model$Probs),times=Predict_survregexp_Model$Times, newtimes=newtimes)
    }


  if (!is.null(models$survregweib_Model)){
    newprobssurvregweib<-survprobMatInterpolator(probsMat=t(Predict_survregweib_Model$Probs),times=Predict_survregweib_Model$Times, newtimes=newtimes)
    }


  if (!is.null(models$bart_Model)){

    newprobsbart_Model<-survprobMatInterpolator(probsMat=t(Predict_bart_Model$Probs),times=Predict_bart_Model$Times, newtimes=newtimes)
    }


  ModelPredictions<-list()
  if (!is.null(models$RuleFit_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsrulefit=newprobsrulefit))

  }

  if (!is.null(models$glmnet_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsglmnet=newprobsglmnet))

  }
  if (!is.null(models$RF_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsRF=newprobsRF))
    }


  if (!is.null(models$RF_Model2)){
    ModelPredictions<-c(ModelPredictions,list(newprobsRF2=newprobsRF2))
    }

  if (!is.null(models$CPH_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobscph=newprobscph))
    }

  if (!is.null(models$xgboost_Model)){
     ModelPredictions<-c(ModelPredictions,list(newprobsxgboost=newprobsxgboost))
    }
  if (!is.null(models$deepsurv_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsdeepsurv_Model=newprobsdeepsurv_Model))
    }
  if (!is.null(models$gam_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsgam_Model=newprobsgam_Model))
  }
  if (!is.null(models$gbm_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsgbm_Model=newprobsgbm_Model))
  }

  if (!is.null(models$survregexp_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobssurvregexp_Model=newprobssurvregexp))
}
  if (!is.null(models$survregweib_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobssurvregweib_Model=newprobssurvregweib ))
}
  if (!is.null(models$bart_Model)){
    ModelPredictions<-c(ModelPredictions,list(newprobsbart_Model=newprobsbart_Model))
  }
  NewProbs<-survprobMatListAveraging(ModelPredictions)
  NewProbs<-NewProbs
  return(list(ModelPredictions=ModelPredictions,NewProbs=NewProbs))
}
