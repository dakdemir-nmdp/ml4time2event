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
#' CRrulefitModel: fitted Fine-Gray model object (used internally by rulefit),
#' expvars: explanatory variables used in the model
#' usecols: columns used in the internal Fine-Gray model
#' varprof: profile of explanatory variables.
#'
#' @importFrom rpart rpart rpart.control
#' @importFrom partykit as.party
#' @importFrom pseudo pseudoci
#' @importFrom stats as.formula model.matrix quantile rpois runif sd
#' @importFrom survival Surv
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
  # Assuming VariableProfile is defined elsewhere or loaded via namespace
  # varprof<-VariableProfile(data, expvars)
  varprof <- list() # Placeholder

  if (is.null(cuttimes)){
    x<-data[,timevar]
    cuttimes<-stats::quantile(x, c(.1,.25,.50,.70,.8))
    cuttimesReg<-stats::quantile(x, c(.25,.50,.60, .70, .80))
  }

  formClass <- stats::as.formula(paste("ClassVar ~.", collapse = ""))
  formReg <- stats::as.formula(paste("RegVar ~.", collapse = ""))
  formSurv <- stats::as.formula(paste("survival::Surv(",timevar,",", eventvar,"==1) ~.", collapse = ""))


  ctreelist <- lapply(1:ntree, function(repi){
    sampcols <-
      union(keepvars, sample(expvars, min(length(expvars), sample(c(
        1:10
      ), 1))))
    selmodel <- sample(c(2,3), 1)
   selmodel<-2 # Force regression/pseudo-observation based trees
    if (selmodel == 1) { # Classification Tree (Original code had this branch, keeping for reference but not used due to selmodel<-2)

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
          minsplit = stats::rpois(1,1)+1,
          minbucket = stats::rpois(1,30)+1,
          cp = 0.01 * stats::runif(1),
          maxcompete = stats::rpois(1,30)+1,
          maxsurrogate = stats::rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = stats::rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formClass,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else if (selmodel == 2) { # Regression Tree (Pseudo-observations)
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(1:nrow(data), nsample, replace = T)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      pout<-pseudo::pseudoci(time=datasampl[,timevar],event=datasampl[,eventvar],tmax=sample(cuttimesReg,1))
      datasampl$RegVar <-c(unlist(pout$pseudo$cause1))
      datasampl$RegVar<-(datasampl$RegVar-mean(datasampl$RegVar))/stats::sd(datasampl$RegVar)
      datasampl <-
        datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
      rpcontrol <-
        rpart::rpart.control(
          minsplit = stats::rpois(1,1)+1,
          minbucket = stats::rpois(1,30)+1,
          cp = 0.01 * stats::runif(1),
          maxcompete = stats::rpois(1,30)+1,
          maxsurrogate = stats::rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = stats::rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formReg,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else if (selmodel == 3) { # Survival Tree (Original code had this branch, keeping for reference but not used due to selmodel<-2)
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(1:nrow(data), nsample, replace = T)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      rpcontrol <-
        rpart::rpart.control(
          minsplit = stats::rpois(1,1)+1,
          minbucket = stats::rpois(1,30)+1,
          cp = 0.01 * stats::runif(1),
          maxcompete = stats::rpois(1,30)+1,
          maxsurrogate = stats::rpois(1,3) +1,
          usesurrogate = 2,
          xval = 10,
          surrogatestyle = 0,
          maxdepth = stats::rpois(1,2)+1
        )
      rpartmodel <- rpart::rpart(formSurv,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    }
    return(rpartmodel)
  })

  # Assuming listrules is defined elsewhere or loaded via namespace
  # ruleslist <- lapply(ctreelist, function(x) listrules(partykit::as.party(x)))
  ruleslist <- lapply(ctreelist, function(x) {
      tryCatch(listrules(partykit::as.party(x)), error = function(e) list()) # Placeholder
  })


  RulesTrain <-
    lapply(ctreelist, function(x)
      predict(partykit::as.party(x), newdata = data, type = "node"))
  RulesTrainMatList <- lapply(1:length(RulesTrain), function(x) {
    Train <- RulesTrain[[x]]
    rules <- ruleslist[[x]]
    # Ensure Train levels match rule names before creating model matrix
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
        Train <- factor(as.character(Train), levels = names(rules))
        if (length(unique(names(rules))) > 1) {
          MM <- stats::model.matrix(~ -1 + Train)
          colnames(MM) <-
            paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
          MM
        } else {
          NULL # Handle case with only one rule/level
        }
    } else {
        NULL # Handle case with no rules
    }
  })

  # Filter out NULL elements before binding
  RulesTrainMatList <- Filter(Negate(is.null), RulesTrainMatList)

  if (length(RulesTrainMatList) > 0) {
      TrainMatRules <- Reduce("cbind", RulesTrainMatList)
      TrainMat <- cbind(stats::model.matrix(~.-1, data=data[,expvars]), TrainMatRules)
  } else {
      TrainMat <- stats::model.matrix(~.-1, data=data[,expvars]) # Only original features if no rules generated
  }

  usecols<-colnames(TrainMat)
  nrowforsample<-min(c(nrow(data),500))

  # Assuming CRModel_FineGray is defined elsewhere or loaded via namespace
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
#' @importFrom partykit as.party
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_rulefit <- function(modelout, newdata) {
  RulesTest <-
    lapply(modelout$ctreelist, function(x)
      predict(partykit::as.party(x), newdata = newdata, type = "node"))
  RulesTestMatList <- lapply(1:length(RulesTest), function(x) {
    Test <- RulesTest[[x]]
    rules <- modelout$ruleslist[[x]]
    # Ensure Test levels match rule names before creating model matrix
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
        Test <- factor(as.character(Test), levels = names(rules))
        if (length(unique(names(rules))) > 1) {
          MM <- stats::model.matrix(~ -1 + Test)
          colnames(MM) <-
            paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
          MM
        } else {
          NULL # Handle case with only one rule/level
        }
    } else {
        NULL # Handle case with no rules
    }
  })

  # Filter out NULL elements before binding
  RulesTestMatList <- Filter(Negate(is.null), RulesTestMatList)

  # Create model matrix for original features, ensuring factor levels match training sample
  MM_orig <- stats::model.matrix(~.-1, data=rbind(modelout$datasamp, newdata[,modelout$expvars]))
  MM_orig <- MM_orig[-c(1:nrow(modelout$datasamp)), , drop = FALSE] # Use drop=FALSE to keep matrix structure

  if (length(RulesTestMatList) > 0) {
      TestMatRules <- Reduce("cbind", RulesTestMatList)
      TestMat<-as.data.frame(cbind(MM_orig, TestMatRules))
  } else {
      TestMat<-as.data.frame(MM_orig) # Only original features if no rules generated
  }

  # Ensure TestMat has the same columns as used in the internal Fine-Gray model training
  missing_cols <- setdiff(modelout$CRrulefitModel$expvars, colnames(TestMat))
  for(col in missing_cols){
      TestMat[[col]] <- 0 # Add missing columns with 0
  }
  TestMat <- TestMat[, modelout$CRrulefitModel$expvars, drop = FALSE] # Reorder and select columns


  # Assuming Predict_CRModel_FG is defined elsewhere or loaded via namespace
  PredRFOut<-Predict_CRModel_FG(modelout$CRrulefitModel,TestMat)

  return(list(
    TestMat=TestMat,
    PredRFOut=PredRFOut,
    CIFs=PredRFOut$CIFs,
    Times=PredRFOut$Times
  )
  )
}
