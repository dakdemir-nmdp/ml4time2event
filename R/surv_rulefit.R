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
#' varprof: profile of explanatory variables,
#' ruleslist: list of the fitted rules,
#' TrainMat: fitted rules matrix (combined with original features),
#' ctreelist: list of tree models,
#' yTrain: survival object for outcome variables,
#' cv.fitRules: fitted glmnet survival models for a range of shrinkage parameters,
#' timesTrain: unique times in the outcome variable from survfit object
#' expvars: explanatory variables used in the model
#'
#' @importFrom rpart rpart rpart.control
#' @importFrom partykit as.party
#' @importFrom glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix quantile rpois runif
#' @importFrom survival Surv survfit
#' @export
SurvModel_rulefit <- function(data,
                              expvars,
                              timevar,
                              eventvar,
                              ntree = 300,
                              nsample = 300,
                              keepvars = NULL,
                              cuttimes= NULL) {
  # Assuming VariableProfile and listrules are loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder if not loaded

  if (is.null(cuttimes)){
    x<-data[,timevar]
    cuttimes<-stats::quantile(x, c(.1,.25,.50,.70,.90))
  }
  formRF1 <-
    stats::as.formula(paste("survival::Surv(", timevar, ",", eventvar, ") ~.", collapse = ""))
  formRF2 <- stats::as.formula(paste("ClassVar ~.", collapse = ""))

  ctreelist <- lapply(1:ntree, function(repi) {
    sampcols <-
      union(keepvars, sample(expvars, min(length(expvars), sample(c(
        1:10
      ), 1))))
    selmodel <- sample(c(1, 1, 1, 1, 1, 1, 2), 1) # Favor survival trees
    formRF <- list(formRF1, formRF2)[[selmodel]]
    if (selmodel == 1) { # Survival Tree
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
      rpartmodel <- rpart::rpart(formRF,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    } else { # Classification Tree
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
      rpartmodel <- rpart::rpart(formRF,
                                 data = datasampl,
                                 control = rpcontrol,
                                 model = TRUE)
    }
    return(rpartmodel)
  })

  # Extract rules
  ruleslist <- lapply(ctreelist, function(x) {
      tryCatch(listrules(partykit::as.party(x)), error = function(e) list()) # Placeholder
  })

  # Create rule matrix for training data
  RulesTrain <-
    lapply(ctreelist, function(x)
      predict(partykit::as.party(x), newdata = data, type = "node"))
  RulesTrainMatList <- lapply(1:length(RulesTrain), function(x) {
    Train <- RulesTrain[[x]]
    rules <- ruleslist[[x]]
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
        Train <- factor(as.character(Train), levels = names(rules))
        if (length(unique(names(rules))) > 1) {
          MM <- stats::model.matrix(~ -1 + Train)
          colnames(MM) <-
            paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
          MM
        } else { NULL }
    } else { NULL }
  })
  RulesTrainMatList <- Filter(Negate(is.null), RulesTrainMatList)

  # Combine original features and rules
  TrainMatOrig <- stats::model.matrix(~.-1, data=data[,expvars])
  if (length(RulesTrainMatList) > 0) {
      TrainMatRules <- Reduce("cbind", RulesTrainMatList)
      TrainMat <- cbind(TrainMatOrig, TrainMatRules)
  } else {
      TrainMat <- TrainMatOrig
  }

  # Prepare outcome for glmnet
  yTrain = survival::Surv(data[, colnames(data) == timevar], as.numeric(data[, colnames(data) == eventvar] == 1))

  # Fit penalized Cox model (glmnet)
  cv.fitRules = glmnet::cv.glmnet(
      x = TrainMat,
      y = yTrain,
      alpha = .5, # Elastic Net
      family = "cox",
      maxit = 1000
    )
  # est.coefRules = stats::coef(cv.fitRules, s = cv.fitRules$lambda.1se) # Not returned but could be useful

  # Get survival fit for training data (needed for prediction times)
  sfitTrainRules <-
      survival::survfit(
        cv.fitRules,
        s = cv.fitRules$lambda.1se, # Or lambda.min
        x = TrainMat,
        y = yTrain,
        newx = TrainMat # Predict on training data itself
      )
  timesTrain <- sfitTrainRules$time

  # Sample data for prediction consistency
  nrowforsample<-min(c(nrow(data),500))

  return(
      list(
        datasamp=data[1:nrowforsample,expvars], # Sample of original predictors
        varprof=varprof,
        ruleslist = ruleslist,
        TrainMat = TrainMat[1:nrowforsample,], # Sample of combined matrix (for prediction structure)
        ctreelist = ctreelist,
        yTrain = yTrain[1:nrowforsample,], # Sample of outcome (for prediction structure)
        cv.fitRules = cv.fitRules,
        timesTrain = timesTrain,
        expvars=expvars # Original predictor names
        )
    )
}


#' @title Predict_SurvModel_rulefit
#'
#' @description Get predictions from a rulefit model for a test dataset
#'
#' @param modelout the output from 'SurvModel_rulefit'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list object with three items. Probs: survival probability predictions matrix (rows=times, cols=observations),
#' Times: the times of the returned probabilities, sfitTestRules: prediction object from survfit
#'
#' @importFrom partykit as.party
#' @importFrom stats model.matrix
#' @importFrom survival survfit
#' @export
Predict_SurvModel_rulefit <- function(modelout, newdata) {

  # Create rule matrix for test data
  RulesTest <-
    lapply(modelout$ctreelist, function(x)
      predict(partykit::as.party(x), newdata = newdata, type = "node"))
  RulesTestMatList <- lapply(1:length(RulesTest), function(x) {
    Test <- RulesTest[[x]]
    rules <- modelout$ruleslist[[x]]
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
        Test <- factor(as.character(Test), levels = names(rules))
        if (length(unique(names(rules))) > 1) {
          MM <- stats::model.matrix(~ -1 + Test)
          colnames(MM) <-
            paste(paste(x, 1:length(rules), sep = "_"), rules, sep = "***")
          MM
        } else { NULL }
    } else { NULL }
  })
  RulesTestMatList <- Filter(Negate(is.null), RulesTestMatList)

  # Combine original features and rules for test data
  # Ensure factor levels match training sample using rbind trick
  MM_orig_test <- stats::model.matrix(~.-1, data=rbind(modelout$datasamp, newdata[,modelout$expvars]))
  MM_orig_test <- MM_orig_test[-c(1:nrow(modelout$datasamp)), , drop = FALSE]

  if (length(RulesTestMatList) > 0) {
      TestMatRules <- Reduce("cbind", RulesTestMatList)
      TestMat <- cbind(MM_orig_test, TestMatRules)
  } else {
      TestMat <- MM_orig_test
  }

  # Ensure TestMat has the same columns as the training matrix used in cv.glmnet
  train_cols <- colnames(modelout$TrainMat) # Use colnames from the saved TrainMat sample
  missing_cols <- setdiff(train_cols, colnames(TestMat))
  for(col in missing_cols){
      TestMat[[col]] <- 0 # Add missing rule columns with 0
  }
  # Ensure correct column order and selection
  TestMat <- TestMat[, train_cols, drop = FALSE]


  # Predict using the fitted glmnet model
  sfitTestRules <-
      survival::survfit(
        modelout$cv.fitRules,
        s = modelout$cv.fitRules$lambda.min, # Use lambda.min for prediction
        x = modelout$TrainMat, # Provide the training matrix used by cv.glmnet
        y = modelout$yTrain,   # Provide the training outcome used by cv.glmnet
        newx = TestMat # Provide the new data matrix
      )
  predSurvsTestRules <- sfitTestRules$surv
  timesTest <- sfitTestRules$time

  # Add time 0 with probability 1 if missing
  if (sum(timesTest==0)==0){
      timesTest<-c(0,timesTest)
      predSurvsTestRules<-rbind(rep(1, ncol(predSurvsTestRules)),predSurvsTestRules)
  }
  return(
      list(
        Probs = predSurvsTestRules, # Rows=times, Cols=observations
        Times = timesTest,
        sfitTestRules=sfitTestRules
      )
    )
}
