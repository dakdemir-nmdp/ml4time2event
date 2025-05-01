#' @title SurvModel_glmnet
#' @description Fit a penalized Cox model for survival outcomes using glmnet.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param alpha elastic net mixing parameter (1=lasso, 0=ridge)
#' @param maxit maximum number of iterations for glmnet
#' @param nfolds number of folds for cross-validation in cv.glmnet
#'
#' @return a list containing the following objects:
#' traindata: a sample of original training predictors (for factor level consistency),
#' TrainMat: a sample of the model matrix used for training (for prediction structure),
#' yTrain: a sample of the Surv object used for training (for prediction structure),
#' expvars: character vector of explanatory variables used,
#' cv.fit: cv.glmnet fit object,
#' est.coef: estimated model coefficients at lambda.min,
#' varprof: profile of explanatory variables.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats model.matrix coef
#' @importFrom survival Surv survfit
#' @export
SurvModel_glmnet<-function(data,expvars, timevar, eventvar, alpha=.5, maxit=5000, nfolds=30){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Create model matrix (handle factors appropriately)
  TrainMat<-stats::model.matrix(~., data=data[,expvars, drop=FALSE])
  # Remove intercept column if present (often added by model.matrix)
  intercept_col <- which(colnames(TrainMat) == "(Intercept)")
  if (length(intercept_col) > 0) {
      TrainMat <- TrainMat[, -intercept_col, drop = FALSE]
  }


  # Create survival object
  yTrain = survival::Surv(data[[timevar]], data[[eventvar]])

  # Fit cv.glmnet for Cox model
  cv.fit = glmnet::cv.glmnet(x = TrainMat,
                     y = yTrain,
                     alpha = alpha,
                     family = "cox",
                     maxit = maxit,
                     nfolds=nfolds
                     )
  # Get coefficients at lambda.min
  est.coef = stats::coef(cv.fit, s = cv.fit$lambda.min)

  # sfitTrain<-survival::survfit(cv.fit, s =cv.fit$lambda.min, x = TrainMat, y=yTrain, newx=TrainMat) # Not returned
  # predSurvsTrain<-sfitTrain$surv # Not returned
  # lpTrain<-TrainMat%*%est.coef # Not returned
  # timesTrain<-sfitTrain$time # Not returned

  # Store samples for prediction consistency
  sample_rows<-sample(1:nrow(data),min(c(500,nrow(data))))
  out<-list(traindata=data[sample_rows, expvars, drop=FALSE], # Sample of original predictors
            TrainMat=TrainMat[sample_rows, , drop=FALSE], # Sample of model matrix
            yTrain=yTrain[sample_rows, ], # Sample of Surv object
            expvars=expvars,
            cv.fit=cv.fit,
            est.coef=est.coef,
            varprof=varprof)
  return(out)
}



#' @title Predict_SurvModel_glmnet
#'
#' @description Get predictions from a glmnet survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_glmnet'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @importFrom stats model.matrix
#' @importFrom survival survfit
#' @export
Predict_SurvModel_glmnet<-function(modelout, newdata){
  # Create test matrix ensuring factor levels and columns match training
  # Use the rbind trick with the sampled training data
  mmdata<-rbind(modelout$traindata, newdata[,modelout$expvars, drop=FALSE])
  TestMat_full<-stats::model.matrix(~., data=mmdata)
  # Remove intercept if present
  intercept_col <- which(colnames(TestMat_full) == "(Intercept)")
  if (length(intercept_col) > 0) {
      TestMat_full <- TestMat_full[, -intercept_col, drop = FALSE]
  }
  # Select rows corresponding to newdata
  TestMat <- TestMat_full[-c(1:nrow(modelout$traindata)), , drop=FALSE]

  # Ensure TestMat has the same columns as the matrix used for training cv.glmnet
  train_cols <- colnames(modelout$TrainMat) # Use colnames from the saved TrainMat sample
  missing_cols <- setdiff(train_cols, colnames(TestMat))
  for(col in missing_cols){
      TestMat[[col]] <- 0 # Add missing columns with 0 (e.g., factor levels not in newdata)
  }
  # Ensure correct column order and selection
  TestMat <- TestMat[, train_cols, drop = FALSE]


  # Predict survival curves using the fitted model
  sfitTest<-survival::survfit(modelout$cv.fit, s=modelout$cv.fit$lambda.min,
                              x = modelout$TrainMat, # Provide training matrix sample
                              y = modelout$yTrain,   # Provide training Surv object sample
                              newx=TestMat)         # Provide new data matrix
  predSurvsTest<-sfitTest$surv # Matrix: rows=times, cols=observations
  timesTest<-sfitTest$time

  # Add time 0 with probability 1 if missing
  if (sum(timesTest==0)==0){
    timesTest<-c(0,timesTest)
    predSurvsTest<-rbind(rep(1, ncol(predSurvsTest)),predSurvsTest)
  }
  return(list(Probs=predSurvsTest, Times=timesTest))
}
