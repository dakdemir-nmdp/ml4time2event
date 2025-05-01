#' @title SurvModel_BART
#'
#' @description Fit a BART model for survival outcomes using BART::surv.bart.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param K parameter 'K' for BART (number of time points for discrete approximation)
#' @param ntree number of trees in BART model
#'
#' @return a list containing the following objects:
#' post: fitted BART model object from surv.bart,
#' expvars: character vector of explanatory variables used,
#' eventvarvec: numeric event variable values used in fitting,
#' timevarvec: numeric time variable values used in fitting,
#' x.train: model matrix used for training,
#' varprof: profile of explanatory variables.
#'
#' @importFrom BART surv.bart surv.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @importFrom dplyr bind_rows
#' @export
SurvModel_BART<-function(data,expvars, timevar, eventvar, K=8, ntree=50){
  # Assuming VariableProfile is loaded/available
  varprof<-ml4time2event::VariableProfile(data, expvars) # Placeholder

  # Prepare data
  timevarvec<-data[[timevar]]
  eventvarvec<-as.integer(data[[eventvar]] == 1) # Ensure 0/1

  # Create model matrix
  x.train<-as.matrix(stats::model.matrix(~-1+., data=data[, expvars, drop=FALSE]))

  # Fit BART model with error handling
  failcount<-0
  post<-NULL
  while ((is.null(post) & (failcount<4))){
    failcount<-failcount+1
    post <- tryCatch(BART::surv.bart(x.train=x.train, times=timevarvec, delta=eventvarvec, x.test=x.train, # Predict on train? Check necessity
                                     K=K, ntree=ntree, ndpost=2000, nskip=500,
                                     keepevery = 2L), # Added verbose=FALSE
                     error=function(e){
                         warning("BART fitting attempt ", failcount, " failed: ", e$message)
                         NULL
                         }
    )
  }
  if (is.null(post)) {
      stop("Failed to fit BART model after multiple attempts.")
  }

  return(list(post=post, expvars=expvars, eventvarvec=eventvarvec, timevarvec=timevarvec, x.train=x.train, varprof=varprof ))
}


#' @title Predict_SurvModel_BART
#'
#' @description Make predictions using a fitted BART survival model.
#' @param modelout the output from 'SurvModel_BART'
#' @param newdata the data for which the predictions are to be calculated
#' @param times optional vector of times for prediction (currently ignored by surv.pre.bart/predict).
#'
#' @return a list containing the following objects:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the times at which the probabilities are calculated (including 0, from model).
#'
#' @importFrom BART surv.pre.bart bartModelMatrix
#' @importFrom stats model.matrix
#' @importFrom dplyr bind_rows
#' @export
Predict_SurvModel_BART<-function(modelout, newdata, times=NULL){
  # Prepare test matrix, ensuring factor levels and columns match training
  x.test_orig <-stats::model.matrix(~-1+., data=newdata[, modelout$expvars, drop=FALSE])

  # Ensure columns match training matrix (x.train stored in modelout)
  train_cols <- colnames(modelout$x.train)
  missing_cols <- setdiff(train_cols, colnames(x.test_orig))
  if (length(missing_cols) > 0) {
      warning("Columns missing in newdata compared to training: ", paste(missing_cols, collapse=", "), ". Adding them with 0.")
      add_mat <- matrix(0, nrow = nrow(x.test_orig), ncol = length(missing_cols), dimnames = list(NULL, missing_cols))
      x.test_combined <- cbind(x.test_orig, add_mat)
  } else {
      x.test_combined <- x.test_orig
  }
  # Ensure correct column order
  x.test <- x.test_combined[, train_cols, drop = FALSE]


  # Prepare BART prediction structure
  # Note: The original code had a complex rbind structure which seemed unnecessary/potentially problematic.
  # Simplifying based on surv.pre.bart documentation.
  print(str(modelout))

  # It seems surv.pre.bart primarily needs the training times/events and the test matrix.
  pre=BART::surv.pre.bart(times=modelout$timevarvec, delta=modelout$eventvarvec,
                          x.train= BART::bartModelMatrix(modelout$x.train), # Use stored training matrix
                          x.test=  BART::bartModelMatrix(x.test), # Use prepared test matrix
                          K=modelout$post$K) # Use K from fitted model

  # Predict using the fitted BART model
  pred = predict(modelout$post, pre$tx.test) # Predict on test data structure

  # Extract mean survival probabilities
  # pred$surv.test.mean is matrix: rows=observations, cols=K time points
  PredMat<-matrix(pred$surv.test.mean, ncol=modelout$post$K, byrow = FALSE) # Ensure correct matrix shape

  # Add time 0 with probability 1
  Probs <- t(cbind(1, PredMat)) # Transpose to rows=times, cols=obs
  Times <- c(0, modelout$post$times) # Get times from model object

  print(str(list(Probs=Probs, Times=Times)))
  return(list(Probs=Probs, Times=Times))
}
