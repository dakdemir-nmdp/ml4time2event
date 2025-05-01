#' @title SurvModel_RF
#'
#' @description Fit a RF model for survival outcomes using randomForestSRC.
#' Includes tuning of nodesize and mtry.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param ntree integer value,  number of trees to grow
#' @param samplesize integer value,  sample size for each grown tree (swor)
#' @param nsplit integer value, maximum number of splits for each tree
#' @param trace logical, trace tuning process or not
#'
#' @return a list containing the fitted randomForestSRC object ('hd.obj') and
#' variable profile ('varprof').
#'
#' @importFrom randomForestSRC tune rfsrc predict.rfsrc
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @export
SurvModel_RF<-function(data,expvars, timevar, eventvar, ntree=300, samplesize=500, nsplit=5, trace=TRUE){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Convert character columns to factors for randomForestSRC
  for (vari in expvars){
    if (is.character(data[[vari]])){
      data[[vari]]<-as.factor(data[[vari]])
    }
  }

  # Define formula
  formRF<-stats::as.formula(paste("Surv(",timevar, ",", eventvar,") ~ .", collapse = "")) # Removed survival::

  # Adjust samplesize if it exceeds 70% of data
  samplesize <- min(ceiling(0.7 * nrow(data)), samplesize)

  # Tune hyperparameters (nodesize, mtry)
  # Using bs.gradient split rule and sampling without replacement (swor) as in original code
  o <- randomForestSRC::tune(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                             splitrule="bs.gradient", samptype = "swor", sampsize = samplesize,
                             trace = trace, nsplit=nsplit, stepFactor = 1.5,
                             mtryStart = 2, # Start tuning mtry from 2
                             nodesizeTry = c(seq(1, 101, by = 10)), # Tune nodesize
                             ntreeTry = ntree) # Use fixed ntree for tuning speed

  # Fit final model with optimal parameters
  hd.obj <- randomForestSRC::rfsrc(formRF, data = data[,c(timevar, eventvar, expvars), drop=FALSE],
                                   nodesize =o$optimal[[1]], ntree=ntree, mtry= o$optimal[[2]],
                                   tree.err = FALSE, importance = TRUE, statistics=TRUE,
                                   do.trace = trace, splitrule="bs.gradient", samptype = "swor",
                                   sampsize = samplesize, nsplit = nsplit)

  return(list(hd.obj=hd.obj, varprof=varprof))
}

#' @title Predict_SurvModel_RF
#'
#' @description Get predictions from a RF survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_RF' (a list containing 'hd.obj')
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#'  Probs: predicted Survival probability matrix (rows=times, cols=observations),
#'  Times: The times at which the probabilities are predicted,
#'  predSurvsTestRF: the raw output object from 'randomForestSRC::predict.rfsrc'.
#'
#' @importFrom randomForestSRC predict.rfsrc
#' @export
Predict_SurvModel_RF<-function(modelout, newdata){
  # Ensure character columns in newdata are factors with same levels as training data
  # This requires access to training data factor levels, often stored in varprof or model object
  # Placeholder: Assuming levels are handled correctly or newdata is preprocessed.
  # A more robust solution would involve storing factor levels from training.
  for (vari in names(modelout$varprof)) {
      if (inherits(modelout$varprof[[vari]], "table")) { # Check if it was a factor/character
          if (vari %in% colnames(newdata) && is.character(newdata[[vari]])) {
              # Attempt to apply levels from training profile
              training_levels <- names(modelout$varprof[[vari]])
              newdata[[vari]] <- factor(newdata[[vari]], levels = training_levels)
          }
      }
  }


  predSurvsTestRF<-randomForestSRC::predict.rfsrc(modelout$hd.obj, newdata = newdata)

  # Extract survival probabilities and times
  # Add time 0 with probability 1
  Probs<-rbind(1, predSurvsTestRF$survival) # rfsrc survival matrix is obs x times
  Times<-c(0, predSurvsTestRF$time.interest)

  return(list(
    Probs = Probs, # Return as rows=times, cols=observations
    Times = Times,
    predSurvsTestRF=predSurvsTestRF) # Return the full prediction object
    )
}
