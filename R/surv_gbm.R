#' @title SurvModel_gbm
#'
#' @description Fit a gbm model for survival outcomes using gbm package.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param ntree number of trees to grow
#' @param max.depth maximum depth for the trees
#' @param bag.fraction fraction of data to sample in each bagging iteration
#' @param train.fraction fraction for training data (used internally by gbm for CV/OOB error)
#' @param learninrate learning rate (shrinkage) for the boosting algorithm
#'
#' @return a list which contains the following objects:
#' gbmmodel: fitted gbm model object,
#' best.iter: best iteration number based on cross-validation performance,
#' expvars: vector of the names of explanatory variables used,
#' survMat: survival probabilities matrix for training data (rows=obs, cols=times),
#' basehaz.cum: baseline cumulative hazard estimated on training data,
#' time.interest: unique observed event times from training data,
#' varprof: profile of explanatory variables.
#'
#' @importFrom gbm gbm gbm.perf basehaz.gbm
#' @importFrom survival Surv
#' @export
SurvModel_gbm<-function(data,expvars, timevar, eventvar, ntree=200, max.depth=3, bag.fraction=.3, train.fraction=.3, learninrate=.01){
  # Assuming VariableProfile is loaded/available
  varprof<-VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[,eventvar]<-as.numeric(data[,eventvar]==1)

  # Convert character columns to factors for gbm
  fctvars<-names(which(sapply(as.data.frame(data[,expvars, drop=FALSE]), function(x){(is.factor(x) | is.character(x))})))
  for (i in fctvars){
    if (is.character(data[[i]])) {
        data[[i]]<-as.factor(data[[i]])
    }
  }

  # Define formula
  formula_gbm <- survival::Surv(data[[timevar]], data[[eventvar]]) ~ .

  # Fit gbm model
  gbmmodel <- gbm::gbm(formula = formula_gbm,
                   data=data[, expvars, drop=FALSE], # Use only predictors here
                   distribution="coxph",
                   n.trees=ntree,
                   shrinkage=learninrate,
                   interaction.depth=max.depth,
                   bag.fraction = bag.fraction,
                   train.fraction = train.fraction, # For internal CV/OOB error estimation
                   cv.folds = 2, # Use 2 folds as in original code
                   # n.minobsinnode = 10, # Use gbm default
                   keep.data = TRUE, # Keep data for prediction/baseline hazard
                   verbose = FALSE) # Suppress verbose output

  # Find best iteration based on CV performance
  best.iter <- gbm::gbm.perf(gbmmodel, method="cv", plot.it = FALSE)

  # Get unique event times from training data
  time.interest <- sort(unique(data[[timevar]][data[[eventvar]]==1]))

  # Predict linear predictor on training data
  pred.train <- predict(gbmmodel, data[, expvars, drop=FALSE], n.trees = best.iter, type="link")

  # Estimate the cumulative baseline hazard function using training data
  basehaz.cum <- gbm::basehaz.gbm(t = data[[timevar]], delta = data[[eventvar]], f.x = pred.train, t.eval = time.interest, cumulative = TRUE)

  # Calculate survival probabilities for training data (optional, but done in original)
  survMat<-NULL
  for (i in 1:length(pred.train)){
    surf.i <- exp(-exp(pred.train[i])*basehaz.cum)
    survMat<-rbind(survMat,surf.i) # Matrix: rows=obs, cols=times
  }

  return(list(gbmmodel=gbmmodel, best.iter=best.iter, expvars=expvars,
              survMat=survMat, basehaz.cum=basehaz.cum, time.interest=time.interest,
              varprof=varprof))
}




#' @title Predict_SurvModel_gbm
#'
#' @description Get predictions from a gbm survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_gbm'
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the unique times for which the probabilities are calculated (including 0).
#'
#' @export
Predict_SurvModel_gbm<-function(modelout, newdata){
  # Prepare newdata: ensure factors have same levels as training data
  data_test <- newdata[, modelout$expvars, drop=FALSE]
  for (vari in modelout$expvars){
    if (is.factor(modelout$gbmmodel$data$x[[vari]])) { # Check if var was factor in training
        train_levels <- levels(modelout$gbmmodel$data$x[[vari]])
        # Ensure the column exists in newdata before attempting to modify
        if (vari %in% colnames(data_test)) {
            # Convert to character first to handle potential new levels, then factor
            data_test[[vari]] <- factor(as.character(data_test[[vari]]), levels = train_levels)
        }
    } else if (is.character(modelout$gbmmodel$data$x[[vari]])) {
         # If original was character but became factor in gbm, treat as factor
         # This case might need careful handling depending on gbm's internal factor conversion
         # Assuming it was treated as factor if is.factor check above failed but it's character
         # This part might need refinement based on how gbm handles character inputs
         if (vari %in% colnames(data_test) && is.character(data_test[[vari]])) {
             # Attempt to find levels if possible (e.g., from varprof if stored)
             # Placeholder: Assume direct conversion works or levels are handled
             # train_levels <- names(modelout$varprof[[vari]]) # Example if levels stored in varprof
             # data_test[[vari]] <- factor(data_test[[vari]], levels = train_levels)
             # Fallback: Convert to factor directly (might cause errors if new levels)
             data_test[[vari]] <- as.factor(data_test[[vari]])
         }
    }
  }


  # Predict linear predictor for newdata
  event_prediction <- suppressWarnings(predict(modelout$gbmmodel, data_test, n.trees=modelout$best.iter, type="link"))

  # Calculate survival probabilities using stored baseline hazard
  survMat<-NULL
  for (i in 1:length(event_prediction)){
    surf.i <- exp(-exp(event_prediction[i])*modelout$basehaz.cum)
    survMat<-rbind(survMat,surf.i) # Matrix: rows=obs, cols=times
  }

  # Add time 0 with probability 1
  Probs <- t(cbind(1, survMat)) # Transpose to rows=times, cols=obs
  Times <- c(0, modelout$time.interest)

  return(list(Probs=Probs, Times=Times))
}
