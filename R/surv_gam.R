#' @title score2proba
#' @description internal function: from linear score to survival probabilities using Cox PH baseline.
#' @param datasurv data frame with 'time' and 'event' columns.
#' @param score numeric vector of linear predictor scores.
#' @param conf.int confidence level for survival curve.
#' @param which.est which estimate to return ("point", "lower", "upper").
#' @return a list containing the Cox model ('model') and survfit object ('sf').
#' @importFrom survival coxph Surv survfit coxph.control
#' @noRd
score2proba <-
  function(datasurv, score, conf.int=0.95, which.est=c("point", "lower", "upper")) {
    which.est <- match.arg(which.est)
    # pred <- rep(NA, length(score)) # Not used
    # names(pred) <- names(score) # Not used
    datacox<-cbind(datasurv, data.frame(score=score))
    # Ensure column names are 'time' and 'event' for Surv formula
    colnames(datacox)[1:2]<-c("time","event")
    # Fit a Cox model with score as the only predictor, fixing coefficient to 1 (offset)
    # This estimates the baseline hazard based on the provided scores
    predm <- survival::coxph(survival::Surv(time, event) ~ score, data=datacox, init=1, control=survival::coxph.control(iter.max = 0))
    # Get survival curve predictions based on the fitted baseline hazard and new scores
    sf <- survival::survfit(predm, newdata=data.frame("score"=score), conf.int=conf.int)
    return(list(model=predm,sf=sf))
  }


#' @title SurvModel_GAM
#'
#' @description Fit a GAM model for survival outcomes using Cox PH family.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param shrinkTreshold integer value, minimum number of
#' factor levels for factor variables to be considered for shrinkage ('re' basis).
#'
#' @return a list of four items: model: fitted GAM model object from mgcv::gam,
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables,
#'  expvars: the explanatory variables used.
#'
#' @importFrom mgcv gam cox.ph s
#' @importFrom stats as.formula predict
#' @export
SurvModel_GAM <-
  function(data,
           expvars,
           timevar,
           eventvar,
           shrinkTreshold = 10) {
    # Assuming VariableProfile is loaded/available
    varprof<-VariableProfile(data, expvars) # Placeholder

    # Ensure event variable is numeric 0/1
    data[, eventvar] <- as.numeric(data[, eventvar] == 1)

    # Identify variable types
    numvars <-expvars[which(sapply(as.data.frame(data[, expvars, drop=FALSE]), is.numeric))]
    fctvars <-expvars[which(sapply(as.data.frame(data[, expvars, drop=FALSE]),function(x) {
        (is.factor(x) | is.character(x))}))]

    # Convert characters to factors
    for (i in fctvars) {
      if (is.character(data[[i]])) {
          data[[i]] <- as.factor(data[[i]])
      }
    }

    # Separate variables based on number of levels for smoothing/shrinkage
    catvarstoshrink <-if (length(fctvars)>0){fctvars[which(sapply(as.data.frame(data[, fctvars, drop=FALSE]), function(x) {
        length(levels(x)) > shrinkTreshold # Use levels() for factors
      }))]} else {c()}
    catvarsnottoshrink <-if (length(fctvars)>0){fctvars[which(sapply(as.data.frame(data[, fctvars, drop=FALSE]), function(x) {
      length(levels(x)) <= shrinkTreshold
    }))]} else {c()}
    numvarstosmooth <-if (length(numvars)>0){numvars[which(sapply(as.data.frame(data[, numvars, drop=FALSE]), function(x) {
      length(unique(x[!is.na(x)])) > shrinkTreshold # Use unique non-NA values for numerics
    }))]} else {c()}
    numvarsnottosmooth <-if (length(numvars)>0){numvars[which(sapply(as.data.frame(data[, numvars, drop=FALSE]), function(x) {
      length(unique(x[!is.na(x)])) <= shrinkTreshold
    }))]} else {c()}

    # Build GAM formula string
    formGAM_terms <- c()
    for (vari in catvarstoshrink) {
      formGAM_terms <- c(formGAM_terms, paste0("s(", vari, ", bs='re')")) # Random effect smooth for high-cardinality factors
    }
    for (vari in numvarstosmooth) {
      formGAM_terms <- c(formGAM_terms, paste0("s(", vari, ")")) # Default smooth for numerics
    }
    # Add linear terms for remaining variables
    formGAM_terms <- c(formGAM_terms, numvarsnottosmooth, catvarsnottoshrink)

    # Combine terms into formula
    if (length(formGAM_terms) > 0) {
        formGAM_rhs <- paste(formGAM_terms, collapse = "+")
    } else {
        formGAM_rhs <- "1" # Intercept only if no predictors
    }
    formGAM <- stats::as.formula(paste(timevar, "~", formGAM_rhs))
    print(formGAM) # Print the final formula used

    # Fit the GAM model with Cox PH family
    b <- mgcv::gam(
      formGAM,
      family = mgcv::cox.ph(),
      data = data,
      weights = data[[eventvar]], # Use event indicator as weights
      select = TRUE # Enable shrinkage via double penalty approach
    )

    # Predict linear scores on training data
    fvTrain <- stats::predict(b,
                       newdata = data[, expvars, drop=FALSE],
                       type = "link", # Get linear predictor
                       se.fit = TRUE) # Request standard errors (though not used here)

    # Get survival data (time and event)
    datasurv <- data[, c(timevar, eventvar)]
    # times <- sort(unique(data[data[[eventvar]] == 1, timevar])) # Not directly used here

    # Convert scores to probabilities using baseline hazard from an internal Cox model
    preds <-
      score2proba(
        datasurv = datasurv,
        score = fvTrain$fit,
        conf.int = 0.95,
        which.est = "point"
      )

    # Store the baseline hazard information in the model
    b$baseline_model <- preds$model
    b$baseline_sf <- preds$sf

    # Get unique event times from training data
    times <- sort(unique(data[data[[eventvar]] == 1, timevar]))

    return(list(
      model = b, # The fitted GAM object with baseline hazard info
      times = times,
      varprof = varprof,
      expvars = expvars
    ))
  }

#' @title Predict_SurvModel_GAM
#'
#' @description Get predictions from a GAM survival model for a test dataset.
#'
#' @param modelout the output from 'SurvModel_GAM' (a list containing 'model', 'times', 'varprof', 'expvars')
#' @param newdata the data for which the predictions are to be calculated
#'
#' @return a list containing the following items:
#'  Probs: predicted Survival probability matrix (rows=times, cols=observations),
#'  Times: The times at which the probabilities are predicted.
#'
#' @importFrom stats predict
#' @importFrom survival survfit
#' @export
Predict_SurvModel_GAM <- function(modelout, newdata, newtimes = NULL) {
  if (missing(modelout)) stop("argument \"modelout\" is missing")
  if (missing(newdata)) stop("argument \"newdata\" is missing")

  # Predict linear scores for new data using the fitted GAM
  fvTest <- stats::predict(modelout$model,
                    newdata = newdata,
                    type = "link", # Get linear predictor
                    se.fit = TRUE) # Request standard errors (though not used here)

  # Use the stored baseline hazard from the training fit
  # and the new scores (fvTest$fit) to get survival probabilities for the new data
  sf <-
    survival::survfit(
      modelout$model$baseline_model, # Use the stored Cox model
      newdata = data.frame("score" = fvTest$fit),
      conf.int = .95
    )
  estSURVTest <- sf$surv # Matrix: rows=times, cols=observations
  time.interest <- sf$time

  # Add time 0 with probability 1 if missing
  if (sum(time.interest==0)==0){
    time.interest <- c(0,time.interest)
    estSURVTest <- rbind(rep(1, ncol(estSURVTest)),estSURVTest)
  }

  Probs <- estSURVTest
  Times <- time.interest

  # If newtimes specified, interpolate to those times
  if (!is.null(newtimes)) {
    Probs <- survprobMatInterpolator(probsMat = Probs, times = Times, newtimes = newtimes)
    Times <- newtimes
  }

  return(list(
    Probs = Probs,
    Times = Times
  ))
}
