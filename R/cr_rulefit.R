#' @title CRModel_rulefit
#'
#' @description Fit a rulefit model for competing risks outcomes using RuleFit algorithm
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0,1,2 where 0=censored, 1=event of interest, 2=competing event)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param ntree number of trees to fit to extract rules (default: 300)
#' @param nsample number of samples for each tree (default: 300)
#' @param keepvars these variables will be used in each bagging iteration
#' @param cuttimes these cut times are used for calculating pseudo observation
#' @param alpha numeric, the elastic net mixing parameter for glmnet (default: 0.5)
#' @param maxit integer, maximum number of iterations for glmnet (default: 2000)
#' @param ... additional parameters passed to glmnet::cv.glmnet
#'
#' @return a list with the following components:
#'   \item{rulefit_model}{the fitted RuleFit model object}
#'   \item{times}{vector of unique event times in the training data}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cr_rulefit"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{ctreelist}{list of fitted tree models}
#'   \item{ruleslist}{list of extracted rules from trees}
#'
#' @importFrom rpart rpart rpart.control
#' @importFrom partykit as.party
#' @importFrom pseudo pseudoci
#' @importFrom stats as.formula model.matrix quantile rpois runif sd
#' @importFrom survival Surv
#' @export
CRModel_rulefit <- function(data, expvars, timevar, eventvar, failcode = 1,
                           ntree = 300, nsample = 300, keepvars = NULL,
                           cuttimes = NULL, alpha = 0.5, maxit = 2000, ...) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("'expvars' must be a non-empty character vector")
  }
  if (!timevar %in% colnames(data)) {
    stop("'timevar' not found in data: ", timevar)
  }
  if (!eventvar %in% colnames(data)) {
    stop("'eventvar' not found in data: ", eventvar)
  }
  if (!is.numeric(failcode) || length(failcode) != 1 || !(failcode %in% c(1, 2))) {
    stop("'failcode' must be 1 or 2")
  }
  if (!all(expvars %in% colnames(data))) {
    missing_vars <- expvars[!expvars %in% colnames(data)]
    stop("expvars not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # ============================================================================
  # Variable Profiling
  # ============================================================================
  varprof <- VariableProfile(data, expvars)

  # ============================================================================
  # Model Fitting Logic
  # ============================================================================
  if (is.null(cuttimes)){
    x <- data[, timevar]
    cuttimes <- stats::quantile(x, c(.1,.25,.50,.70,.8))
    cuttimesReg <- stats::quantile(x, c(.25,.50,.60, .70, .80))
  }

  formClass <- stats::as.formula(paste("ClassVar ~.", collapse = ""))
  formReg <- stats::as.formula(paste("RegVar ~.", collapse = ""))
  formSurv <- stats::as.formula(paste("survival::Surv(",timevar,",", eventvar,"==", failcode, ") ~.", collapse = ""))

  ctreelist <- lapply(1:ntree, function(repi){
    sampcols <- union(keepvars, sample(expvars, min(length(expvars), sample(c(1:10), 1))))
    selmodel <- sample(c(2,3), 1)
    selmodel <- 2 # Force regression/pseudo-observation based trees

    if (selmodel == 1) { # Classification Tree (Original code had this branch, keeping for reference but not used due to selmodel<-2)
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)
      datasampl <- data[samprows, colnames(data) %in% usevars]

      datasampl$ClassVar <- as.character(datasampl[, colnames(datasampl) %in% timevar] < sample(cuttimes, 1) &
                                         datasampl[, colnames(datasampl) %in% eventvar] == failcode)
      datasampl <- datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
      rpcontrol <- rpart::rpart.control(
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
      rpartmodel <- rpart::rpart(formClass, data = datasampl, control = rpcontrol, model = TRUE)
    } else if (selmodel == 2) { # Regression Tree (Pseudo-observations)
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      pout <- pseudo::pseudoci(time=datasampl[,timevar], event=datasampl[,eventvar], tmax=sample(cuttimesReg,1))
      datasampl$RegVar <- c(unlist(pout$pseudo[[paste0("cause", failcode)]]))
      datasampl$RegVar <- (datasampl$RegVar - mean(datasampl$RegVar)) / stats::sd(datasampl$RegVar)
      datasampl <- datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
      rpcontrol <- rpart::rpart.control(
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
      rpartmodel <- rpart::rpart(formReg, data = datasampl, control = rpcontrol, model = TRUE)
    } else if (selmodel == 3) { # Survival Tree (Original code had this branch, keeping for reference but not used due to selmodel<-2)
      usevars <- c(timevar, eventvar, sampcols)
      samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)
      datasampl <- data[samprows, colnames(data) %in% usevars]
      rpcontrol <- rpart::rpart.control(
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
      rpartmodel <- rpart::rpart(formSurv, data = datasampl, control = rpcontrol, model = TRUE)
    }
    return(rpartmodel)
  })

  # Extract rules from trees
  ruleslist <- lapply(ctreelist, function(x) {
    tryCatch(extract_rules_from_party(partykit::as.party(x)), error = function(e) list())
  })

  # Create training matrix with rules
  RulesTrain <- lapply(ctreelist, function(x)
    predict(partykit::as.party(x), newdata = data, type = "node"))
  RulesTrainMatList <- lapply(seq_along(RulesTrain), function(x) {
    Train <- RulesTrain[[x]]
    rules <- ruleslist[[x]]
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
      Train <- factor(as.character(Train), levels = names(rules))
      if (length(unique(names(rules))) > 1) {
        MM <- stats::model.matrix(~ -1 + Train)
        colnames(MM) <- paste(paste(x, seq_along(rules), sep = "_"), rules, sep = "***")
        MM
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  # Filter out NULL elements before binding
  RulesTrainMatList <- Filter(Negate(is.null), RulesTrainMatList)

  if (length(RulesTrainMatList) > 0) {
    TrainMatRules <- Reduce("cbind", RulesTrainMatList)
    TrainMat <- cbind(stats::model.matrix(~.-1, data=data[,expvars]), TrainMatRules)
  } else {
    TrainMat <- stats::model.matrix(~.-1, data=data[,expvars])
  }

  usecols <- colnames(TrainMat)
  nrowforsample <- min(c(nrow(data), 500))

  # Fit the underlying Fine-Gray model
  CRrulefitModel <- CRModel_FineGray(
    data = cbind(TrainMat, data[, colnames(data) %in% c(eventvar, timevar)]),
    expvars = usecols,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode
  )

  # Get unique event times from training data
  times <- sort(unique(data[data[[eventvar]] != 0, timevar]))

  # Get time range
  time_range <- range(data[data[[eventvar]] != 0, timevar])

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    rulefit_model = list(
      ctreelist = ctreelist,
      ruleslist = ruleslist,
      CRrulefitModel = CRrulefitModel,
      usecols = usecols,
      datasamp = data[1:nrowforsample, expvars]
    ),
    times = times,
    varprof = varprof,
    model_type = "cr_rulefit",
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    failcode = failcode,
    time_range = time_range
  )

  class(result) <- "ml4t2e_cr_rulefit"
  return(result)
}


#' @title Predict_CRModel_rulefit
#'
#' @description Get predictions from a CR RuleFit model for a test dataset.
#'
#' @param modelout the output from 'CRModel_rulefit' (a list containing model and metadata)
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the model's native time points.
#'   Can be any positive values - interpolation handles all time points.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated
#'     (always includes time 0)}
#'
#' @importFrom partykit as.party
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_rulefit <- function(modelout, newdata, newtimes = NULL, failcode = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_cr_rulefit")) {
    stop("'modelout' must be output from CRModel_rulefit")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame")
  }

  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop("The following variables missing in newdata: ",
         paste(missing_vars, collapse = ", "))
  }

  # ============================================================================
  # Prepare newdata
  # ============================================================================
  # Ensure factor levels match training data
  newdata_prepared <- newdata[, modelout$expvars, drop = FALSE]

  for (var in modelout$expvars) {
    if (var %in% names(modelout$varprof)) {
      varinfo <- modelout$varprof[[var]]
      # Handle factors
      if (is.table(varinfo)) {
        training_levels <- names(varinfo)
        if (is.factor(newdata_prepared[[var]])) {
          # Ensure factor has same levels as training
          new_levels <- levels(newdata_prepared[[var]])
          extra_levels <- setdiff(new_levels, training_levels)
          if (length(extra_levels) > 0) {
            warning("Factor '", var, "' has new levels in newdata: ",
                    paste(extra_levels, collapse = ", "),
                    ". These will be set to NA.")
          }
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
                                            levels = training_levels)
        } else if (is.character(newdata_prepared[[var]])) {
          # Convert character to factor with training levels
          newdata_prepared[[var]] <- factor(newdata_prepared[[var]],
                                            levels = training_levels)
        }
      }
    }
  }

  # ============================================================================
  # Generate Rules for Test Data
  # ============================================================================
  RulesTest <- lapply(modelout$rulefit_model$ctreelist, function(x)
    predict(partykit::as.party(x), newdata = newdata_prepared, type = "node"))
  RulesTestMatList <- lapply(seq_along(RulesTest), function(x) {
    Test <- RulesTest[[x]]
    rules <- modelout$rulefit_model$ruleslist[[x]]
    if (length(rules) > 0 && length(unique(names(rules))) > 0) {
      Test <- factor(as.character(Test), levels = names(rules))
      if (length(unique(names(rules))) > 1) {
        MM <- stats::model.matrix(~ -1 + Test)
        colnames(MM) <- paste(paste(x, seq_along(rules), sep = "_"), rules, sep = "***")
        MM
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  # Filter out NULL elements before binding
  RulesTestMatList <- Filter(Negate(is.null), RulesTestMatList)

  # Create model matrix for original features, ensuring factor levels match training sample
  MM_orig <- stats::model.matrix(~.-1, data=rbind(modelout$rulefit_model$datasamp, newdata_prepared))
  MM_orig <- MM_orig[-seq_len(nrow(modelout$rulefit_model$datasamp)), , drop = FALSE] # Use drop=FALSE to keep matrix structure

  if (length(RulesTestMatList) > 0) {
    TestMatRules <- Reduce("cbind", RulesTestMatList)
    TestMat <- as.data.frame(cbind(MM_orig, TestMatRules))
  } else {
    TestMat <- as.data.frame(MM_orig) # Only original features if no rules generated
  }

  # Ensure TestMat has the same columns as used in the internal Fine-Gray model training
  missing_cols <- setdiff(modelout$rulefit_model$CRrulefitModel$expvars, colnames(TestMat))
  for(col in missing_cols){
    TestMat[[col]] <- 0 # Add missing columns with 0
  }
  TestMat <- TestMat[, modelout$rulefit_model$CRrulefitModel$expvars, drop = FALSE] # Reorder and select columns

  # ============================================================================
  # Make Predictions using underlying Fine-Gray model
  # ============================================================================
  pred_fg <- Predict_CRModel_FineGray(modelout$rulefit_model$CRrulefitModel, TestMat, newtimes = newtimes)

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    CIFs = pred_fg$CIFs,
    Times = pred_fg$Times
  )

  return(result)
}
