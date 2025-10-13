#' @title (Helper) Extract Rules from a party Object
#'
#' @description This helper function traverses a party tree
#' and extracts the decision rules for each terminal node.
#' @param x A party object.
#' @return A named character vector where names are the node IDs and values are the rules.
#' @keywords internal
extract_rules_from_party <- function(x) {
  if (is.null(x) || is.null(x$node)) return(list())
  # Get rule definitions from partykit's internal functions
  rule_strings <- partykit:::.list.rules.party(x)
  # Get the ID of the terminal node for each rule
  node_ids <- names(partykit::nodeids(x, terminal = TRUE))
  # Name the rules by their terminal node ID
  names(rule_strings) <- node_ids
  return(rule_strings)
}

#' @title Fit a RuleFit Model for Survival Outcomes
#'
#' @description Implements the RuleFit algorithm by generating rules from an ensemble of
#' decision trees and fitting a penalized Cox model (via glmnet) on the original
#' features plus the derived rules.
#'
#' @param data A data frame with explanatory and outcome variables.
#' @param expvars Character vector of names of explanatory variables in `data`.
#' @param timevar Character name of the time-to-event variable in `data`.
#' @param eventvar Character name of the event indicator variable in `data`.
#' @param ntree Integer, the number of trees to generate for rule extraction.
#' @param nsample Integer, the number of samples to draw (with replacement) for training each tree.
#' @param keepvars Character vector of variable names to be included in every tree's bagging iteration.
#' @param cuttimes Numeric vector of time points used to create binary outcomes for classification trees.
#'   If NULL, quantiles of the time variable are used.
#' @param alpha Numeric, the elastic net mixing parameter (default: 0.5).
#' @param maxit Integer, maximum number of iterations for glmnet (default: 2000).
#' @param ... Additional parameters passed to glmnet::cv.glmnet.
#'
#' @return A list with the following components:
#' \item{model}{The fitted RuleFit model object containing all components.}
#' \item{times}{Unique event times observed in the training data.}
#' \item{varprof}{A profile of the explanatory variables.}
#' \item{expvars}{The names of the original explanatory variables.}
#' \item{factor_levels}{A list containing factor levels for each variable to ensure consistent prediction.}
#'
#' @importFrom rpart rpart rpart.control
#' @importFrom partykit as.party
#' @importFrom glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix quantile rpois runif coef
#' @importFrom survival Surv survfit
#'
#' @export
#' @family RuleFit Survival Modeling
SurvModel_rulefit <- function(data,
                              expvars,
                              timevar,
                              eventvar,
                              ntree = 300,
                              nsample = 300,
                              keepvars = NULL,
                              cuttimes = NULL,
                              alpha = 0.5,
                              maxit = 2000,
                              ...) {

  if (missing(data)) stop("argument \"data\" is missing")
  if (missing(expvars)) stop("argument \"expvars\" is missing")
  if (missing(timevar)) stop("argument \"timevar\" is missing")
  if (missing(eventvar)) stop("argument \"eventvar\" is missing")

  # --- 1. Initial Setup ---
  varprof <- VariableProfile(data, expvars)

  # Store factor levels for prediction
  factor_levels <- lapply(data[, expvars, drop=FALSE], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })

  if (is.null(cuttimes)) {
    x <- data[[timevar]]
    cuttimes <- stats::quantile(x, c(.1, .25, .50, .70, .90), na.rm = TRUE)
  }

  # Define model formulas
  formRF1 <- stats::as.formula(paste("survival::Surv(", timevar, ",", eventvar, ") ~ ."))
  formRF2 <- stats::as.formula("ClassVar ~ .")


  # --- 2. Generate Ensemble of Trees ---
  # This loop generates a diverse set of shallow trees by randomizing hyperparameters.
  # Both survival and classification trees are used to capture different patterns.
  ctreelist <- lapply(1:ntree, function(i) {
    # Bagging: sample columns and rows
    sampcols <- union(keepvars, sample(expvars, min(length(expvars), sample(1:10, 1))))
    samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)

    # Randomly choose between a survival tree (favored) and a classification tree
    selmodel <- sample(c(1, 1, 1, 1, 2), 1)
    formRF <- if (selmodel == 1) formRF1 else formRF2

    usevars <- c(timevar, eventvar, sampcols)
    datasampl <- data[samprows, usevars, drop = FALSE]

    # For classification trees, create a binary outcome based on a random cut time
    if (selmodel == 2) {
      datasampl$ClassVar <- as.factor(
        datasampl[[timevar]] < sample(cuttimes, 1) & datasampl[[eventvar]] == 1
      )
      datasampl <- datasampl[, !colnames(datasampl) %in% c(timevar, eventvar)]
    }

    # Randomized hyperparameters to foster tree diversity
    rpcontrol <- rpart::rpart.control(
      minsplit = stats::rpois(1, 2) + 1,
      minbucket = stats::rpois(1, 20) + 5,
      cp = 0.01 * stats::runif(1),
      maxdepth = stats::rpois(1, 2) + 1,
      usesurrogate = 0, # Turn off surrogates for simplicity and speed
      xval = 0 # No cross-validation needed for individual trees
    )

    rpart::rpart(formRF, data = datasampl, control = rpcontrol, model = TRUE)
  })


  # --- 3. Extract Rules and Create Rule Matrix ---
  ruleslist <- lapply(ctreelist, function(x) {
    tryCatch(extract_rules_from_party(partykit::as.party(x)), error = function(e) list())
  })

  RulesTrainMatList <- lapply(seq_along(ctreelist), function(i) {
    tree <- ctreelist[[i]]
    rules <- ruleslist[[i]]
    if (length(rules) == 0 || length(unique(names(rules))) <= 1) {
      return(NULL)
    }

    # Get terminal node predictions and one-hot encode them to form the rule matrix
    node_preds <- predict(partykit::as.party(tree), newdata = data, type = "node")
    train_factor <- factor(as.character(node_preds), levels = names(rules))

    MM <- stats::model.matrix(~ -1 + train_factor)
    colnames(MM) <- paste(paste(i, seq_along(rules), sep = "_"), rules, sep = "***")
    MM
  })
  RulesTrainMatList <- Filter(Negate(is.null), RulesTrainMatList)


  # --- 4. Combine Original Features and Rules ---
  TrainMatOrig <- stats::model.matrix(~ . - 1, data = data[, expvars, drop = FALSE])

  if (length(RulesTrainMatList) > 0) {
    TrainMatRules <- Reduce("cbind", RulesTrainMatList)
    TrainMat <- cbind(TrainMatOrig, TrainMatRules)
  } else {
    TrainMat <- TrainMatOrig
  }


  # --- 5. Fit Penalized Cox Model ---
  yTrain <- survival::Surv(data[[timevar]], data[[eventvar]] == 1)

  cv.fitRules <- glmnet::cv.glmnet(
    x = TrainMat,
    y = yTrain,
    alpha = alpha, # User-provided Elastic Net penalty
    family = "cox",
    maxit = maxit, # User-provided max iterations
    ...
  )

  # Get unique event times from the training data for prediction horizon
  sfitTrain <- survival::survfit(yTrain ~ 1)
  timesTrain <- sfitTrain$time


  # --- 6. Prepare and Return Model Object ---
  # Sample of original data is kept to ensure factor levels are consistent during prediction
  nrowforsample <- min(c(nrow(data), 500))
  datasamp <- data[sample(nrow(data), nrowforsample), expvars, drop = FALSE]

  model_output <- list(
    datasamp = datasamp,
    varprof = varprof,
    ruleslist = ruleslist,
    ctreelist = ctreelist,
    TrainMat = TrainMat, # Return the full matrix for correct prediction
    yTrain = yTrain,     # Return the full outcome vector for correct prediction
    cv.fitRules = cv.fitRules,
    timesTrain = timesTrain,
    expvars = expvars
  )

  # Return standardized output
  result <- list(model = model_output, times = timesTrain, varprof = varprof, expvars = expvars, factor_levels = factor_levels)
  class(result) <- c("ml4t2e_surv_rulefit", "list")
  return(result)
}


#' @title Predict Survival Probabilities from a RuleFit Model
#'
#' @description Generates survival probability predictions for new data using a
#' previously trained RuleFit model.
#'
#' @param modelout The output from `SurvModel_rulefit` (a list containing 'model', 'times', 'varprof', 'expvars', 'factor_levels').
#' @param newdata A new data frame for which to generate predictions. It must contain
#'   all variables specified in `modelout$expvars`.
#'
#' @return A list with two items:
#' \item{Probs}{A matrix of survival probability predictions, where rows are time points and columns are observations.}
#' \item{Times}{The vector of time points corresponding to the rows of `Probs`.}
#'
#' @importFrom partykit as.party
#' @importFrom stats model.matrix
#' @importFrom survival survfit
#'
#' @export
Predict_SurvModel_rulefit <- function(modelout, newdata, newtimes = NULL) {

  if (missing(modelout)) stop("argument \"modelout\" is missing")
  if (missing(newdata)) stop("argument \"newdata\" is missing")

  # Prepare newdata: ensure factors have same levels as training data
  data_test <- newdata[, modelout$expvars, drop=FALSE]
  for (vari in modelout$expvars){
    if (!is.null(modelout$factor_levels[[vari]])) { # Check if var was factor in training
        train_levels <- modelout$factor_levels[[vari]]
        # Ensure the column exists in newdata before attempting to modify
        if (vari %in% colnames(data_test)) {
            # Convert to character first to handle potential new levels, then factor
            data_test[[vari]] <- factor(as.character(data_test[[vari]]), levels = train_levels)
        }
    }
  }

  # --- 1. Create Rule Matrix for New Data ---
  RulesTestMatList <- lapply(seq_along(modelout$model$ctreelist), function(i) {
    tree <- modelout$model$ctreelist[[i]]
    rules <- modelout$model$ruleslist[[i]]
    if (length(rules) == 0 || length(unique(names(rules))) <= 1) {
      return(NULL)
    }

    node_preds <- predict(partykit::as.party(tree), newdata = data_test, type = "node")
    test_factor <- factor(as.character(node_preds), levels = names(rules))

    MM <- stats::model.matrix(~ -1 + test_factor)
    colnames(MM) <- paste(paste(i, seq_along(rules), sep = "_"), rules, sep = "***")
    MM
  })
  RulesTestMatList <- Filter(Negate(is.null), RulesTestMatList)


  # --- 2. Create Original Feature Matrix for New Data ---
  # The rbind trick ensures that factor levels in newdata are handled identically
  # to the training data.
  full_data <- rbind(modelout$model$datasamp, data_test)
  MM_orig_full <- stats::model.matrix(~ . - 1, data = full_data)
  MM_orig_test <- MM_orig_full[-seq_len(nrow(modelout$model$datasamp)), , drop = FALSE]


  # --- 3. Combine and Align Feature Matrices ---
  if (length(RulesTestMatList) > 0) {
    TestMatRules <- Reduce("cbind", RulesTestMatList)
    TestMat <- cbind(MM_orig_test, TestMatRules)
  } else {
    TestMat <- MM_orig_test
  }

  # Ensure the new data matrix has exactly the same columns in the same order
  # as the training matrix. Add missing columns (rules not triggered by new data)
  # and remove extra columns (should not happen but is a safeguard).
  train_cols <- colnames(modelout$model$TrainMat)
  test_cols <- colnames(TestMat)

  missing_cols <- setdiff(train_cols, test_cols)
  for (col in missing_cols) {
    TestMat <- cbind(TestMat, 0)
    colnames(TestMat)[ncol(TestMat)] <- col
  }
  TestMat <- TestMat[, train_cols, drop = FALSE]


  # --- 4. Predict Survival Curves ---
  # It is CRITICAL to provide the original training data (x and y) to survfit.
  # This is because it needs them to compute the baseline hazard function, which
  # is then adjusted for each new observation based on its predictor values.
  sfitTestRules <- survival::survfit(
    modelout$model$cv.fitRules$glmnet.fit,
    s = modelout$model$cv.fitRules$lambda.min,
    x = modelout$model$TrainMat, # Full training matrix
    y = modelout$model$yTrain,   # Full training outcome
    newx = TestMat        # New data matrix for prediction
  )

  predSurvsTestRules <- sfitTestRules$surv
  timesTest <- sfitTestRules$time

  # Ensure time 0 with probability 1 is included for complete curves
  if (sum(timesTest == 0) == 0) {
    timesTest <- c(0, timesTest)
    predSurvsTestRules <- rbind(rep(1, ncol(predSurvsTestRules)), predSurvsTestRules)
  }

  Probs <- predSurvsTestRules
  Times <- timesTest

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
