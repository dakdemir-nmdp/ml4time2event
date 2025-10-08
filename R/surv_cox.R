#' @title SurvModel_Cox
#'
#' @description Fit a Cox proportional hazards model for survival outcomes with optional variable selection.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param varsel character string specifying variable selection method:
#'   "none" (default, no selection),
#'   "backward" (backward elimination),
#'   "forward" (forward selection),
#'   "both" (stepwise selection),
#'   "penalized" (elastic net penalized Cox via glmnet)
#' @param penalty character string specifying penalty criterion for stepwise methods: "AIC" (default) or "BIC"
#' @param alpha numeric value in [0,1] for elastic net mixing parameter when varsel="penalized".
#'   alpha=1 is lasso, alpha=0 is ridge, alpha=0.5 (default) is elastic net.
#' @param nfolds integer, number of cross-validation folds for penalized Cox (default: 10)
#' @param ntimes integer, number of time points to use for prediction grid (default: 50)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{cph_model}{the fitted Cox model object (coxph or cv.glmnet)}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cox_standard" or "cox_penalized"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{varsel_method}{character string indicating variable selection method used}
#'
#' @note Predictions can be made at any time points via interpolation.
#'   No fixed time grid is stored in the model.
#'
#' @importFrom survival coxph Surv survfit
#' @importFrom stats as.formula model.matrix coef
#' @importFrom glmnet cv.glmnet
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate survival data
#' set.seed(123)
#' n <- 200
#' data <- data.frame(
#'   time = rexp(n, 0.1),
#'   event = rbinom(n, 1, 0.7),
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = factor(sample(c("A","B","C"), n, replace=TRUE))
#' )
#'
#' # Fit Cox model without variable selection
#' model1 <- SurvModel_Cox(data, expvars=c("x1","x2","x3"),
#'                         timevar="time", eventvar="event")
#'
#' # Fit Cox model with backward selection (AIC)
#' model2 <- SurvModel_Cox(data, expvars=c("x1","x2","x3"),
#'                         timevar="time", eventvar="event",
#'                         varsel="backward", penalty="AIC")
#'
#' # Fit penalized Cox model (elastic net)
#' model3 <- SurvModel_Cox(data, expvars=c("x1","x2","x3"),
#'                         timevar="time", eventvar="event",
#'                         varsel="penalized", alpha=0.5)
#' }
SurvModel_Cox <- function(data, expvars, timevar, eventvar,
                          varsel = "none",
                          penalty = "AIC",
                          alpha = 0.5,
                          nfolds = 10,
                          ntimes = 50,
                          verbose = FALSE) {

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
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following expvars not found in data: ", paste(missing_vars, collapse=", "))
  }

  varsel <- match.arg(varsel, c("none", "backward", "forward", "both", "penalized"))
  penalty <- match.arg(penalty, c("AIC", "BIC"))

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value between 0 and 1")
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1 || nfolds < 2) {
    stop("'nfolds' must be an integer >= 2")
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================
  if (verbose) cat("Creating variable profile...\n")
  varprof <- VariableProfile(data, expvars)

  # Ensure event variable is numeric 0/1
  data[[eventvar]] <- as.numeric(data[[eventvar]] == 1)

  # Prepare data subset with complete cases
  XYTrain <- data[, c(timevar, eventvar, expvars), drop = FALSE]
  n_original <- nrow(XYTrain)
  XYTrain <- XYTrain[complete.cases(XYTrain), , drop = FALSE]
  n_complete <- nrow(XYTrain)

  if (n_complete < n_original) {
    warning(sprintf("Removed %d rows with missing values. %d rows remaining.",
                    n_original - n_complete, n_complete))
  }
  if (n_complete < 10) {
    stop("Insufficient data after removing missing values. Need at least 10 observations.")
  }

  # Get unique event times (for reference, not used for fixed grid)
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == 1]
  if (length(event_times) == 0) {
    stop("No events in training data. Cannot fit survival model.")
  }

  # Store event time range for reference
  # Note: Predictions can be made at ANY time points via interpolation
  time_range <- c(0, max(event_times))

  # ============================================================================
  # Model Fitting
  # ============================================================================
  if (varsel == "penalized") {
    # ------------------------------------------------------------------------
    # Penalized Cox Model (glmnet)
    # ------------------------------------------------------------------------
    if (verbose) cat("Fitting penalized Cox model with elastic net (alpha=", alpha, ")...\n", sep="")

    # Create model matrix
    TrainMat <- stats::model.matrix(~ ., data = XYTrain[, expvars, drop = FALSE])
    # Remove intercept column if present
    intercept_col <- which(colnames(TrainMat) == "(Intercept)")
    if (length(intercept_col) > 0) {
      TrainMat <- TrainMat[, -intercept_col, drop = FALSE]
    }

    # Create survival object
    y_surv <- survival::Surv(XYTrain[[timevar]], XYTrain[[eventvar]])

    # Fit CV glmnet
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = TrainMat,
        y = y_surv,
        family = "cox",
        alpha = alpha,
        nfolds = nfolds,
        type.measure = "deviance"
      ),
      error = function(e) {
        stop("Failed to fit penalized Cox model: ", e$message)
      }
    )

    # Store sample of training data and matrix for prediction
    sample_size <- min(500, nrow(XYTrain))
    sample_idx <- sample(seq_len(nrow(XYTrain)), sample_size)

    cph_model <- list(
      cv_fit = cv_fit,
      train_sample = XYTrain[sample_idx, expvars, drop = FALSE],
      train_matrix = TrainMat[sample_idx, , drop = FALSE],
      y_train = y_surv[sample_idx],
      alpha = alpha
    )
    class(cph_model) <- c("ml4t2e_cox_penalized", "list")
    model_type <- "cox_penalized"

  } else {
    # ------------------------------------------------------------------------
    # Standard Cox Model (with optional stepwise selection)
    # ------------------------------------------------------------------------
    # Define formula
    form <- stats::as.formula(
      paste0("survival::Surv(", timevar, ", ", eventvar, ") ~ ",
             paste(expvars, collapse = " + "))
    )

    if (verbose) cat("Fitting Cox model...\n")

    # Fit initial Cox model
    cph_model <- tryCatch(
      survival::coxph(form, data = XYTrain, x = TRUE, y = TRUE),
      error = function(e) {
        stop("Failed to fit Cox model: ", e$message)
      }
    )

    # Variable selection
    if (varsel != "none") {
      if (verbose) cat("Performing variable selection (", varsel, ", ", penalty, ")...\n", sep="")

      # Set penalty parameter k for step()
      k_penalty <- switch(penalty, "AIC" = 2, "BIC" = log(nrow(XYTrain)))

      # Define scope for forward/both methods
      if (varsel %in% c("forward", "both")) {
        null_formula <- stats::as.formula(
          paste0("survival::Surv(", timevar, ", ", eventvar, ") ~ 1")
        )
        null_model <- survival::coxph(null_formula, data = XYTrain, x = TRUE, y = TRUE)
        lower_scope <- null_formula
        upper_scope <- form
      }

      # Perform stepwise selection
      cph_model <- tryCatch(
        {
          if (varsel == "backward") {
            stats::step(cph_model, direction = "backward", k = k_penalty, trace = as.numeric(verbose))
          } else if (varsel == "forward") {
            stats::step(null_model, direction = "forward", scope = list(lower = lower_scope, upper = upper_scope),
                       k = k_penalty, trace = as.numeric(verbose))
          } else if (varsel == "both") {
            stats::step(null_model, direction = "both", scope = list(lower = lower_scope, upper = upper_scope),
                       k = k_penalty, trace = as.numeric(verbose))
          }
        },
        error = function(e) {
          warning("Variable selection failed: ", e$message, ". Using full model.")
          cph_model # Return original full model
        }
      )
    }

    model_type <- "cox_standard"
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cox model fitting complete.\n")

  result <- list(
    cph_model = cph_model,
    time_range = time_range,  # Store time range instead of fixed grid
    varprof = varprof,
    model_type = model_type,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    varsel_method = varsel
  )

  class(result) <- c("ml4t2e_surv_cox", "list")
  return(result)
}

#' @title Predict_SurvModel_Cox
#'
#' @description Get predictions from a fitted Cox survival model for new data.
#'
#' @param modelout the output from 'SurvModel_Cox'
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), generates 50 equally-spaced points from 0 to max observed time.
#'   Can be any positive values - interpolation handles all time points.
#'
#' @return a list containing:
#'   \item{Probs}{predicted survival probability matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which probabilities are calculated
#'     (always includes time 0)}
#'   \item{survfit_obj}{the raw survfit object from prediction
#'     (for diagnostics)}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit model
#' model <- SurvModel_Cox(train_data, expvars = c("x1", "x2"),
#'                        timevar = "time", eventvar = "event")
#'
#' # Predict on test data
#' preds <- Predict_SurvModel_Cox(model, test_data)
#'
#' # Predict at specific times
#' preds_custom <- Predict_SurvModel_Cox(model, test_data,
#'                                       newtimes = c(30, 60, 90, 180, 365))
#' }
Predict_SurvModel_Cox <- function(modelout, newdata, newtimes = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!inherits(modelout, "ml4t2e_surv_cox")) {
    stop("'modelout' must be output from SurvModel_Cox")
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

  # Generate default times if not specified
  if (is.null(newtimes)) {
    # Create a reasonable default grid: 50 points from 0 to max observed time
    max_time <- modelout$time_range[2]
    newtimes <- seq(0, max_time, length.out = 50)
  } else {
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))
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

  # Check for NAs
  if (any(is.na(newdata_prepared))) {
    warning("Missing values in newdata will result in NA predictions")
  }

  # ============================================================================
  # Make Predictions
  # ============================================================================
  if (modelout$model_type == "cox_penalized") {
    # ------------------------------------------------------------------------
    # Penalized Cox Predictions
    # ------------------------------------------------------------------------
    # Create model matrix ensuring factor levels match
    combined_data <- rbind(
      modelout$cph_model$train_sample,
      newdata_prepared
    )
    combined_mat <- stats::model.matrix(~ ., data = combined_data)
    intercept_col <- which(colnames(combined_mat) == "(Intercept)")
    if (length(intercept_col) > 0) {
      combined_mat <- combined_mat[, -intercept_col, drop = FALSE]
    }

    # Extract test portion
    n_train <- nrow(modelout$cph_model$train_sample)
    test_mat <- combined_mat[-seq_len(n_train), , drop = FALSE]

    # Ensure columns match training
    train_cols <- colnames(modelout$cph_model$train_matrix)
    missing_cols <- setdiff(train_cols, colnames(test_mat))
    for (col in missing_cols) {
      test_mat <- cbind(test_mat, 0)
      colnames(test_mat)[ncol(test_mat)] <- col
    }
    test_mat <- test_mat[, train_cols, drop = FALSE]

    # Get predictions using survfit
    survfit_obj <- tryCatch(
      survival::survfit(
        modelout$cph_model$cv_fit,
        s = modelout$cph_model$cv_fit$lambda.min,
        x = modelout$cph_model$train_matrix,
        y = modelout$cph_model$y_train,
        newx = test_mat
      ),
      error = function(e) {
        stop("Prediction failed: ", e$message)
      }
    )

  } else {
    # ------------------------------------------------------------------------
    # Standard Cox Predictions
    # ------------------------------------------------------------------------
    survfit_obj <- tryCatch(
      survival::survfit(
        modelout$cph_model,
        newdata = newdata_prepared
      ),
      error = function(e) {
        stop("Prediction failed: ", e$message)
      }
    )
  }

  # ============================================================================
  # Extract and Format Predictions
  # ============================================================================
  # Get survival probabilities and times
  # survfit returns: surv matrix is times x observations
  pred_probs <- survfit_obj$surv
  pred_times <- survfit_obj$time

  # Ensure matrix format (handle single observation case)
  if (!is.matrix(pred_probs)) {
    pred_probs <- matrix(pred_probs, ncol = 1)
  }

  # Use the standard interpolation utility function
  # This handles: time 0 addition, interpolation to newtimes, and monotonicity
  pred_probs <- survprobMatInterpolator(
    probsMat = pred_probs,
    times = pred_times,
    newtimes = newtimes
  )
  pred_times <- newtimes

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    Probs = pred_probs,
    Times = pred_times,
    survfit_obj = survfit_obj
  )

  return(result)
}
