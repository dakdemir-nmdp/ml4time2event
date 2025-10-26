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
#'   \item{times}{vector of time points for prediction grid (equally-spaced from 0 to max time)}
#'   \item{time_range}{vector with min and max observed event times}
#'   \item{varprof}{variable profile list containing factor levels and numeric ranges}
#'   \item{model_type}{character string "cox_standard" or "cox_penalized"}
#'   \item{expvars}{character vector of explanatory variables used}
#'   \item{timevar}{character name of time variable}
#'   \item{eventvar}{character name of event variable}
#'   \item{varsel_method}{character string indicating variable selection method used}
#'   \item{alpha}{numeric value of elastic net mixing parameter used (if varsel="penalized")}
#'   \item{nfolds}{integer number of cross-validation folds used (if varsel="penalized")}
#'
#' @note Predictions can be made at any time points via interpolation.
#'   The times vector is used as default when predicting, but can be overridden.
#'
#' @importFrom survival coxph Surv survfit
#' @importFrom stats as.formula model.matrix coef complete.cases
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
  varprof <- varprof[expvars]

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
  
  # Generate default times grid for predictions
  times <- seq(time_range[1], time_range[2], length.out = ntimes)

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

    # Check for at least one nonzero coefficient
    coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
    if (all(coefs == 0 | is.na(coefs))) {
      stop("All coefficients are zero after penalized Cox fitting. The final model has no predictors. This is not a valid Cox model. Try using a different variable selection method or check your data.")
    }

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
      # Always keep the single best variable if none selected (robust variable selection)
      coefs <- coef(cph_model)
      if (length(coefs) == 0 || all(coefs == 0 | is.na(coefs))) {
        # Try all single-variable models and pick the best by AIC/BIC
        best_var <- NULL
        best_aic <- Inf
        for (var in expvars) {
          single_form <- stats::as.formula(paste0("survival::Surv(", timevar, ", ", eventvar, ") ~ ", var))
          single_model <- tryCatch(
            survival::coxph(single_form, data = XYTrain, x = TRUE, y = TRUE),
            error = function(e) NULL
          )
          if (!is.null(single_model)) {
            aic_val <- switch(penalty, "AIC" = AIC(single_model), "BIC" = AIC(single_model, k = log(nrow(XYTrain))))
            if (aic_val < best_aic) {
              best_aic <- aic_val
              best_var <- var
            }
          }
        }
        if (!is.null(best_var)) {
          if (verbose) cat("No variable improved over intercept-only, keeping best single variable:", best_var, "AIC/BIC:", round(best_aic, 2), "\n")
          cph_model <- survival::coxph(
            stats::as.formula(paste0("survival::Surv(", timevar, ", ", eventvar, ") ~ ", best_var)),
            data = XYTrain, x = TRUE, y = TRUE
          )
        } else {
          stop("No valid single-variable Cox model could be fit after variable selection.")
        }
      }
    }

    model_type <- "cox_standard"
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("Cox model fitting complete.\n")

  result <- list(
    cph_model = cph_model,
    times = times,
    time_range = time_range,
    varprof = varprof,
    model_type = model_type,
    expvars = expvars,
    timevar = timevar,
    eventvar = eventvar,
    varsel_method = varsel,
    alpha = alpha,
    nfolds = nfolds
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
#' @param new_times optional numeric vector of time points for prediction.
#'   If NULL (default), uses the times generated during model training
#'   (50 equally-spaced points from 0 to max observed time).
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
#'                                       new_times = c(30, 60, 90, 180, 365))
#' }
Predict_SurvModel_Cox <- function(modelout, newdata, new_times = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("'modelout' is missing")
  }
  if (!is.list(modelout) || !all(c("expvars", "time_range") %in% names(modelout))) {
    stop("must be output from SurvModel_Cox")
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
  if (is.null(new_times)) {
    # Create a reasonable default grid: 50 points from 0 to max observed time
    max_time <- modelout$time_range[2]
    new_times <- seq(0, max_time, length.out = 50)
  } else {
    if (!is.numeric(new_times) || any(new_times < 0)) {
      stop("'new_times' must be a numeric vector of non-negative values")
    }
    new_times <- sort(unique(new_times))
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

  # Check for NAs and count complete cases
  has_missing <- any(is.na(newdata_prepared))
  n_complete_cases <- sum(complete.cases(newdata_prepared))

  if (has_missing) {
    warning(sprintf(
      "Missing values detected in %d observation(s). ",
      nrow(newdata_prepared) - n_complete_cases),
      "survfit will drop incomplete cases."
    )
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
      glmnet:::survfit.coxnet(
        modelout$cph_model$cv_fit,
        s = "lambda.min",
        newx = test_mat,
        x = modelout$cph_model$train_matrix,
        y = modelout$cph_model$y_train
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

  # Determine expected number of observations
  n_expected_obs <- nrow(newdata_prepared)

  # Ensure matrix format
  if (!is.matrix(pred_probs)) {
    # survfit returned a vector - this happens when:
    # 1. Single observation (expected)
    # 2. Multiple observations but survfit collapsed to single curve (unexpected)

    if (n_expected_obs == 1) {
      # Single observation case - shape as column vector
      pred_probs <- matrix(pred_probs, ncol = 1)
    } else {
      # Check for all-identical covariate patterns
      all_identical <- all(apply(newdata_prepared, 2, function(col) length(unique(col)) == 1))
      # Check for excessive missingness
      all_missing_but_one <- (n_complete_cases == 1)
      if (all_identical) {
        # No warning: all covariate patterns are identical, expected behavior
      } else if (all_missing_but_one) {
        # No warning: all but one row had missing values, expected behavior
      } else {
        warning(sprintf(
          paste0(
            "survfit returned a vector for %d observations. This is unexpected. ",
            "Possible causes: identical risk scores, excessive missingness, or a bug. ",
            "First few rows of newdata:\n%s\n",
            "Treating as single survival curve replicated across all observations."),
          n_expected_obs,
          paste(capture.output(print(utils::head(newdata_prepared, 3))), collapse = "\n")
        ))
      }
      # Create matrix: replicate the single curve for all observations
      pred_probs <- matrix(rep(pred_probs, n_expected_obs),
                          nrow = length(pred_probs),
                          ncol = n_expected_obs)
    }
  }

  # Validate dimensions
  n_returned_obs <- ncol(pred_probs)

  if (n_returned_obs != n_expected_obs) {
    # Check if the mismatch is due to missing values (expected)
    if (has_missing && n_returned_obs == n_complete_cases) {
      # This is expected - survfit dropped incomplete cases
      # Expand predictions to match original newdata by filling with NAs
      full_pred_probs <- matrix(NA, nrow = nrow(pred_probs), ncol = n_expected_obs)
      complete_idx <- which(complete.cases(newdata_prepared))
      full_pred_probs[, complete_idx] <- pred_probs
      pred_probs <- full_pred_probs
    } else {
      # Unexpected mismatch
      stop(sprintf(
        "Dimension mismatch: survfit returned %d observation(s) but expected %d.\n",
        n_returned_obs, n_expected_obs),
        sprintf("Complete cases: %d. ", n_complete_cases),
        "This indicates an unexpected issue in survfit call or data preparation."
      )
    }
  }

  # Use the standard interpolation utility function
  # This handles: time 0 addition, interpolation to new_times, and monotonicity
  pred_probs <- survprobMatInterpolator(
    probsMat = pred_probs,
    times = pred_times,
    new_times = new_times
  )
  pred_times <- new_times

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
