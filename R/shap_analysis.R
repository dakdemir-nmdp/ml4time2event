#' @title Create SHAP Prediction Function
#'
#' @description Creates a prediction wrapper function that takes raw data,
#' applies the pipeline's preprocessing, and returns scalar predictions
#' (expected time lost) suitable for SHAP analysis.
#'
#' @param pipeline A fitted ml4t2e_pipeline object from [ml4t2e_fit_pipeline()].
#' @param time_horizon Numeric value specifying the upper limit for expected
#'   time lost calculation. Should be a clinically meaningful time point
#'   (e.g., 365 for 1-year, 730 for 2-year).
#' @param ensemble_method Character string specifying ensemble method for
#'   predictions. Default is "average". Options: "average", "weighted", "stacking".
#'
#' @return A function that takes a data frame of raw (unprocessed) data and
#'   returns a numeric vector of expected time lost predictions.
#'
#' @details
#' The returned prediction function:
#' \itemize{
#'   \item Accepts raw data (before preprocessing)
#'   \item Applies the pipeline's preprocessing recipe
#'   \item Generates ensemble predictions (survival curves or CIFs)
#'   \item Calculates expected time lost up to time_horizon
#'   \item Returns scalar values suitable for SHAP analysis
#' }
#'
#' Expected time lost is calculated as the integral of the event probability
#' from 0 to time_horizon. For survival models, this is the integral of
#' \eqn{1 - S(t)}. For competing risks, it's the integral of the CIF.
#'
#' @examples
#' \dontrun{
#' library(ml4time2event)
#'
#' # Fit a pipeline
#' lung_df <- get_lung_survival_data()
#' pipeline <- ml4t2e_fit_pipeline(
#'   data = lung_df,
#'   analysis_type = "survival",
#'   timevar = "time",
#'   eventvar = "status",
#'   models = c("coxph", "glmnet")
#' )
#'
#' # Create prediction function for 1-year expected time lost
#' pred_fn <- ml4t2e_shap_predict_fn(pipeline, time_horizon = 365)
#'
#' # Use with new data
#' predictions <- pred_fn(lung_df[1:10, ])
#' }
#'
#' @export
ml4t2e_shap_predict_fn <- function(pipeline,
                                    time_horizon,
                                    ensemble_method = "average") {

  # Validate inputs
  if (!inherits(pipeline, "ml4t2e_pipeline")) {
    stop("'pipeline' must be an ml4t2e_pipeline object created with ml4t2e_fit_pipeline().",
         call. = FALSE)
  }

  if (!is.numeric(time_horizon) || length(time_horizon) != 1 ||
      is.na(time_horizon) || time_horizon <= 0) {
    stop("'time_horizon' must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.character(ensemble_method) || length(ensemble_method) != 1) {
    stop("'ensemble_method' must be a character string.", call. = FALSE)
  }

  # Extract pipeline components
  analysis_type <- pipeline$analysis_type
  recipe <- pipeline$recipe
  model <- pipeline$model
  timevar <- pipeline$timevar
  eventvar <- pipeline$eventvar

  # Create prediction grid that includes the time_horizon
  # Use a fine grid for accurate integration
  max_time <- max(time_horizon, max(pipeline$prediction_grid, na.rm = TRUE))
  prediction_times <- sort(unique(c(0, seq(0, max_time, length.out = 100), time_horizon)))

  # Return the prediction function
  function(newdata) {

    if (!is.data.frame(newdata)) {
      stop("Input to prediction function must be a data.frame.", call. = FALSE)
    }

    # Prepare newdata to have all required columns
    required_cols <- unique(c(timevar, eventvar, pipeline$original_expvars, pipeline$idvars))
    missing_cols <- setdiff(required_cols, colnames(newdata))

    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        newdata[[col]] <- NA
      }
    }

    # Apply preprocessing using the pipeline's recipe
    tryCatch({
      baked_data <- bake_data_recipe(recipe, data = newdata)
    }, error = function(e) {
      stop("Error applying preprocessing recipe: ", e$message, call. = FALSE)
    })

    # Generate predictions using the pipeline's model
    if (identical(analysis_type, "survival")) {
      preds <- tryCatch({
        suppressMessages(
          PredictSurvModels(
            models = model,
            newdata = baked_data,
            new_times = prediction_times,
            ensemble_method = ensemble_method
          )
        )
      }, error = function(e) {
        stop("Error generating survival predictions: ", e$message, call. = FALSE)
      })

      # preds should contain NewProbs (survival probabilities)
      if (is.null(preds$NewProbs)) {
        stop("Prediction output missing 'NewProbs' component.", call. = FALSE)
      }

      survival_probs <- preds$NewProbs  # Matrix: rows = times, cols = observations
      pred_times <- preds$NewTimes

      # Calculate expected time lost for each observation
      # ETL = integral of (1 - S(t)) from 0 to time_horizon
      event_probs <- 1 - survival_probs

    } else {
      # Competing risks
      preds <- tryCatch({
        suppressMessages(
          PredictCRModels(
            models = model,
            newdata = baked_data,
            new_times = prediction_times,
            ensemble_method = ensemble_method
          )
        )
      }, error = function(e) {
        stop("Error generating competing risks predictions: ", e$message, call. = FALSE)
      })

      if (is.null(preds$NewProbs)) {
        stop("Prediction output missing 'NewProbs' component.", call. = FALSE)
      }

      # For CR, NewProbs contains CIF (already event probability)
      event_probs <- preds$NewProbs
      pred_times <- preds$NewTimes
    }

    # Calculate expected time lost using integration
    # Use the Integrator function with time_horizon as upper limit
    time_lost <- apply(event_probs, 2, function(obs_probs) {
      Integrator(
        times = pred_times,
        scores = obs_probs,
        minmax = c(0, time_horizon),
        scale = FALSE
      )
    })

    # Return scalar predictions
    as.numeric(time_lost)
  }
}


#' @title Calculate SHAP Values for ml4time2event Pipeline
#'
#' @description Calculates SHAP (SHapley Additive exPlanations) values for
#' a fitted ml4time2event pipeline. SHAP values explain individual predictions
#' by showing the contribution of each feature.
#'
#' @param pipeline A fitted ml4t2e_pipeline object from [ml4t2e_fit_pipeline()].
#' @param data Data frame containing observations to explain. Should contain
#'   the same variables used to train the pipeline (raw data, before preprocessing).
#' @param time_horizon Numeric value specifying the time horizon for expected
#'   time lost calculation.
#' @param background Optional data frame to use as background/reference for
#'   SHAP calculation. If NULL (default), uses a random sample from `data`.
#' @param nsim Number of Monte Carlo simulations for SHAP calculation.
#'   Higher values give more accurate estimates but take longer. Default is 100.
#' @param ... Additional arguments passed to the SHAP calculation function.
#'
#' @return A list of class "ml4t2e_shap" containing:
#' \describe{
#'   \item{shap_values}{Matrix of SHAP values (rows = observations, cols = features)}
#'   \item{baseline}{Baseline prediction (expected value over background data)}
#'   \item{predictions}{Actual predictions for each observation}
#'   \item{feature_values}{Original feature values for each observation}
#'   \item{feature_names}{Names of features}
#'   \item{time_horizon}{Time horizon used for expected time lost}
#'   \item{analysis_type}{Type of analysis ("survival" or "competing_risks")}
#' }
#'
#' @details
#' This function uses Kernel SHAP (or fast SHAP if available) to calculate
#' feature contributions. SHAP values satisfy:
#' \deqn{prediction = baseline + sum(SHAP values)}
#'
#' The calculation works on raw (unprocessed) features, so SHAP values
#' directly reflect the contribution of original variables.
#'
#' @examples
#' \dontrun{
#' library(ml4time2event)
#'
#' # Fit pipeline
#' lung_df <- get_lung_survival_data()
#' pipeline <- ml4t2e_fit_pipeline(
#'   data = lung_df,
#'   analysis_type = "survival",
#'   timevar = "time",
#'   eventvar = "status",
#'   models = c("coxph", "glmnet")
#' )
#'
#' # Calculate SHAP values
#' shap_result <- ml4t2e_calculate_shap(
#'   pipeline = pipeline,
#'   data = lung_df[1:50, ],
#'   time_horizon = 365,
#'   nsim = 100
#' )
#'
#' # Examine results
#' print(shap_result)
#' }
#'
#' @export
ml4t2e_calculate_shap <- function(pipeline,
                                   data,
                                   time_horizon,
                                   background = NULL,
                                   nsim = 100,
                                   ...) {

  # Validate inputs
  if (!inherits(pipeline, "ml4t2e_pipeline")) {
    stop("'pipeline' must be an ml4t2e_pipeline object created with ml4t2e_fit_pipeline().",
         call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

  if (!is.numeric(time_horizon) || length(time_horizon) != 1 ||
      is.na(time_horizon) || time_horizon <= 0) {
    stop("'time_horizon' must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(nsim) || length(nsim) != 1 || nsim < 1) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }

  # Check for required packages
  if (!requireNamespace("fastshap", quietly = TRUE) &&
      !requireNamespace("kernelshap", quietly = TRUE)) {
    stop("Either 'fastshap' or 'kernelshap' package must be installed for SHAP calculation.\n",
         "Install with: install.packages('fastshap') or install.packages('kernelshap')",
         call. = FALSE)
  }

  # Extract feature columns (exclude time, event, and id variables)
  feature_cols <- pipeline$original_expvars

  # Ensure all features are present in data
  missing_features <- setdiff(feature_cols, colnames(data))
  if (length(missing_features) > 0) {
    stop("Missing features in data: ", paste(missing_features, collapse = ", "),
         call. = FALSE)
  }

  # Extract feature matrix
  X <- data[, feature_cols, drop = FALSE]

  # Set up background data if not provided
  if (is.null(background)) {
    # Use a random sample from the data as background
    n_background <- min(100, nrow(data))
    if (nrow(data) > n_background) {
      background_idx <- sample(nrow(data), n_background)
      background <- data[background_idx, feature_cols, drop = FALSE]
    } else {
      background <- X
    }
  } else {
    background <- background[, feature_cols, drop = FALSE]
  }

  # Create prediction wrapper
  pred_fn <- ml4t2e_shap_predict_fn(
    pipeline = pipeline,
    time_horizon = time_horizon
  )

  # Get baseline prediction (average over background)
  baseline_pred <- mean(pred_fn(background), na.rm = TRUE)

  # Get actual predictions for the data
  predictions <- pred_fn(data)

  # Calculate SHAP values
  message("Calculating SHAP values using ", nsim, " simulations...")

  if (requireNamespace("fastshap", quietly = TRUE)) {
    # Use fastshap if available
    shap_obj <- tryCatch({
      fastshap::explain(
        object = pred_fn,
        X = X,
        pred_wrapper = function(object, newdata) {
          object(newdata)
        },
        nsim = nsim,
        newdata = X,
        ...
      )
    }, error = function(e) {
      stop("Error calculating SHAP values with fastshap: ", e$message, call. = FALSE)
    })

    shap_values <- as.matrix(shap_obj)

  } else if (requireNamespace("kernelshap", quietly = TRUE)) {
    # Fall back to kernelshap
    shap_obj <- tryCatch({
      kernelshap::kernelshap(
        object = pred_fn,
        X = X,
        bg_X = background,
        ...
      )
    }, error = function(e) {
      stop("Error calculating SHAP values with kernelshap: ", e$message, call. = FALSE)
    })

    shap_values <- shap_obj$S
  } else {
    stop("Neither fastshap nor kernelshap could be loaded.", call. = FALSE)
  }

  # Ensure SHAP values have feature names
  if (is.null(colnames(shap_values))) {
    colnames(shap_values) <- feature_cols
  }

  # Create result object
  result <- list(
    shap_values = shap_values,
    baseline = baseline_pred,
    predictions = predictions,
    feature_values = as.data.frame(X),
    feature_names = feature_cols,
    time_horizon = time_horizon,
    analysis_type = pipeline$analysis_type,
    pipeline = pipeline
  )

  class(result) <- c("ml4t2e_shap", "list")

  message("SHAP calculation complete.")
  return(result)
}


#' @export
print.ml4t2e_shap <- function(x, ...) {
  cat("ml4time2event SHAP Analysis\n")
  cat("============================\n")
  cat("Analysis type  :", x$analysis_type, "\n")
  cat("Time horizon   :", x$time_horizon, "\n")
  cat("Observations   :", nrow(x$shap_values), "\n")
  cat("Features       :", ncol(x$shap_values), "\n")
  cat("Baseline pred  :", round(x$baseline, 3), "\n")
  cat("\nFeatures:\n")
  cat(paste(" ", x$feature_names, collapse = "\n"), "\n")
  invisible(x)
}
