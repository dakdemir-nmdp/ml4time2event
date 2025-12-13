#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @importFrom rlang sym
#' @importFrom tidyr pivot_wider
ml4t2e_shap_predict_fn <- function(object,
                                   time_horizon,
                                   ensemble_method = "average") {
  if (inherits(object, "T2EPipeline")) {
    pipeline <- object
    outcome_type <- pipeline$outcome$type
    fit_object <- pipeline$fit_object
    model_names <- pipeline$models
  } else if (inherits(object, "t2e_fit")) {
    pipeline <- NULL
    outcome_type <- object$outcome_type
    fit_object <- object
    model_names <- object$model_names
  } else {
    stop("'object' must be a `T2EPipeline` or `t2e_fit` object.", call. = FALSE)
  }

  if (!is.numeric(time_horizon) || length(time_horizon) != 1 ||
    is.na(time_horizon) || time_horizon <= 0) {
    stop("'time_horizon' must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(fit_object)) {
    stop("Object must be fitted before SHAP predictions can be generated.", call. = FALSE)
  }

  base_grid <- fit_object$time_grid %||% numeric(0)
  prediction_times <- sort(unique(c(0, base_grid, time_horizon)))
  model_candidates <- fit_object$model_names %||% character(0)
  preferred_model <- if ("ensemble" %in% model_candidates) {
    "ensemble"
  } else if (length(model_names) > 0) {
    model_names[1]
  } else if (length(model_candidates) > 0) {
    model_candidates[1]
  } else {
    stop("No fitted models available.", call. = FALSE)
  }

  function(newdata) {
    if (!is.data.frame(newdata)) {
      stop("Input to prediction function must be a data.frame.", call. = FALSE)
    }

    if (!is.null(pipeline)) {
      processed <- .pipeline_process_new_data(
        pipeline = pipeline,
        newdata = newdata,
        require_outcomes = FALSE
      )
    } else {
      processed <- newdata
    }

    prediction_type <- if (identical(outcome_type, "survival")) "survival" else "cif"
    pred_df <- tryCatch(
      {
        predict(
          fit_object,
          newdata = processed,
          times = prediction_times,
          type = prediction_type,
          include = "all"
        )
      },
      error = function(e) {
        stop("Error generating predictions for SHAP calculation: ", e$message, call. = FALSE)
      }
    )

    pred_df <- pred_df[pred_df$set != "train", , drop = FALSE]
    if (nrow(pred_df) == 0) {
      return(rep(NA_real_, nrow(processed)))
    }

    available_models <- unique(pred_df$model)
    model_choice <- if (preferred_model %in% available_models) {
      preferred_model
    } else {
      available_models[1]
    }

    value_col <- if (identical(outcome_type, "survival")) "surv" else "cif"
    pred_selected <- pred_df[pred_df$model == model_choice, , drop = FALSE]
    if (nrow(pred_selected) == 0) {
      return(rep(NA_real_, nrow(processed)))
    }

    val_sym <- rlang::sym(value_col)

    if (!identical(outcome_type, "survival")) {
      pred_selected <- pred_selected |>
        dplyr::group_by(id, time) |>
        dplyr::summarise(
          !!value_col := sum(!!val_sym, na.rm = TRUE),
          .groups = "drop"
        )
    }

    id_levels <- unique(pred_selected$id)
    pred_selected <- pred_selected |>
      dplyr::mutate(id = factor(id, levels = id_levels)) |>
      dplyr::arrange(time, id)

    wide <- pred_selected |>
      dplyr::select(time, id, value = !!val_sym) |>
      tidyr::pivot_wider(names_from = id, values_from = value) |>
      dplyr::arrange(time)

    time_vec <- wide$time
    value_mat <- as.matrix(wide[, -1, drop = FALSE])

    if (identical(outcome_type, "survival")) {
      value_mat[is.na(value_mat)] <- 1
      event_probs <- 1 - value_mat
    } else {
      value_mat[is.na(value_mat)] <- 0
      value_mat[value_mat < 0] <- 0
      value_mat[value_mat > 1] <- 1
      event_probs <- value_mat
    }

    time_lost <- apply(event_probs, 2, function(obs_probs) {
      Integrator(
        times = time_vec,
        scores = obs_probs,
        minmax = c(0, time_horizon),
        scale = FALSE
      )
    })

    ordered_cols <- colnames(event_probs)
    as.numeric(time_lost[ordered_cols])
  }
}


#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' Calculate SHAP values for ml4time2event pipelines
#'
#' Generates expected-time-lost SHAP values for a fitted `ml4time2event`
#' pipeline by wrapping the tidy survival predictions in a prediction function
#' suitable for `fastshap` or `kernelshap`.
#'
#' @param pipeline A fitted pipeline produced by `ml4t2e_fit_pipeline()`.
#' @param data Data frame of observations to explain.
#' @param time_horizon Single positive numeric time horizon used for expected
#'   time lost calculations.
#' @param background Optional background data frame; defaults to a sample from
#'   `data`.
#' @param nsim Number of Monte Carlo simulations passed to `fastshap::explain()`
#'   when available.
#' @param ... Additional arguments forwarded to the underlying SHAP engine.
#'
#' @return An object of class `ml4t2e_shap` containing SHAP values, baseline
#'   prediction, raw predictions, and feature values.
#' @export
ml4t2e_calculate_shap <- function(object,
                                  data,
                                  time_horizon,
                                  background = NULL,
                                  nsim = 100,
                                  ...) {
  # Determine object type and properties
  if (inherits(object, "T2EPipeline")) {
    pipeline <- object
    feature_cols <- pipeline$features %||% pipeline$outcome$features
    analysis_type <- pipeline$outcome$type
  } else if (inherits(object, "t2e_fit")) {
    pipeline <- NULL
    feature_cols <- object$task$features
    analysis_type <- object$outcome_type
  } else {
    stop("'object' must be a `T2EPipeline` or `t2e_fit` object.", call. = FALSE)
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
      call. = FALSE
    )
  }

  # Ensure all features are present in data
  missing_features <- setdiff(feature_cols, colnames(data))
  if (length(missing_features) > 0) {
    stop("Missing features in data: ", paste(missing_features, collapse = ", "),
      call. = FALSE
    )
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
    object = object,
    time_horizon = time_horizon
  )

  # Get baseline prediction (average over background)
  baseline_pred <- mean(pred_fn(background), na.rm = TRUE)

  # Get actual predictions for the data
  predictions <- pred_fn(data)

  # Calculate SHAP values
  message("Calculating SHAP values using ", nsim, " simulations...")

  if (requireNamespace("kernelshap", quietly = TRUE)) {
    # Use kernelshap if available (preferred for consistency)
    shap_obj <- tryCatch(
      {
        kernelshap::kernelshap(
          object = object,
          X = X,
          bg_X = background,
          pred_fun = function(m, d) pred_fn(d),
          ...
        )
      },
      error = function(e) {
        stop("Error calculating SHAP values with kernelshap: ", e$message, call. = FALSE)
      }
    )

    shap_values <- shap_obj$S
  } else if (requireNamespace("fastshap", quietly = TRUE)) {
    # Fall back to fastshap
    shap_obj <- tryCatch(
      {
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
      },
      error = function(e) {
        stop("Error calculating SHAP values with fastshap: ", e$message, call. = FALSE)
      }
    )

    shap_values <- as.matrix(shap_obj)
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
    analysis_type = analysis_type,
    source_object = object
  )

  class(result) <- c("ml4t2e_shap", "list")

  message("SHAP calculation complete.")
  return(result)
}


#' @noRd
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
