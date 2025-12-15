#' @importFrom stats predict
#' @importFrom utils modifyList
NULL

#' Make predictions from a fitted time-to-event model
#'
#' @param object A `t2e_fit` object.
#' @param newdata Optional data frame for prediction. If missing, attempts to use
#'   training data if cached.
#' @param times Numeric vector of times to predict. If `NULL`, defaults to the
#'   time grid used during training.
#' @param type type of prediction:
#'   * `"survival"`: Survival probabilities (default for survival tasks).
#'   * `"risk"`: Risk score (e.g., log-hazard).
#'   * `"cif"`: Cumulative incidence (default for competing risk tasks).
#' @param include Character vector of models to include. `"all"` (default)
#'   includes all individual models + ensemble.
#' @param conformal_alpha Optional alpha for conformal intervals (e.g. 0.05 for 95%).
#' @param ... Additional arguments passed to methods.
#'
#' @return A tidy tibble of predictions (class `t2e_pred`).
#' @export
predict.t2e_fit <- function(object,
                            newdata = NULL,
                            times = NULL,
                            type = NULL,
                            include = "fit_default",
                            conformal_alpha = NULL,
                            ...) {
    if (!inherits(object, "t2e_fit")) {
        rlang::abort("`object` must be a `t2e_fit`.")
    }
    outcome_type <- object[["outcome_type"]]
    if (is.null(type)) {
        type <- if (identical(outcome_type, "survival")) "survival" else "cif"
    }
    type_choices <- if (identical(outcome_type, "survival")) {
        c("survival", "risk")
    } else {
        c("cif")
    }
    type <- match.arg(type, type_choices)

    data <- newdata
    set_label <- "test"
    if (is.null(data)) {
        task <- object[["task"]]
        data <- task$data
        if (is.null(data)) {
            rlang::abort("Training data was not saved with the model. Please provide `newdata`.")
        }
        set_label <- "train"
    } else {
        data <- dplyr::as_tibble(data)
        .verify_new_data_columns(data, object[["task"]])
    }

    time_grid <- object[["time_grid"]]
    times <- .resolve_prediction_times(times, time_grid)

    include <- .normalize_include(include, object[["model_names"]])
    model_store <- object[["models"]]
    predictions <- list()

    # Check for conformal request
    compute_bands <- !is.null(conformal_alpha)
    if (compute_bands) {
        if (is.null(object$conformal)) {
            rlang::warn("`conformal_alpha` provided but model was not fitted with `conformal_calibration`. Bands will not be computed.")
            compute_bands <- FALSE
        } else {
            if (!is.numeric(conformal_alpha) || conformal_alpha <= 0 || conformal_alpha >= 1) {
                rlang::abort("`conformal_alpha` must be between 0 and 1.")
            }
        }
    }

    for (engine in include) {
        if (identical(engine, "ensemble")) {
            ensemble_obj <- object$ensemble
            if (is.null(ensemble_obj)) {
                rlang::warn("No ensemble available; skipping ensemble predictions.")
                next
            }
            if (identical(outcome_type, "survival")) {
                if (identical(type, "survival")) {
                    pred <- ensemble_obj$predict_survival(newdata = data, times = times, set = set_label)
                } else {
                    pred <- ensemble_obj$predict_risk(newdata = data, times = times, set = set_label)
                }
            } else {
                pred <- ensemble_obj$predict_cif(newdata = data, times = times, set = set_label)
            }
        } else {
            model_obj <- model_store[[engine]]
            if (is.null(model_obj)) {
                rlang::warn(glue::glue("Model '{engine}' not found; skipping."))
                next
            }
            if (identical(outcome_type, "survival")) {
                if (identical(type, "survival")) {
                    pred <- model_obj$predict_survival(newdata = data, times = times, set = set_label)
                } else {
                    pred <- model_obj$predict_risk(newdata = data, times = times, set = set_label)
                }
            } else {
                pred <- model_obj$predict_cif(newdata = data, times = times, set = set_label)
            }
        }

        # Compute bands if requested
        if (compute_bands) {
            scores_obj <- object$conformal$scores[[engine]]
            if (!is.null(scores_obj)) {
                if (identical(outcome_type, "survival") && identical(type, "survival")) {
                    # Pre-compute quantiles for all grid points
                    grid_times <- object$time_grid
                    n_grid <- length(grid_times)
                    q_grid <- numeric(n_grid)

                    found_scores <- FALSE
                    if (!is.null(scores_obj$scores) && !is.null(scores_obj$weights)) {
                        found_scores <- TRUE
                        for (g_idx in seq_len(n_grid)) {
                            s <- scores_obj$scores[, g_idx]
                            w <- scores_obj$weights[, g_idx]
                            q_grid[g_idx] <- ml4t2e_weighted_quantile(s, w, conformal_alpha)
                        }
                    }

                    if (found_scores) {
                        # Interpolate to requested times
                        # Use step-function (constant) or linear?
                        # Conformal scores Q(t) are often smooth-ish. Linear is reasonable.
                        # Use rule=2 to clamp to range (carry forward/backward)
                        q_interp <- stats::approx(grid_times, q_grid, xout = times, rule = 2)$y

                        q_map <- data.frame(time = times, q = q_interp)

                        pred <- pred %>%
                            dplyr::left_join(q_map, by = "time") %>%
                            dplyr::mutate(
                                lower = pmax(0, surv - q),
                                upper = pmin(1, surv + q)
                            ) %>%
                            dplyr::select(-q)
                    } else {
                        rlang::warn(glue::glue("No scores found for engine {engine}"))
                    }
                } else if (identical(outcome_type, "competing_risk") && identical(type, "cif")) {
                    # Competing Risks Interpolation
                    grid_times <- object$time_grid
                    n_grid <- length(grid_times)

                    q_df_list <- list()

                    # scores_obj$scores is a list of causes? No, scores_obj IS a list wrapper:
                    # scores_obj = list(scores = list("1" = ..., "2" = ...))
                    # Wait, look at ml4t2e_calibrate in conformal_prediction.R
                    # scores_list[[cause_val]] <- .compute_scores_core(...)
                    # scores <- list(scores = scores_list)
                    # So scores_obj$scores is a list of cause-specific score objects.

                    for (cause_val in names(scores_obj$scores)) {
                        cause_scores_struct <- scores_obj$scores[[cause_val]]
                        # cause_scores_struct has $scores (matrix) and $weights (matrix)

                        if (!is.null(cause_scores_struct$scores)) {
                            q_grid <- numeric(n_grid)
                            for (g_idx in seq_len(n_grid)) {
                                s <- cause_scores_struct$scores[, g_idx]
                                w <- cause_scores_struct$weights[, g_idx]
                                q_grid[g_idx] <- ml4t2e_weighted_quantile(s, w, conformal_alpha)
                            }

                            q_interp <- stats::approx(grid_times, q_grid, xout = times, rule = 2)$y
                            q_df_list[[cause_val]] <- data.frame(time = times, cause = cause_val, q = q_interp, stringsAsFactors = FALSE)
                        }
                    }

                    if (length(q_df_list) > 0) {
                        q_map <- dplyr::bind_rows(q_df_list)

                        if (is.numeric(pred$cause) && !is.numeric(q_map$cause)) {
                            q_map$cause <- as.numeric(q_map$cause)
                        } else if (is.factor(pred$cause)) {
                            # Need to match factor levels or convert to char for joining
                            # Safest to join on character
                            pred$cause_join <- as.character(pred$cause)
                            q_map$cause_join <- as.character(q_map$cause)
                        } else {
                            pred$cause_join <- as.character(pred$cause)
                            q_map$cause_join <- as.character(q_map$cause)
                        }

                        # Fix lint
                        cause_join <- NULL

                        # Join
                        pred <- pred %>%
                            dplyr::left_join(q_map %>% dplyr::select(-cause), by = c("time", "cause_join")) %>%
                            dplyr::mutate(
                                lower = pmax(0, cif - q),
                                upper = pmin(1, cif + q)
                            ) %>%
                            dplyr::select(-q, -cause_join)
                    }
                }
            }
        }

        predictions[[engine]] <- pred
    }

    if (length(predictions) == 0) {
        rlang::abort("No predictions were generated.")
    }

    result <- dplyr::bind_rows(predictions)

    # Clamp values to [0, 1] to ensure validity
    if ("cif" %in% colnames(result)) {
        result$cif <- pmax(0, pmin(1, result$cif))
    }
    if ("surv" %in% colnames(result)) {
        result$surv <- pmax(0, pmin(1, result$surv))
    }

    attr(result, "pred_type") <- type
    attr(result, "time_grid") <- times
    class(result) <- unique(c("t2e_pred", class(result)))
    result
}

.normalize_include <- function(include, models) {
    # Default behavior: If 'fit_default' or missing/NULL, prioritize ensemble
    if (length(include) == 1 && include == "fit_default") {
        if ("ensemble" %in% models) {
            return("ensemble")
        } else {
            return(models)
        }
    }

    if (identical(include, "all")) {
        return(models)
    }
    include <- unique(include)
    valid <- include %in% models
    if (!all(valid)) {
        invalid <- include[!valid]
        rlang::warn(paste0("Ignoring unknown models in include: ", paste(invalid, collapse = ", ")))
        include <- include[valid]
    }
    if (length(include) == 0) {
        rlang::abort("`include` did not match any fitted models.")
    }
    include
}

.resolve_prediction_times <- function(times, default) {
    if (is.null(times)) {
        return(default)
    }
    times <- sort(unique(as.numeric(times)))
    if (anyNA(times)) {
        rlang::abort("`times` must be numeric without NA values.")
    }
    times
}

.verify_new_data_columns <- function(newdata, task) {
    required <- unique(task$features)
    missing <- setdiff(required, colnames(newdata))
    if (length(missing) > 0) {
        rlang::abort(paste0("`newdata` is missing columns: ", paste(missing, collapse = ", ")))
    }
}
