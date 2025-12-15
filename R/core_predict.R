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
                            include = "all",
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
                    time_indices <- match(times, object$time_grid)
                    valid_t_idx <- !is.na(time_indices)
                    if (!all(valid_t_idx)) {
                        rlang::warn("Some requested times are not in the calibration time grid. Bands will be NA for those times.")
                    }

                    q_values <- rep(NA_real_, length(times))

                    for (k in seq_along(times)) {
                        if (valid_t_idx[k]) {
                            idx_grid <- time_indices[k]
                            s <- scores_obj$scores[, idx_grid]
                            w <- scores_obj$weights[, idx_grid]
                            q_values[k] <- ml4t2e_weighted_quantile(s, w, conformal_alpha)
                        }
                    }

                    q_map <- data.frame(time = times, q = q_values)
                    pred <- pred %>%
                        dplyr::left_join(q_map, by = "time") %>%
                        dplyr::mutate(
                            lower = pmax(0, surv - q),
                            upper = pmin(1, surv + q)
                        ) %>%
                        dplyr::select(-q)
                } else if (identical(outcome_type, "competing_risk") && identical(type, "cif")) {
                    q_df_list <- list()
                    for (cause_val in names(scores_obj$scores)) {
                        cause_scores <- scores_obj$scores[[cause_val]]
                        time_indices <- match(times, object$time_grid)
                        valid_t_idx <- !is.na(time_indices)

                        q_vals <- rep(NA_real_, length(times))
                        for (k in seq_along(times)) {
                            if (valid_t_idx[k]) {
                                idx_grid <- time_indices[k]
                                s <- cause_scores$scores[idx_grid, ]
                                w <- cause_scores$weights[idx_grid, ]
                                q_vals[k] <- ml4t2e_weighted_quantile(s, w, conformal_alpha)
                            }
                        }
                        q_df_list[[cause_val]] <- data.frame(time = times, cause = cause_val, q = q_vals, stringsAsFactors = FALSE)
                    }

                    q_map <- dplyr::bind_rows(q_df_list)

                    if (is.numeric(pred$cause) && !is.numeric(q_map$cause)) {
                        q_map$cause <- as.numeric(q_map$cause)
                    } else if (is.factor(pred$cause)) {
                        q_map$cause <- as.factor(q_map$cause)
                    }

                    pred <- pred %>%
                        dplyr::left_join(q_map, by = c("time", "cause")) %>%
                        dplyr::mutate(
                            lower = pmax(0, cif - q),
                            upper = pmin(1, cif + q)
                        ) %>%
                        dplyr::select(-q)
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
    if (identical(include, "all")) {
        return(models)
    }
    include <- unique(include)
    include <- include[include %in% models]
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
