#' Evaluate fitted models on metrics
#'
#' Computes prediction performance metrics on a task. Currently supports:
#' - **Concordance Index (C-index)**: Harrell's C-statistic measuring discrimination.
#'   For survival: \eqn{C = \frac{\sum_{i,j} I(t_i < t_j) \cdot I(\hat{S}_i > \hat{S}_j) \cdot \delta_j}{\sum_{i,j} I(t_i < t_j) \cdot \delta_j}}
#'   where \eqn{\delta_j} is event indicator, \eqn{\hat{S}_i} is predicted survival probability.
#'   Range: 0.5 (random) to 1.0 (perfect).
#'
#' - **Integrated Brier Score (IBS)**: Mean squared difference between predicted
#'   and observed outcomes over time.
#'   \eqn{IBS = \frac{1}{t_{max} - t_{min}} \int_0^{t_{max}} BS(t) dt}
#'   where \eqn{BS(t) = E[(S(t|X) - Y(t))^2]} for survival data.
#'   Range: 0 (perfect) to 1 (worst).
#'
#' For competing risks, C-index is computed for each cause separately using
#' time-dependent concordance, and Brier score uses cumulative incidence functions.
#'
#' @param fit_or_preds A `t2e_fit` object or predictions returned by
#'   `predict.t2e_fit()`.
#' @param task Optional task providing ground truth if `fit_or_preds` contains
#'   predictions on new data.
#' @param metrics Character vector of metrics to compute. Default: `c("c_index", "ibs")`.
#' @param include Which models to include (same semantics as `predict()`).
#'
#' @return A tibble with columns `model`, `metric`, and `value`.
#'
#' @references
#' Harrell, F. E., Lee, K. L., & Mark, D. B. (1996).
#' "Multivariable prognostic models: Issues in developing models, evaluating
#' assumptions and adequacy, and measuring and reducing errors."
#' *Statistics in Medicine*, 15(4), 361–387.
#'
#' @keywords internal
#' @export
ml4t2e_evaluate <- function(fit_or_preds,
                            task = NULL,
                            metrics = c("c_index", "ibs"),
                            include = "all") {
    predictions <- fit_or_preds
    truth_task <- task

    if (inherits(fit_or_preds, "t2e_fit")) {
        outcome_type <- fit_or_preds$outcome_type
        pred_type <- if (identical(outcome_type, "competing_risk")) "cif" else "survival"

        # Check if task data is available
        if (is.null(fit_or_preds$task$data) && is.null(task)) {
            rlang::abort("Model fit does not contain training data. You must provide `task` (or `test data` task) for evaluation.")
        }

        # If using internal task, check if it has data. If not, error.
        # Note: predict() will check for data.
        predictions <- predict(fit_or_preds, include = include, type = pred_type)

        if (is.null(truth_task)) {
            truth_task <- fit_or_preds[["task"]]
        }
    }

    if (!inherits(predictions, "t2e_pred")) {
        rlang::abort("`fit_or_preds` must be a `t2e_fit` or predictions from `predict()`.")
    }
    if (is.null(truth_task)) {
        rlang::abort("`task` must be supplied when evaluating detached predictions.")
    }
    if (is.null(truth_task$data)) {
        rlang::abort("The supplied `task` must contain data.")
    }

    task_type <- attr(truth_task, "task_type")
    if (identical(task_type, "competing_risk")) {
        return(.evaluate_competing_risk(predictions, truth_task, metrics))
    }

    metrics <- unique(metrics)
    allowed_metrics <- c("c_index", "ibs")
    bad <- setdiff(metrics, allowed_metrics)
    if (length(bad) > 0) {
        rlang::abort(paste0("Unsupported metrics: ", paste(bad, collapse = ", ")))
    }

    task_df <- truth_task[["data"]]
    id_col <- truth_task[["id_col"]]
    time_col <- truth_task[["time_col"]]
    event_col <- truth_task[["event_col"]]

    results <- list()
    for (model_name in unique(predictions[["model"]])) {
        pred_subset <- predictions[predictions[["model"]] == model_name, , drop = FALSE]
        metric_values <- list()
        if ("c_index" %in% metrics) {
            lp <- .extract_etl_from_survival(pred_subset, task_df, id_col)
            valid_idx <- !is.na(lp)
            if (sum(valid_idx) == 0) {
                metric_values[["c_index"]] <- NA_real_
            } else {
                # Calculate concordance (Risk score 'lp' is inversely related to survival time, so use -lp)
                # We allow this to fail (return NA) rather than flipping it.
                c_val <- tryCatch(
                    {
                        # Create a safe data frame for evaluation to avoid scoping issues with formula
                        eval_df <- data.frame(
                            time = task_df[[time_col]][valid_idx],
                            status = task_df[[event_col]][valid_idx],
                            score = -1 * lp[valid_idx]
                        )
                        cc <- survival::concordance(survival::Surv(time, status) ~ score, data = eval_df)
                        unname(cc$concordance)
                    },
                    error = function(e) {
                        # Debug: print(e)
                        NA_real_
                    }
                )
                metric_values[["c_index"]] <- c_val
            }
        }
        if ("ibs" %in% metrics) {
            ibs <- .integrated_brier(pred_subset, task_df, id_col, time_col, event_col)
            metric_values[["ibs"]] <- ibs
        }
        for (metric_name in names(metric_values)) {
            results[[length(results) + 1]] <- data.frame(
                model = model_name,
                metric = metric_name,
                value = metric_values[[metric_name]],
                stringsAsFactors = FALSE
            )
        }
    }
    dplyr::bind_rows(results)
}

.extract_risk_from_predictions <- function(risk_predictions, task_df, id_col) {
    risk_scores <- numeric(nrow(task_df))
    risk_scores[] <- NA_real_

    id_map <- match(task_df[[id_col]], risk_predictions$id)
    valid_idx <- !is.na(id_map)
    risk_scores[valid_idx] <- risk_predictions$risk[id_map[valid_idx]]

    risk_scores
}

.extract_etl_from_survival <- function(predictions, task_df, id_col) {
    base_preds <- predictions[, c("id", "time", "surv")]
    # Ensure max_time and times are numeric
    unique_times <- sort(unique(base_preds$time))
    max_time <- max(unique_times, na.rm = TRUE)

    # Use dplyr for vectorized grouping
    # Calculate ETL for each ID
    etl_df <- base_preds %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(
            etl = calculate_expected_time_lost(
                times = time,
                event_probs = 1 - surv,
                upper_limit = max_time,
                lower_limit = 0
            ),
            .groups = "drop"
        )

    # Join back to ensure alignment with task_df
    # Convert IDs to character for safe matching if needed, or rely on type match
    task_ids <- data.frame(id = task_df[[id_col]])
    # Ensure type compatibility
    if (is.character(etl_df$id) != is.character(task_ids$id)) {
        etl_df$id <- as.character(etl_df$id)
        task_ids$id <- as.character(task_ids$id)
    }

    result_df <- dplyr::left_join(task_ids, etl_df, by = "id")
    return(result_df$etl)
}

.extract_etl_from_cif <- function(predictions, task_df, id_col, cause_label = NULL) {
    if (!is.null(cause_label)) {
        base_preds <- predictions[predictions$cause == cause_label, c("id", "time", "cif")]
    } else {
        base_preds <- predictions[, c("id", "time", "cif")]
    }

    if (nrow(base_preds) == 0) {
        return(rep(NA_real_, nrow(task_df)))
    }

    unique_times <- sort(unique(base_preds$time))
    max_time <- max(unique_times, na.rm = TRUE)

    ids <- task_df[[id_col]]
    unique_ids <- unique(ids)

    etl_scores <- numeric(length(unique_ids))
    names(etl_scores) <- unique_ids

    for (id_val in unique_ids) {
        id_preds <- base_preds[base_preds$id == id_val, ]
        if (nrow(id_preds) == 0) {
            etl_scores[as.character(id_val)] <- NA_real_
            next
        }

        id_preds <- id_preds[order(id_preds$time), ]
        cif_curve <- id_preds$cif
        time_curve <- id_preds$time

        cif_curve <- pmax(0, pmin(1, cif_curve))

        etl <- calculate_expected_time_lost(
            times = time_curve,
            event_probs = cif_curve,
            upper_limit = max_time,
            lower_limit = 0
        )

        etl_scores[as.character(id_val)] <- etl
    }

    risk_scores <- numeric(nrow(task_df))
    risk_scores[] <- NA_real_

    id_map <- match(as.character(ids), names(etl_scores))
    valid_idx <- !is.na(id_map)
    risk_scores[valid_idx] <- etl_scores[id_map[valid_idx]]

    risk_scores
}

.extract_linear_predictor <- function(predictions, task_df, id_col, time_col, event_col = NULL) {
    if ("risk" %in% colnames(predictions)) {
        return(.extract_risk_from_predictions(predictions, task_df, id_col))
    }

    if ("surv" %in% colnames(predictions)) {
        return(.extract_etl_from_survival(predictions, task_df, id_col))
    }

    rlang::abort("Unable to extract risk scores from predictions.")
}

.integrated_brier <- function(predictions, task_df, id_col, time_col, event_col) {
    times <- sort(unique(predictions$time))
    base_preds <- as.data.frame(predictions[, c("id", "time", "surv")])

    # Reshape to wide: One row per ID, one column per time
    surv_mat_wide <- base_preds %>%
        tidyr::pivot_wider(
            id_cols = "id",
            names_from = "time",
            values_from = "surv"
        )

    # Ensure alignment with task_df IDs
    id_map <- match(task_df[[id_col]], surv_mat_wide$id)
    surv_mat_wide <- surv_mat_wide[id_map, , drop = FALSE]

    # Remove ID column
    # Identify numeric columns corresponding to times
    # (reshape keeps id column)
    surv_mat <- as.matrix(surv_mat_wide[, -1, drop = FALSE])

    # BrierScore expects: rows = times, cols = observations
    # Current: rows = observations, cols = times
    # Transpose
    surv_mat_t <- t(surv_mat)

    # Extract observed data
    obstimes <- task_df[[time_col]]
    obsevents <- task_df[[event_col]]

    # Call the IPCW-corrected Integrated Brier Score
    # We use the predicted times as evaluation times for integration
    integratedBrier(
        predsurv = surv_mat_t,
        pred_times = times,
        obstimes = obstimes,
        obsevents = obsevents,
        eval_times = times
    )
}

.evaluate_competing_risk <- function(predictions, task, metrics) {
    metrics <- unique(metrics)
    allowed_metrics <- c("c_index", "ibs")
    bad <- setdiff(metrics, allowed_metrics)
    if (length(bad) > 0) {
        rlang::abort(paste0("Unsupported metrics: ", paste(bad, collapse = ", ")))
    }

    task_df <- task[["data"]]
    id_col <- task[["id_col"]]
    time_col <- task[["time_col"]]
    event_col <- task[["event_col"]]
    cause_map <- task[["metadata"]][["cause_map"]]
    if (is.null(cause_map) || nrow(cause_map) == 0) {
        rlang::abort("Competing risk task metadata must include a cause map.")
    }

    ids <- task_df[[id_col]]
    time_grid <- attr(predictions, "time_grid")
    if (is.null(time_grid)) {
        time_grid <- sort(unique(predictions$time))
    }
    surv_obj <- survival::Surv(task_df[[time_col]], task_df[[event_col]], type = "mstate")

    results <- list()
    model_names <- unique(predictions$model)
    for (model_name in model_names) {
        model_pred <- predictions[predictions$model == model_name, , drop = FALSE]
        for (i in seq_len(nrow(cause_map))) {
            cause_label <- cause_map$cause[i]
            cause_code <- cause_map$code[i]
            cause_pred <- model_pred[model_pred$cause == cause_label, , drop = FALSE]
            if (nrow(cause_pred) == 0) {
                next
            }
            matrix <- .cif_prediction_matrix(cause_pred, ids, time_grid)

            if ("ibs" %in% metrics) {
                ibs_val <- integratedBrierCR(
                    SurvObj = surv_obj,
                    Predictions = matrix,
                    eval_times = time_grid,
                    cause = cause_code,
                    pred_times = time_grid
                )
                results[[length(results) + 1]] <- data.frame(
                    model = model_name,
                    cause = cause_label,
                    metric = "ibs",
                    value = ibs_val,
                    stringsAsFactors = FALSE
                )
            }
            if ("c_index" %in% metrics) {
                cause_pred_subset <- cause_pred[cause_pred$cause == cause_label, ]
                if (nrow(cause_pred_subset) > 0) {
                    etl_scores <- .extract_etl_from_cif(
                        predictions = cause_pred_subset,
                        task_df = task_df,
                        id_col = id_col,
                        cause_label = cause_label
                    )
                    event_code_col <- ".event_code"
                    if (event_code_col %in% colnames(task_df)) {
                        event_indicator <- ifelse(task_df[[event_code_col]] == cause_code, 1L, 0L)
                    } else {
                        event_indicator <- ifelse(task_df[[event_col]] == cause_code, 1L, 0L)
                    }

                    valid_idx <- !is.na(etl_scores) & is.finite(etl_scores)
                    if (sum(valid_idx) > 0 && sum(event_indicator[valid_idx]) > 0) {
                        surv_obj_cause <- survival::Surv(task_df[[time_col]], event_indicator)
                        c_val <- tryCatch(
                            {
                                concordance <- survival::concordance(surv_obj_cause[valid_idx] ~ -etl_scores[valid_idx])
                                unname(concordance$concordance)
                            },
                            error = function(e) {
                                concordance <- survival::concordance(surv_obj_cause[valid_idx] ~ etl_scores[valid_idx])
                                1 - unname(concordance$concordance)
                            }
                        )
                    } else {
                        c_val <- NA_real_
                    }
                } else {
                    c_val <- NA_real_
                }
                results[[length(results) + 1]] <- data.frame(
                    model = model_name,
                    cause = cause_label,
                    metric = "c_index",
                    value = c_val,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(results) == 0) {
        rlang::abort("No competing risk predictions available for evaluation.")
    }
    dplyr::bind_rows(results)
}

.cif_prediction_matrix <- function(pred_subset, ids, time_grid) {
    pred_clean <- pred_subset %>%
        dplyr::group_by(id, time) %>%
        dplyr::summarise(cif = mean(cif, na.rm = TRUE), .groups = "drop")

    mat <- matrix(NA_real_, nrow = length(time_grid), ncol = length(ids))
    colnames(mat) <- ids
    rownames(mat) <- time_grid

    id_index <- match(pred_clean$id, ids)
    time_index <- match(pred_clean$time, time_grid)
    valid <- !is.na(id_index) & !is.na(time_index)
    mat[cbind(time_index[valid], id_index[valid])] <- pred_clean$cif[valid]
    mat
}
