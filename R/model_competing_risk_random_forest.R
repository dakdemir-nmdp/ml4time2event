#' Random Forest Model for Competing Risks (R6)
#'
#' Fits a Random Forest for competing risks, potentially tuning for each cause separately.
#'
#' @keywords internal
#' @noRd
RandomForestCompetingRiskModel <- R6::R6Class(
    classname = "RandomForestCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of models, one per cause
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL,
        varprof = NULL,
        factor_levels = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)

            data <- as.data.frame(task$data)
            cause_map <- task$metadata$cause_map
            if (is.null(cause_map) || nrow(cause_map) == 0) {
                rlang::abort("Competing-risk tasks must contain a non-empty `cause_map` in metadata.")
            }

            self$task <- task
            self$time_grid <- time_grid
            self$cause_codes <- setNames(cause_map$code, as.character(cause_map$cause))
            self$varprof <- VariableProfile(data, task$features)

            spec <- self$spec
            # Extract control parameters
            ntree <- spec$ntree %||% 300
            nodesize <- spec$nodesize
            mtry <- spec$mtry
            splitrule <- spec$splitrule %||% "logrankCR"

            cause_labels <- names(self$cause_codes)
            models_list <- list()

            timevar <- task$time_col
            eventcodevar <- task$event_col # 0, 1, 2...
            features <- task$features

            # Prepare data for rfsrc
            # rfsrc requires numeric/integer event code
            # and factors for character columns
            # We must subset data to only relevant columns for formula usage
            data_fit <- data[, c(timevar, eventcodevar, features), drop = FALSE]

            for (col in features) {
                if (is.character(data_fit[[col]])) {
                    data_fit[[col]] <- as.factor(data_fit[[col]])
                }
            }

            self$factor_levels <- lapply(data_fit[, features, drop = FALSE], function(x) {
                if (is.factor(x)) levels(x) else NULL
            })

            # Using standard Surv formula with data_fit
            # Note: variable names in formula must match data
            # Use `get(timevar)` style or constructing formula string carefully.
            # rfsrc typically handles string formulas well if syntax is clean.
            # "Surv(time, event) ~ ."
            form_str_rf <- paste0("Surv(", timevar, ", ", eventcodevar, ") ~ .")
            form_rf <- stats::as.formula(form_str_rf)

            unique_events <- sort(unique(data[[eventcodevar]][data[[eventcodevar]] != 0]))

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]
                cause_numeric <- as.numeric(cause_code)

                cause_arg <- c(cause_numeric, setdiff(unique_events, cause_numeric))

                run_mtry <- mtry
                if (is.null(run_mtry)) {
                    run_mtry <- max(1, floor(sqrt(length(features))))
                }
                run_nodesize <- nodesize
                if (is.null(run_nodesize)) {
                    run_nodesize <- 15
                }

                model_obj <- tryCatch(
                    {
                        randomForestSRC::rfsrc(
                            formula = form_rf,
                            data = data_fit,
                            ntree = as.integer(ntree),
                            nodesize = as.integer(run_nodesize),
                            mtry = as.integer(run_mtry),
                            splitrule = splitrule,
                            cause = cause_arg,
                            forest = TRUE,
                            importance = FALSE,
                            statistics = TRUE
                        )
                    },
                    error = function(e) {
                        # Try removing survival:: namespace or pre-loading logic if needed?
                        # Usually rfsrc needs Surv available.
                        rlang::warn(glue::glue("Failed to fit RF for cause {cause_label}: {e$message}"))
                        NULL
                    }
                )

                if (!is.null(model_obj)) {
                    models_list[[cause_label]] <- model_obj
                }
            }

            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()

            complete_data <- .ensure_prediction_data(newdata, self$task)

            if (is.null(times)) {
                target_times <- self$time_grid
            } else {
                target_times <- sort(unique(as.numeric(times)))
            }

            # Handle factors in newdata
            # (Use basic matching, or rely on rfsrc handling)
            # rfsrc is picky, best to ensure levels match training.
            # self$varprof or manual level check would be good.
            # Simplistic approach: cast chars to factor
            for (col in self$task$features) {
                if (is.character(complete_data[[col]])) {
                    complete_data[[col]] <- as.factor(complete_data[[col]])
                }
                levels_train <- self$factor_levels[[col]]
                if (!is.null(levels_train) && col %in% colnames(complete_data)) {
                    complete_data[[col]] <- factor(as.character(complete_data[[col]]), levels = levels_train)
                }
            }

            feature_frame <- complete_data[, self$task$features, drop = FALSE]
            complete_idx <- stats::complete.cases(feature_frame)

            ids_all <- complete_data[[self$task$id_col]]
            ids_complete <- ids_all[complete_idx]
            pred_data <- complete_data[complete_idx, , drop = FALSE]

            cause_labels <- names(self$cause_codes)
            preds_list <- list()

            for (cause_label in cause_labels) {
                model <- self$model[[cause_label]]
                cause_code <- self$cause_codes[[cause_label]]
                cif_complete <- matrix(
                    NA_real_,
                    nrow = length(ids_complete),
                    ncol = length(target_times)
                )

                if (!is.null(model) && nrow(pred_data) > 0) {
                    pred_obj <- randomForestSRC::predict.rfsrc(model, newdata = pred_data)

                    n_pred <- dim(pred_obj$cif)[1]
                    if (!is.null(n_pred) && n_pred > 0) {

                        # Extract CIF for specific cause
                        # pred_obj$cif is [obs, times, cause]
                        # Cause index:
                        # rfsrc usually names causes as "CIF.1", "CIF.2".
                        # We need to find the column corresponding to cause_code

                        cif_col_name <- paste0("CIF.", cause_code)
                        if (!cif_col_name %in% dimnames(pred_obj$cif)[[3]]) {
                            # Fallback: if cause dimension names are unavailable,
                            # return NA for this cause.
                            cif_vals <- matrix(NA, nrow = n_pred, ncol = length(pred_obj$time.interest))
                        } else {
                            cif_vals <- pred_obj$cif[, , cif_col_name] # [obs, times]
                        }

                        n_use <- min(nrow(cif_vals), length(ids_complete))
                        cif_vals <- cif_vals[seq_len(n_use), , drop = FALSE]

                        # Interpolate to target_times
                        base_times <- c(0, pred_obj$time.interest)
                        # Add 0 to cif_vals (time 0 = 0)
                        cif_vals <- cbind(0, cif_vals) # now [obs, n_times+1]

                        # We need [times, obs] for interpolator
                        cif_vals_t <- t(cif_vals) # [times, obs]

                        cif_interp <- cifMatInterpolator(cif_vals_t, base_times, target_times)
                        cif_complete[seq_len(n_use), ] <- t(cif_interp[, seq_len(n_use), drop = FALSE])
                    }
                }

                cif_vec <- as.vector(t(cif_complete))

                preds_complete <- new_cif_prediction(
                    id = rep(ids_complete, each = length(target_times)),
                    time = rep(target_times, times = length(ids_complete)),
                    cause = rep(cause_label, length(ids_complete) * length(target_times)),
                    cif = cif_vec,
                    model = rep("random_forest", length(ids_complete) * length(target_times)),
                    ensemble = FALSE,
                    set = set
                )

                if (!all(complete_idx)) {
                    missing_ids <- ids_all[!complete_idx]
                    preds_missing <- new_cif_prediction(
                        id = rep(missing_ids, each = length(target_times)),
                        time = rep(target_times, times = length(missing_ids)),
                        cause = rep(cause_label, length(missing_ids) * length(target_times)),
                        cif = rep(NA_real_, length(missing_ids) * length(target_times)),
                        model = rep("random_forest", length(missing_ids) * length(target_times)),
                        ensemble = FALSE,
                        set = set
                    )
                    preds_list[[length(preds_list) + 1]] <- dplyr::bind_rows(preds_complete, preds_missing)
                } else {
                    preds_list[[length(preds_list) + 1]] <- preds_complete
                }
            }

            result <- dplyr::bind_rows(preds_list)
            # Clamp
            result$cif <- pmin(pmax(result$cif, 0), 1)
            result
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "Random Forest Competing Risk"
            info
        },
        required_packages = function() {
            c("randomForestSRC")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) {
                rlang::abort("Model must be fitted before predictions can be generated.")
            }
        }
    )
)

# Register the model
.register_time_to_event_model(
    engine = "cr_random_forest",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        RandomForestCompetingRiskModel$new(spec = modifyList(list(engine = "cr_random_forest"), spec))
    },
    packages = "randomForestSRC",
    label = "Random Forest Competing Risk",
    tags = c("random-forest", "competing-risk")
)
