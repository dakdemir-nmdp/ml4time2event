#' XGBoost Competing Risk Model (R6)
#'
#' Fits a cause-specific Cox model using XGBoost for each competing risk.
#'
#' @keywords internal
#' @noRd
XGBoostCompetingRiskModel <- R6::R6Class(
    classname = "XGBoostCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of {model, baseline_hazard} per cause
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL, # map cause name -> code
        varprof = NULL,
        feature_names = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            cause_map <- task$metadata$cause_map
            if (is.null(cause_map) || nrow(cause_map) == 0) {
                rlang::abort("Competing-risk tasks must contain a non-empty `cause_map` in metadata.")
            }
            self$cause_codes <- setNames(cause_map$code, as.character(cause_map$cause))
            self$varprof <- VariableProfile(data, task$features)

            spec <- self$spec
            nrounds <- spec$nrounds %||% 100
            eta <- spec$eta %||% 0.01
            max_depth <- spec$max_depth %||% 5
            verbose <- spec$verbose %||% 0

            # Model Matrix
            X <- stats::model.matrix(~ -1 + ., data[, task$features, drop = FALSE])
            self$feature_names <- colnames(X)

            models_list <- list()
            cause_labels <- names(self$cause_codes)

            time_col <- task$time_col
            event_code_col <- task$event_col # numeric codes 0, 1, 2...

            y_times <- data[[time_col]]
            y_codes <- data[[event_code_col]]

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]

                # Binary event for this cause: 1 if code matches, 0 otherwise (censored/competing)
                # XGBoost Cox: label = time if event, -time if censored

                is_event <- (y_codes == cause_code)
                y_label <- ifelse(is_event, y_times, -y_times)

                if (sum(is_event) < 1) {
                    rlang::warn(glue::glue("No events for cause {cause_label}. Skipping."))
                    next
                }

                dtrain <- xgboost::xgb.DMatrix(data = X, label = y_label)

                params <- list(
                    objective = "survival:cox",
                    eval_metric = "cox-nloglik",
                    eta = eta,
                    max_depth = max_depth
                )

                bst <- tryCatch(
                    {
                        xgboost::xgb.train(
                            params = params,
                            data = dtrain,
                            nrounds = nrounds,
                            verbose = verbose
                        )
                    },
                    error = function(e) {
                        rlang::warn(glue::glue("Failed to fit XGBoost for cause {cause_label}: {e$message}"))
                        NULL
                    }
                )

                if (!is.null(bst)) {
                    # Estimate Baseline Hazard (Breslow)
                    lp_train <- predict(bst, X)
                    df_train <- data.frame(
                        time = y_times,
                        status = as.numeric(is_event),
                        exp_lp = exp(lp_train)
                    )

                    unique_event_times <- sort(unique(df_train$time[df_train$status == 1]))

                    if (length(unique_event_times) == 0) {
                        bh <- data.frame(time = c(0, max(y_times)), hazard = c(0, 0.01))
                    } else {
                        # Vectorized Breslow
                        df_train <- df_train[order(df_train$time), ]

                        # Denominator: sum exp_lp where time >= t
                        risk_denominators <- vapply(unique_event_times, function(t) {
                            sum(df_train$exp_lp[df_train$time >= t])
                        }, numeric(1))

                        events_at_t <- as.vector(table(factor(df_train$time[df_train$status == 1], levels = unique_event_times)))

                        h0_vals <- ifelse(risk_denominators > 0, events_at_t / risk_denominators, 0)

                        bh <- data.frame(time = unique_event_times, hazard = h0_vals)

                        if (bh$time[1] != 0) {
                            bh <- rbind(data.frame(time = 0, hazard = 0), bh)
                        }
                    }

                    models_list[[cause_label]] <- list(
                        model = bst,
                        baseline_hazard = bh
                    )
                }
            }

            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()

            complete_data <- .ensure_prediction_data(newdata, self$task)

            # Prepare Matrix
            X_new <- stats::model.matrix(~ -1 + ., complete_data[, self$task$features, drop = FALSE])

            # Alignment
            missing_cols <- setdiff(self$feature_names, colnames(X_new))
            if (length(missing_cols) > 0) {
                match_mat <- matrix(0, nrow = nrow(X_new), ncol = length(missing_cols))
                colnames(match_mat) <- missing_cols
                X_new <- cbind(X_new, match_mat)
            }
            extra_cols <- setdiff(colnames(X_new), self$feature_names)
            if (length(extra_cols) > 0) {
                X_new <- X_new[, !colnames(X_new) %in% extra_cols, drop = FALSE]
            }
            X_new <- X_new[, self$feature_names, drop = FALSE]

            dtest <- xgboost::xgb.DMatrix(data = X_new) # Actually predict accepts matrix too usually

            # Determine Eval Times
            if (is.null(times)) {
                eval_times <- self$time_grid
            } else {
                eval_times <- sort(unique(c(0, as.numeric(times))))
            }

            n_obs <- nrow(complete_data)
            n_t <- length(eval_times)
            cause_labels <- names(self$cause_codes)

            # Compute H_k(t) for each cause
            hazards_list <- list()
            overall_cumhaz <- matrix(0, nrow = n_obs, ncol = n_t) # H_all(t) = sum H_k(t)

            for (cause in cause_labels) {
                obj <- self$model[[cause]]

                if (is.null(obj)) {
                    hazards_list[[cause]] <- matrix(0, nrow = n_obs, ncol = n_t)
                    next
                }

                bst <- obj$model
                bh <- obj$baseline_hazard

                lp <- predict(bst, X_new) # Vector n_obs
                exp_lp <- exp(lp)

                # H0(t) at eval_times
                bh$cumhaz <- cumsum(bh$hazard)
                H0_t <- stats::approx(bh$time, bh$cumhaz, xout = eval_times, method = "constant", rule = 2, f = 1)$y

                # H_k(t) = exp(lp) * H0(t)
                # Matrix: rows=obs, cols=times (Wait, match orientation)
                # Outer: exp_lp (n) * H0_t (m) -> n x m
                H_k_mat <- outer(exp_lp, H0_t) # [n_obs, n_t]

                hazards_list[[cause]] <- H_k_mat
                overall_cumhaz <- overall_cumhaz + H_k_mat
            }

            # Overall Survival S(t) = exp(-H_all(t))
            overall_surv <- exp(-overall_cumhaz) # [n_obs, n_t]

            # CIF Calculation (Aalen-Johansen)
            # CIF_k(t) = integral S(u-) dH_k(u)

            preds_list <- list()
            id_complete <- complete_data[[self$task$id_col]]

            # We need final requested times
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

            for (cause in cause_labels) {
                H_k <- hazards_list[[cause]] # [n_obs, n_t]

                # dH_k(t) = H_k(t_j) - H_k(t_{j-1})
                # Columns are time.
                # Diff along columns.

                # apply(H_k, 1, diff) returns [n_t-1, n_obs].
                # We want [n_obs, n_t]. First col is H_k(t1) (since t0=0, H(0)=0 usually)
                # If eval_times[1]==0, H(0)=0.

                # Helper: t(apply(...))
                dH_k <- t(apply(H_k, 1, function(row) c(row[1], diff(row))))

                # S(u-) = S(t_{j-1}) for integral at t_j
                # Shift S right
                S_prev <- cbind(rep(1, n_obs), overall_surv[, -ncol(overall_surv), drop = FALSE])

                # Integrate
                # apply(..., 1, cumsum) -> [n_t, n_obs] (rows=times)
                cif_k_t <- apply(S_prev * dH_k, 1, cumsum) # [times, obs]

                # Interpolate to req_times
                # cif_k_t is [times, obs]. Interpolator expects [times, obs]?
                # Check cifMatInterpolator arguments: probsMat (rows=times, cols=obs).
                # Yes.

                if (!identical(eval_times, req_times)) {
                    cif_interp <- cifMatInterpolator(cif_k_t, eval_times, req_times)
                    times_out <- req_times
                } else {
                    cif_interp <- cif_k_t
                    times_out <- eval_times
                }

                # Flatten
                new_cif_prediction(
                    id = rep(id_complete, each = length(times_out)),
                    time = rep(times_out, times = length(id_complete)),
                    cause = rep(cause, length(id_complete) * length(times_out)),
                    cif = as.vector(cif_interp),
                    model = rep("xgboost", length(id_complete) * length(times_out)),
                    ensemble = FALSE,
                    set = set
                ) -> pred_obj

                preds_list[[length(preds_list) + 1]] <- pred_obj
            }

            result <- dplyr::bind_rows(preds_list)
            result$cif <- pmin(pmax(result$cif, 0), 1)
            result
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "XGBoost Competing Risk"
            info
        },
        required_packages = function() {
            c("xgboost")
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

# Register
.register_time_to_event_model(
    engine = "cr_xgboost",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        XGBoostCompetingRiskModel$new(spec = modifyList(list(engine = "cr_xgboost"), spec))
    },
    packages = "xgboost",
    label = "XGBoost Competing Risk",
    tags = c("xgboost", "competing-risk")
)
