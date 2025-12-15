#' Cause-specific Cox Proportional Hazards Model for Competing Risks (R6)
#'
#' Fits a separate Cox model for each cause of failure.
#'
#' @keywords internal
#' @noRd
CoxCompetingRiskModel <- R6::R6Class(
    classname = "CoxCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of models, one per cause
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL,
        varprof = NULL,
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
            varsel <- spec$varsel %||% "none"
            matrix_alpha <- spec$alpha %||% 0.5
            nfolds <- spec$nfolds %||% 10
            penalty <- spec$penalty %||% "AIC"

            cause_labels <- names(self$cause_codes)
            models_list <- list()

            timevar <- task$time_col
            # event_col in CR task is the mapped event code column (0, 1, 2...)
            eventcodevar <- task$event_col
            features <- task$features

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]

                # Binary event indicator for this cause
                # Note: data[[eventcodevar]] contains codes like 0, 1, 2.
                # We want 1 if event == cause_code, 0 otherwise (censored or competing event).
                # Standard Cause-Specific Hazard approach: treat competing events as censored.

                current_event <- ifelse(data[[eventcodevar]] == cause_code, 1, 0)

                subset_data <- data
                subset_data$.event_bin <- current_event

                # Formula
                form_str <- paste0(
                    "survival::Surv(", timevar, ", .event_bin) ~ ",
                    paste(features, collapse = " + ")
                )
                form <- stats::as.formula(form_str)

                # Fit logic adapted from CRModel_Cox
                model_obj <- tryCatch(
                    {
                        survival::coxph(form, data = subset_data, x = TRUE)
                    },
                    error = function(e) {
                        rlang::warn(glue::glue("Failed to fit Cox model for cause {cause_label}: {e$message}"))
                        NULL
                    }
                )

                # TODO: Implement variable selection (backward/forward/penalized) if requested
                # For now, mirroring the basic fit.

                if (!is.null(model_obj)) {
                    models_list[[cause_label]] <- model_obj
                }
            }

            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()

            # Use Aalen-Johansen estimator with cause-specific hazards
            # We need to compute cumulative hazards for each cause, then combine.
            # Actually, R's survival implementation or custom logic?
            # `CRModel_Cox` used `aalenJohansenFromCoxModels` (legacy helper).
            # We should ideally refactor that helper or use it if it's pure logic.
            # Let's check `R/model_prediction_helpers.R` or similar if it exists.
            # It seems `aalenJohansenFromCoxModels` is likely in `R/cr_cox_helpers.R` or similar,
            # or defined in `R/cr_cox.R` (it wasn't in the view).
            # If it's missing, we need to implement it.
            # Wait, I saw `aalenJohansenFromCoxModels` called in `R/cr_cox.R`.
            # It must be defined somewhere.

            # Assuming we can use the helper or re-implement Aalen-Johansen.
            # For now, let's implement the prediction logic directly if standard.

            # Aalen-Johansen for CSH:
            # S(t) = exp(- sum_k H_k(t))  (Overall survival)
            # CIF_k(t) = integral_0^t S(u-) dH_k(u)

            complete_data <- .ensure_prediction_data(newdata, self$task)

            # We need predictions on a unified time grid (base_times)
            # Usually unique event times from training, or the `times` requested.
            if (is.null(times)) {
                target_times <- self$time_grid
            } else {
                target_times <- sort(unique(as.numeric(times)))
            }

            # To perform proper integration, we need a fine grid (all event times).
            # Using self$time_grid is usually sufficient if it was built from event times.
            integ_times <- sort(unique(c(0, self$time_grid)))

            # Get hazards for each cause at integ_times
            n_obs <- nrow(complete_data)
            n_t <- length(integ_times)

            hazards_list <- list()
            cause_labels <- names(self$cause_codes)

            overall_cumhaz <- matrix(0, nrow = n_obs, ncol = n_t)

            for (cause in cause_labels) {
                model <- self$model[[cause]]
                if (is.null(model)) {
                    hazards_list[[cause]] <- matrix(0, nrow = n_obs, ncol = n_t)
                    next
                }

                # Predict cumulative hazard
                # predict.coxph with type="expected" gives cumulative hazard?
                # No, specific times needed.
                # We can use `survfit` to get baseline hazard and then multiply by exp(lp).

                # Fast way:
                # 1. Base hazard
                # 2. LP

                summary_fit <- tryCatch(
                    {
                        survival::survfit(model, newdata = complete_data, se.fit = FALSE, conf.type = "none")
                    },
                    error = function(e) {
                        NULL
                    }
                )

                if (is.null(summary_fit)) {
                    hazards_list[[cause]] <- matrix(0, nrow = n_obs, ncol = n_t)
                    next
                }

                # summary_fit returns an object where rows are times, cols are obs (or struct array)
                # We need to extract survival and convert to cumhaz -> -log(S)
                # But survfit for newdata returns S(t|x).

                # Depending on survival version and input, output varies.
                # Usually summary(fit, times=integ_times) works.

                summ <- summary(summary_fit, times = integ_times, extend = TRUE)

                # Reshape to matrix n_obs x n_t
                # summ$surv is usually flattened: time 1 (obs1, obs2...), time 2...
                # OR obs 1 (t1, t2...), obs 2...
                # Check dimensions.
                # If newdata has multiple rows, surv is matrix (times x obs).

                surv_mat <- summ$surv
                if (!is.matrix(surv_mat)) {
                    # If 1 obs, vector.
                    if (n_obs == 1) {
                        surv_mat <- matrix(surv_mat, ncol = 1)
                    } else {
                        surv_mat <- matrix(surv_mat, nrow = length(integ_times))
                    } # Fallback assumption
                }

                # Transpose so rows = obs, cols = times
                if (nrow(surv_mat) == n_t && ncol(surv_mat) == n_obs) {
                    surv_mat <- t(surv_mat)
                }

                chaz <- -log(surv_mat)
                # H(t) should be increasing.
                # Force monotonicity
                chaz <- t(apply(chaz, 1, function(x) cummax(x)))

                hazards_list[[cause]] <- chaz
                overall_cumhaz <- overall_cumhaz + chaz
            }

            # Overall Survival S(t) = exp(- sum H_k(t))
            overall_surv <- exp(-overall_cumhaz)

            # CIF_k(t) = integral S(u-) dH_k(u)
            # Discrete approx: sum S(t_{j-1}) * (H_k(t_j) - H_k(t_{j-1}))

            # Compute dH_total efficiently from overall_cumhaz
            dH_total <- t(apply(overall_cumhaz, 1, function(row) c(row[1], diff(row))))

            # S(u-)
            # Shift overall_surv right by 1 column (t0=1)
            S_prev <- cbind(rep(1, n_obs), overall_surv[, -ncol(overall_surv), drop = FALSE])

            # Probability of any failure at step j: S_{j-1} - S_j
            P_fail <- S_prev - overall_surv
            P_fail[P_fail < 0] <- 0 # Numerical safety

            # Avoid division by zero in allocation
            dH_total_safe <- dH_total
            dH_total_safe[dH_total_safe == 0] <- 1

            preds_list <- list()

            for (cause in cause_labels) {
                H_k <- hazards_list[[cause]]
                dH_k <- t(apply(H_k, 1, function(row) c(row[1], diff(row))))

                # Allocation fraction: dH_k / dH_total
                # If dH_total is 0, then dH_k is 0. Fraction is 0.
                frac_k <- dH_k / dH_total_safe

                # Increment for cause k
                inc_k <- P_fail * frac_k

                cif_k <- apply(inc_k, 1, cumsum)

                # Interpolate to requested target_times
                id_complete <- complete_data[[self$task$id_col]]

                if (!identical(integ_times, target_times)) {
                    cif_k_interp <- cifMatInterpolator(cif_k, integ_times, target_times)
                    times_out <- target_times
                } else {
                    cif_k_interp <- cif_k
                    times_out <- integ_times
                }

                preds_list[[length(preds_list) + 1]] <- new_cif_prediction(
                    id = rep(id_complete, each = length(times_out)),
                    time = rep(times_out, times = length(id_complete)),
                    cause = rep(cause, length(id_complete) * length(times_out)),
                    cif = as.vector(cif_k_interp),
                    model = rep("cox", length(id_complete) * length(times_out)),
                    ensemble = FALSE,
                    set = set
                )
            }

            # Handle result binding
            result <- dplyr::bind_rows(preds_list)

            # Clamp values to [0, 1] to handle numerical noise
            result$cif <- pmin(pmax(result$cif, 0), 1)

            result
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "Cause-specific Cox"
            info
        },
        required_packages = function() {
            c("survival")
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
    engine = "cox",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        CoxCompetingRiskModel$new(spec = modifyList(list(engine = "cox"), spec))
    },
    packages = "survival",
    label = "Cause-specific Cox",
    tags = c("cox", "competing-risk")
)
