#' Parametric Competing Risk Model (R6)
#'
#' Fits a cause-specific Parametric model (survreg) for each competing risk.
#' Assumes Exponential distribution for closed-form CIF.
#'
#' @keywords internal
#' @noRd
SurvRegCompetingRiskModel <- R6::R6Class(
    classname = "SurvRegCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of models
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL,
        feature_names = NULL,
        varprof = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            cause_map <- task$metadata$cause_map
            self$cause_codes <- setNames(cause_map$code, as.character(cause_map$cause))
            self$varprof <- VariableProfile(data, task$features)

            spec <- self$spec
            dist <- spec$dist %||% "exponential"
            if (dist != "exponential") {
                rlang::warn("SurvRegCompetingRiskModel currently only supports 'exponential' distribution for valid CIF calculation. Defaulting to 'exponential'.")
                dist <- "exponential"
            }

            timevar <- task$time_col
            eventcodevar <- task$event_col

            models_list <- list()
            cause_labels <- names(self$cause_codes)

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]

                # Cause-Specific Hazard:
                # Event k is event.
                # All others + censoring -> 0 (censored).

                y_codes <- data[[eventcodevar]]
                status_cs <- ifelse(y_codes == cause_code, 1, 0)

                if (sum(status_cs) < 1) {
                    rlang::warn(glue::glue("No events for cause {cause_label}. Skipping."))
                    next
                }

                # Construct data for this cause
                # We can use the original dataframe but need a new status col
                # Or construct formula with I()

                # Formula: Surv(time, status_cs) ~ ...
                # But status_cs is not in data.
                data_cs <- data
                data_cs$.status_cs <- status_cs

                formula_str <- paste("survival::Surv(", timevar, ", .status_cs) ~", paste(task$features, collapse = "+"))
                formula <- stats::as.formula(formula_str)

                fit_obj <- tryCatch(
                    {
                        survival::survreg(formula, data = data_cs, dist = dist, x = TRUE, y = TRUE)
                    },
                    error = function(e) {
                        rlang::warn(glue::glue("Failed to fit model for cause {cause_label}: {e$message}"))
                        NULL
                    }
                )

                models_list[[cause_label]] <- fit_obj
            }

            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()

            complete_data <- .ensure_prediction_data(newdata, self$task)
            n_obs <- nrow(complete_data)
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

            cause_labels <- names(self$cause_codes)

            # Calculate Rates (Lambda) for each cause
            rates_list <- list()
            total_rate <- rep(0, n_obs)

            for (cause in cause_labels) {
                model <- self$model[[cause]]
                if (is.null(model)) {
                    rates_list[[cause]] <- rep(0, n_obs)
                } else {
                    # Predict LP
                    # SurvReg Exponential: rate = exp(-lp)
                    lp <- predict(model, newdata = complete_data, type = "lp")
                    rate <- exp(-lp)
                    rates_list[[cause]] <- rate
                    total_rate <- total_rate + rate
                }
            }

            preds_list <- list()
            id_complete <- complete_data[[self$task$id_col]]

            # CIF_k(t) = (lambda_k / lambda_total) * (1 - exp(-lambda_total * t))

            for (cause in cause_labels) {
                rate_k <- rates_list[[cause]]

                # Matrix [obs, times]
                # rate_k [n]
                # total_rate [n]
                # times [m]

                # ratio R_k = rate_k / total_rate (handle 0/0 -> 0)
                ratio_k <- ifelse(total_rate > 0, rate_k / total_rate, 0)

                # term T(t) = 1 - exp(-total_rate * t)
                # outer(total_rate, req_times) -> [n, m]
                exponent_mat <- outer(total_rate, req_times)
                term_mat <- 1 - exp(-exponent_mat)

                cif_mat <- ratio_k * term_mat # Broadcasting?
                # ratio_k is vector length n. term_mat is n x m.
                # R recycles column wise? NO.
                # term_mat * ratio_k (vector) -> multiplies column-wise.
                # Vector is length n. Matrix has n rows.
                # R matrices are stored column-major (m columns of length n).
                # So ratio_k aligns with columns. Perfect.

                # Flatten
                # We want [t1_id1, t2_id1...]?
                # CIF_mat rows=obs. Cols=times.
                # To get id1_t1, id1_t2... we need transpose to [times, obs] then flatten.

                cif_vec <- as.vector(t(cif_mat))

                new_cif_prediction(
                    id = rep(id_complete, each = length(req_times)),
                    time = rep(req_times, times = length(id_complete)),
                    cause = rep(cause, length(id_complete) * length(req_times)),
                    cif = cif_vec,
                    model = rep("survreg", length(id_complete) * length(req_times)),
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
            info$label <- "Parametric Competing Risk (Exponential)"
            info
        },
        required_packages = function() {
            c("survival")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) {
                rlang::abort("Model not fitted")
            }
        }
    )
)

.register_time_to_event_model(
    engine = "cr_survreg",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        SurvRegCompetingRiskModel$new(spec = modifyList(list(engine = "cr_survreg"), spec))
    },
    packages = "survival",
    tags = c("parametric"),
    label = "Parametric Competing Risk"
)
