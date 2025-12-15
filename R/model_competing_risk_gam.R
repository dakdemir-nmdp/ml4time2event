#' Competing Risk Generalized Additive Model (R6)
#'
#' Fits a cause-specific GAM model for each competing risk using `mgcv`.
#'
#' @keywords internal
#' @noRd
GamCompetingRiskModel <- R6::R6Class(
    classname = "GamCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of models
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL,
        varprof = NULL,
        # Need to store training time grids per cause for efficient Aalen-Johansen?
        # Or just store the fitted models and use predict.

        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            self$varprof <- VariableProfile(data, task$features)
            cause_map <- task$metadata$cause_map
            self$cause_codes <- setNames(cause_map$code, as.character(cause_map$cause))

            spec <- self$spec
            shrinkTreshold <- spec$shrinkTreshold %||% 10

            timevar <- task$time_col
            eventcodevar <- task$event_col
            expvars <- task$features

            # Determine variable types for formula construction (same as Survival GAM)
            df_vars <- data[, expvars, drop = FALSE]
            is_num <- sapply(df_vars, is.numeric)
            is_fac <- sapply(df_vars, function(x) is.factor(x) || is.character(x))
            fctvars <- expvars[is_fac]
            numvars <- expvars[is_num]

            # Factor conversion
            for (v in fctvars) if (is.character(data[[v]])) data[[v]] <- as.factor(data[[v]])

            # Formula helpers
            cat_shrink <- character(0)
            cat_noshrink <- character(0)
            if (length(fctvars) > 0) {
                lvls <- sapply(data[, fctvars, drop = FALSE], function(x) length(levels(as.factor(x))))
                cat_shrink <- fctvars[lvls > shrinkTreshold]
                cat_noshrink <- fctvars[lvls <= shrinkTreshold]
            }
            num_smooth <- character(0)
            num_linear <- character(0)
            if (length(numvars) > 0) {
                uniques <- sapply(data[, numvars, drop = FALSE], function(x) length(unique(x[!is.na(x)])))
                num_smooth <- numvars[uniques > shrinkTreshold]
                num_linear <- numvars[uniques <= shrinkTreshold]
            }

            terms <- c()
            if (length(cat_shrink) > 0) terms <- c(terms, paste0("s(", cat_shrink, ", bs='re')"))
            if (length(num_smooth) > 0) terms <- c(terms, paste0("s(", num_smooth, ")"))
            terms <- c(terms, num_linear, cat_noshrink)
            rhs <- if (length(terms) > 0) paste(terms, collapse = "+") else "1"

            # Models per cause
            models_list <- list()
            cause_labels <- names(self$cause_codes)

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]

                # Cause-specific status
                status_cs <- ifelse(data[[eventcodevar]] == cause_code, 1L, 0L)
                if (sum(status_cs) < 5) {
                    rlang::warn(glue::glue("Fewer than 5 events used for cause {cause_label} GAM fit."))
                }
                if (sum(status_cs) == 0) {
                    rlang::warn(glue::glue("No events for cause {cause_label}. Skipping."))
                    next
                }

                # Construct data with .status_cs
                data$.status_cs <- status_cs

                formula <- stats::as.formula(paste0("survival::Surv(", timevar, ", .status_cs) ~ ", rhs))

                gam_fit <- tryCatch(
                    {
                        mgcv::gam(
                            formula,
                            family = mgcv::cox.ph(),
                            data = data,
                            weights = .status_cs, # weights
                            select = TRUE
                        )
                    },
                    error = function(e) {
                        rlang::warn(glue::glue("GAM fit failed for cause {cause_label}: {e$message}"))
                        NULL
                    }
                )

                if (!is.null(gam_fit)) {
                    # Estimate baseline hazard
                    lp <- stats::predict(gam_fit, newdata = data, type = "link")
                    # Use the score2proba trick (Cox model with fixed coef)
                    base_cox <- survival::coxph(
                        survival::Surv(time, status) ~ score,
                        data = data.frame(time = data[[timevar]], status = status_cs, score = lp),
                        init = 1,
                        control = survival::coxph.control(iter.max = 0)
                    )
                    gam_fit$baseline_model <- base_cox
                    models_list[[cause_label]] <- gam_fit
                }
            }

            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()
            complete_data <- .ensure_prediction_data(newdata, self$task)
            n_obs <- nrow(complete_data)
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))
            n_times <- length(req_times)

            cause_labels <- names(self$cause_codes)

            # Calculate Cumulative Hazards for each cause
            # Lambda_k(t, x) = H0_k(t) * exp(lp_k)

            cum_hazards_array <- array(0, dim = c(n_times, n_obs, length(cause_labels)))
            dimnames(cum_hazards_array)[[3]] <- cause_labels

            for (cause in cause_labels) {
                model <- self$model[[cause]]
                if (is.null(model)) next

                # LP
                lp <- stats::predict(model, newdata = complete_data, type = "link")

                # Baseline H0(t)
                # survfit(base_cox, newdata=0) -> S0(t). H0(t) = -log(S0(t))
                sf <- survival::survfit(model$baseline_model, newdata = data.frame(score = 0))

                # Interpolate H0 to req_times
                if (length(sf$time) > 0) {
                    # H0 at event times
                    H0_raw <- -log(pmax(sf$surv, 1e-10))
                    # Interpolate
                    H0_interp <- stats::approx(
                        x = c(0, sf$time),
                        y = c(0, H0_raw),
                        xout = req_times,
                        method = "constant",
                        rule = 2
                    )$y

                    # Compute Lambda_k
                    # [times] * [obs] -> [times, obs]
                    haz_mat <- outer(H0_interp, exp(lp), "*")
                    cum_hazards_array[, , cause] <- haz_mat
                }
            }

            # Total Cum Hazard and Overall Survival
            overall_cum_haz <- apply(cum_hazards_array, c(1, 2), sum) # [n_times, n_obs]
            S_overall <- exp(-overall_cum_haz)

            # Compute CIF using Aalen-Johansen
            # CIF_k(t) = Integral S(u-) dLambda_k(u)
            # Discrete approx: CIF_k(t_i) = CIF_k(t_{i-1}) + S(t_{i-1}) * (Lambda_k(t_i) - Lambda_k(t_{i-1}))

            preds_list <- list()
            id_complete <- complete_data[[self$task$id_col]]

            for (cause in cause_labels) {
                cif_mat <- matrix(0, nrow = n_times, ncol = n_obs)
                cause_haz <- cum_hazards_array[, , cause]

                # Loop over time steps
                # t=1
                cif_mat[1, ] <- S_overall[1, ] * cause_haz[1, ] # S(t_0)=1 approx?
                # Better: use S(t_{i-1}) ~ 1 for first step if t0=0.
                # Usually AJ: sum S(tj-1) * dHaz(t).
                # dHaz(1) = Haz(1) - Haz(0) = Haz(1).
                # S(0) = 1.
                cif_mat[1, ] <- 1.0 * cause_haz[1, ]

                if (n_times > 1) {
                    for (t in 2:n_times) {
                        dHaz <- cause_haz[t, ] - cause_haz[t - 1, ]
                        dHaz <- pmax(dHaz, 0) # monotone
                        cif_mat[t, ] <- cif_mat[t - 1, ] + S_overall[t - 1, ] * dHaz
                    }
                }

                # Flatten [times, obs] -> [times * obs] but we need [id1, ...]
                # Transpose to [times, obs] is already current shape (rows=times).
                # But earlier I said [times, obs] -> as.vector gives t1_obs1... NO.
                # as.vector gives column-wise: t1_obs1, t2_obs1...
                # BUT R matrices [obs, times] -> [obs1_t1, obs2_t1...] is usually what we handle or vice versa.
                # Let's check new_cif_prediction expectation.
                # id = rep(id, each=n_times) -> id1, id1, id1...
                # So we want values: val(id1, t1), val(id1, t2)...
                # Matrix is [times, obs].
                # Col 1 is obs1: [t1, t2...].
                # Col 2 is obs2: [t1, t2...].
                # as.vector(Matrix) -> col1 (obs1 all times), then col2 (obs2 all times).
                # This matches id1, id1... id2, id2...
                # YES. So we just vectorize cif_mat.

                cif_vec <- as.vector(cif_mat)

                new_cif_prediction(
                    id = rep(id_complete, each = n_times),
                    time = rep(req_times, times = length(id_complete)),
                    cause = rep(cause, length(id_complete) * n_times),
                    cif = cif_vec,
                    model = rep("gam", length(id_complete) * n_times),
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
            info$label <- "GAM Competing Risk"
            info
        },
        required_packages = function() {
            c("mgcv", "survival")
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
    engine = "cr_gam",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        GamCompetingRiskModel$new(spec = modifyList(list(engine = "cr_gam"), spec))
    },
    packages = "mgcv",
    tags = c("smooth", "semiparametric"),
    label = "GAM Competing Risk"
)
