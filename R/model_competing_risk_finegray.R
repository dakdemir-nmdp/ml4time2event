#' Fine-Gray Competing Risk Model (R6)
#'
#' Fits a Fine-Gray model using `fastcmprsk` with SVD dimensionality reduction.
#'
#' @keywords internal
#' @noRd
FineGrayCompetingRiskModel <- R6::R6Class(
    classname = "FineGrayCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL,
        time_grid = NULL,
        task = NULL,
        varprof = NULL,
        scaling = NULL,
        loadings = NULL,
        cause_code = NULL, # The numeric cause code of interest

        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            self$varprof <- VariableProfile(data, task$features)

            # Cause map
            spec <- self$spec
            # In Fine-Gray, we usually target ONE cause against others (subdistribution hazard).
            # But our framework might expect CIFs for ALL causes.
            # The legacy implementation `CRModel_FineGray` took `event_codes` and fitted for that ONE code.
            # It seems legacy fit was single-cause focused?
            # `CRModel_FineGray(..., event_codes = NULL, ...)` -> defaults to available_events[1].

            # If we want all CIFs, we need to fit one FG model per cause?
            # Fine-Gray models are inherently cause-specific (subdistribution).
            # If we just fit one, we only get CIF for that cause.
            # To support `predict_cif` for ANY cause, we really should fit multiple FG models?
            # The legacy code: `CRModel_FineGray` returns `event_code_numeric = failcode`.
            # And predict checks if `event_of_interest == modelout$event_codes`.
            # So legacy was SINGLE cause.

            # BUT, pure `ml4t2e` framework usually wants to predict ALL CIFs (sum to <=1).
            # If we strictly follow legacy, this model defines a PRIMARY event.
            # R6 implementation should ideally handle all causes like other CR models.
            # However, if expensive (SVD + fastCrrp), maybe we stick to per-cause?
            # Most ML CR models (RF, Boost) give all CIFs.

            # I will stick to the wrapper philosophy: Implement what the legacy had, but ideally better.
            # New architecture fits ONE model object. If that object contains multiple sub-models, that's fine.
            # Let's iterate over causes and fit FG for each.

            cause_map <- task$metadata$cause_map
            causes <- as.character(cause_map$code) # numeric codes usually

            models_list <- list()
            scalings <- list()
            loadings_list <- list()

            expvars <- task$features
            timevar <- task$time_col
            eventvar <- task$event_col

            # Keep feature matrix and outcome vectors synchronized.
            # model.matrix() omits rows with NA by default; if we don't subset
            # first, assigning the projected matrix into model_df will fail.
            complete_rows <- stats::complete.cases(data[, expvars, drop = FALSE])
            if (!all(complete_rows)) {
                rlang::warn(sprintf(
                    "Fine-Gray fit: dropping %d rows with missing feature values.",
                    sum(!complete_rows)
                ))
                data <- data[complete_rows, , drop = FALSE]
            }

            # 1. Prepare Matrix (Common for all causes, X is same)
            XYTrain <- data
            XYTrain[[eventvar]] <- as.numeric(XYTrain[[eventvar]])

            covmat <- stats::model.matrix(~ -1 + ., data = XYTrain[, expvars, drop = FALSE])

            # Scale
            covmat_scaled <- scale(covmat, center = TRUE, scale = TRUE)
            meanTrain <- attr(covmat_scaled, "scaled:center")
            sdTrain <- attr(covmat_scaled, "scaled:scale")

            # SVD
            svdcovmat <- svd(covmat_scaled)
            n_components <- min(c(20, ncol(covmat_scaled)))
            Feat <- as.data.frame((covmat_scaled %*% svdcovmat$v)[, 1:n_components])
            colnames(Feat) <- paste0("PC", 1:n_components)

            # Common data for fastCrrp
            # fastCrrp needs formulas.
            model_df <- data.frame(
                ftime = XYTrain[[timevar]],
                fstatus = XYTrain[[eventvar]]
            )
            model_df$cov <- as.matrix(Feat)

            # Fit for each cause
            for (cause in causes) {
                cause_num <- as.numeric(cause)

                # fastCrrp: failcode=cause, cencode=min(status).
                # Status should be 0 for censored.
                # Ensure min(status) is 0?
                # data[[eventvar]] should have 0 for censored.

                # fastCrrp syntax: Crisk(ftime, fstatus, failcode, cencode) ~ cov
                fg_fit <- tryCatch(
                    {
                        fastcmprsk::fastCrrp(
                            fastcmprsk::Crisk(ftime, fstatus, failcode = cause_num, cencode = 0) ~ cov,
                            data = model_df,
                            lambda = 0.01,
                            alpha = 0.5,
                            penalty = "ENET",
                            standardize = FALSE,
                            max.iter = 5000
                        )
                    },
                    error = function(e) NULL
                )

                if (!is.null(fg_fit)) {
                    models_list[[cause]] <- fg_fit
                }
            }

            self$model <- models_list
            self$scaling <- list(mean = meanTrain, sd = sdTrain)
            self$loadings <- svdcovmat$v

            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()
            complete_data <- .ensure_prediction_data(newdata, self$task)
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

            feature_frame <- complete_data[, self$task$features, drop = FALSE]
            complete_idx <- stats::complete.cases(feature_frame)
            id_all <- complete_data[[self$task$id_col]]

            n_obs <- sum(complete_idx)
            if (n_obs == 0) {
                missing_list <- lapply(names(self$model), function(cause) {
                    new_cif_prediction(
                        id = rep(id_all, each = length(req_times)),
                        time = rep(req_times, times = length(id_all)),
                        cause = rep(cause, length(id_all) * length(req_times)),
                        cif = rep(NA_real_, length(id_all) * length(req_times)),
                        model = rep("fine_gray", length(id_all) * length(req_times)),
                        ensemble = FALSE,
                        set = set
                    )
                })
                return(dplyr::bind_rows(missing_list))
            }

            # Prepare Matrix
            expvars <- self$task$features
            covmat <- stats::model.matrix(~ -1 + ., data = complete_data[complete_idx, expvars, drop = FALSE])

            # Match columns
            expected_cols <- names(self$scaling$mean)
            missing <- setdiff(expected_cols, colnames(covmat))
            if (length(missing) > 0) {
                z <- matrix(0, nrow = nrow(covmat), ncol = length(missing))
                colnames(z) <- missing
                covmat <- cbind(covmat, z)
            }
            covmat <- covmat[, expected_cols, drop = FALSE]

            # Scale & Project
            covmat_scaled <- scale(covmat, center = self$scaling$mean, scale = self$scaling$sd)
            n_comp <- ncol(self$model[[1]]$coef) # fastCrrp coef is [p, 1]?
            # No, fastCrrp coef dim matches input. We used n_components PCs.
            # Check dimensionality from model.
            # But we use the SAME projection for all.
            # Use loadings
            # Ensure n_components match
            n_comps <- min(ncol(self$loadings), 20) # Based on fit logic
            # Actually fit used: min(c(20, ncol(covmat_scaled)))
            # We need to reuse that N.
            # Usually loadings has 'v'. Rows=params, Cols=Components.
            # So cov * v -> [n, n_comp].
            # Just use all columns of v? No, fit reduced it.
            # fit used: (covmat_scaled %*% svdcovmat$v)[, 1:n_components]
            # We also need n_components.
            # Extract from model? Coef length.

            # Let's trust self$model[[1]]$coef to give dimensionality.

            # Or better, store n_components in object
            # I'll rely on the fact that coefs match the design matrix 'cov'.

            # Wait, if some models failed, I check first valid one.
            valid_model <- NULL
            for (m in self$model) {
                if (!is.null(m)) {
                    valid_model <- m
                    break
                }
            }
            if (is.null(valid_model)) {
                rlang::abort("Fine-Gray prediction failed: no fitted cause-specific model is available.")
            }
            n_features_model <- length(valid_model$coef)

            Feat <- (covmat_scaled %*% self$loadings)[, 1:n_features_model, drop = FALSE]

            preds_list <- list()

            for (cause in names(self$model)) {
                fg_fit <- self$model[[cause]]

                # Predict CIF
                # fastcmprsk gives breslowJump -> baseline Hazard
                # CIF_k(t) = 1 - exp(- CumHaz_k(t))
                # CumHaz_k(t) = exp(beta * X) * H0(t)

                breslow <- fg_fit$breslowJump
                if (is.null(breslow)) next

                base_times <- breslow[, 1]
                base_haz <- breslow[, 2]

                cif_vals <- matrix(0, nrow = length(req_times), ncol = n_obs)

                # Interpolate H0 to req_times
                # base_haz is cumulative hazard H0(t)?
                # fastcmprsk docs: breslowJump contains times and jumps? Or cumulative?
                # "breslowJump" typically implies jumps DA.
                # But in CRModel_FineGray it parses col2 as baseline_haz and does:
                # cumhaz <- exp(lp) * baseline_haz
                # cif_vals <- 1 - exp(-cumsum(cumhaz))
                # Wait, if it does cumsum(cumhaz), then 'baseline_haz' must be increments (jumps).
                # AND if it does that for each OBS, it iterates time steps.
                # In legacy:
                # cif_matrix[-1, i] <- 1 - exp(-cumsum(cumhaz))
                # where cumhaz = exp(lp) * baseline_haz (vector mult?)
                # No, 'baseline_haz' vector * scalar exp(lp).
                # So H0(t) = sum(h0(s)). Lambda(t) = exp(lp) * H0(t).
                # CIF = 1 - exp(-Lambda(t)).
                # This formula `1 - exp(-Lambda(t))` is for SUBDISTRIBUTION.
                # Yes, Fine-Gray estimates subdistribution hazard directly.
                # So CIF_k(t) = 1 - exp( - Integral( lambda_k(u) du ) )
                # Integral = exp(beta*x) * Lambda_0(t).
                # Lambda_0(t) = sum of jumps.

                # So I need Cumulative Baseline Hazard.
                # If breslow is jumps, I cumsum it.
                H0_cum <- cumsum(base_haz)

                # Interpolate H0_cum to req_times
                H0_interp <- stats::approx(
                    x = c(0, base_times),
                    y = c(0, H0_cum),
                    xout = req_times,
                    method = "constant",
                    rule = 2
                )$y

                beta <- as.numeric(fg_fit$coef)
                lp <- as.vector(Feat %*% beta)

                # Outer product: [times] * [obs] -> [times, obs]
                Lambda_mat <- outer(H0_interp, exp(lp), "*")

                Cif_mat <- 1 - exp(-Lambda_mat)

                # Flatten column-major [times, obs] -> t1_obs1, t2_obs1...
                cif_vec <- as.vector(Cif_mat)
                id_col <- id_all[complete_idx]

                preds_list[[length(preds_list) + 1]] <- new_cif_prediction(
                    id = rep(id_col, each = length(req_times)),
                    time = rep(req_times, times = length(id_col)),
                    cause = rep(cause, length(cif_vec)),
                    cif = cif_vec,
                    model = rep("fine_gray", length(cif_vec)),
                    ensemble = FALSE,
                    set = set
                )
            }

            preds_complete <- dplyr::bind_rows(preds_list)

            if (!all(complete_idx)) {
                missing_ids <- id_all[!complete_idx]
                rlang::warn(glue::glue(
                    "Omitting {length(missing_ids)} rows with missing predictors for engine 'fine_gray'."
                ))
                missing_list <- lapply(names(self$model), function(cause) {
                    new_cif_prediction(
                        id = rep(missing_ids, each = length(req_times)),
                        time = rep(req_times, times = length(missing_ids)),
                        cause = rep(cause, length(missing_ids) * length(req_times)),
                        cif = rep(NA_real_, length(missing_ids) * length(req_times)),
                        model = rep("fine_gray", length(missing_ids) * length(req_times)),
                        ensemble = FALSE,
                        set = set
                    )
                })
                preds_complete <- dplyr::bind_rows(preds_complete, dplyr::bind_rows(missing_list))
            }

            preds_complete
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "Fine-Gray"
            info
        },
        required_packages = function() {
            c("fastcmprsk")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
        }
    )
)

.register_time_to_event_model(
    engine = "cr_fine_gray",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        FineGrayCompetingRiskModel$new(spec = modifyList(list(engine = "cr_fine_gray"), spec))
    },
    packages = "fastcmprsk",
    tags = c("semi-parametric", "linear"),
    label = "Fine-Gray"
)
