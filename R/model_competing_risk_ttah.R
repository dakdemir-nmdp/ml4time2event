#' TTAH Competing Risk Model (R6)
#'
#' Fits a discrete-time hazard model for competing risks using multinomial logistic regression
#' with time-varying coefficients and latent feature interactions.
#'
#' @keywords internal
#' @noRd
TtahCompetingRiskModel <- R6::R6Class(
    classname = "TtahCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL,
        time_grid = NULL, # The internal grid used for discretization
        task = NULL,
        varprof = NULL,
        factor_levels = NULL,
        basis_specs = NULL,
        latent_projection = NULL,
        cause_codes = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            # We ignore external time_grid for model building structure (we build our own TTAH grid)
            # or we treat it as suggestion.

            data <- as.data.frame(task$data)
            expvars <- task$features
            timevar <- task$time_col
            # eventvar in task has 0, 1, 2 codes.
            # ml4t2e_task_cr stores 'cause' column?
            # task$event_col is the numeric code column
            eventvar <- task$event_col

            spec <- self$spec
            n_time <- spec$n_time %||% 50
            spline_knots <- spec$spline_knots %||% 5
            latent_dim <- spec$latent_dim %||% 6
            time_basis_df <- spec$time_basis_df %||% 4
            lambda <- spec$lambda %||% 1e-3
            maxit <- spec$maxit %||% 200
            verbose <- FALSE

            # 1. Profile
            self$varprof <- VariableProfile(data, expvars)

            self$factor_levels <- list()
            for (v in expvars) {
                if (is.factor(data[[v]]) || is.character(data[[v]])) {
                    self$factor_levels[[v]] <- levels(as.factor(data[[v]]))
                }
            }

            obs_times <- data[[timevar]]
            obs_events <- data[[eventvar]] # 0, 1, 2...

            available_events <- sort(unique(obs_events[obs_events != 0]))
            if (length(available_events) == 0) rlang::abort("No events.")

            self$cause_codes <- as.character(available_events)

            # Grid
            tg_pass <- if (!is.null(time_grid)) time_grid[time_grid > 0] else NULL
            if (!is.null(tg_pass) && length(tg_pass) == 0) tg_pass <- NULL
            grid <- ttah_build_time_grid(obs_times, time_grid = tg_pass, n_time = n_time)
            self$time_grid <- grid

            # Features
            prep <- ttah_prepare_features(
                data = data[, expvars, drop = FALSE],
                expvars = expvars,
                factor_levels = self$factor_levels,
                spline_knots = spline_knots
            )
            phi <- prep$phi
            self$basis_specs <- prep$basis_specs
            self$factor_levels <- prep$factor_levels

            # Latent
            self$latent_projection <- ttah_compute_latent_projection(phi, latent_dim = latent_dim)
            phi_latent <- if (ncol(self$latent_projection) > 0) {
                tmp <- phi %*% self$latent_projection
                colnames(tmp) <- colnames(self$latent_projection)
                tmp
            } else {
                NULL
            }

            # Time Basis
            time_basis <- ttah_time_basis(grid, df = time_basis_df)
            time_basis_matrix <- time_basis$matrix
            if (ncol(time_basis_matrix) > 0) {
                colnames(time_basis_matrix) <- paste0("time_basis", seq_len(ncol(time_basis_matrix)))
            }

            # Multiclass Design
            interval_index <- ttah_assign_intervals(obs_times, grid)

            multiclass_design <- ttah_build_multiclass_design(
                phi = phi,
                phi_latent = phi_latent,
                time_basis = time_basis_matrix,
                interval_index = interval_index,
                event = obs_events,
                cause_codes = self$cause_codes
            )

            if (length(unique(multiclass_design$target)) < 2) {
                rlang::abort("Insufficient event diversity for multinomial fit.")
            }

            # Fit Multinom
            design_df <- as.data.frame(multiclass_design$X)
            design_df$target <- multiclass_design$target

            fit_res <- nnet::multinom(
                target ~ .,
                data = design_df,
                decay = lambda,
                maxit = maxit,
                trace = verbose,
                MaxNWts = 20000
            )

            self$model <- list(
                multinom = fit_res,
                feature_columns = setdiff(colnames(design_df), "target"),
                time_basis_specs = time_basis$specs,
                time_basis_names = colnames(time_basis_matrix),
                cause_levels = multiclass_design$levels, # Includes no_event
                lambda = lambda
            )

            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()
            complete_data <- .ensure_prediction_data(newdata, self$task)
            expvars <- self$task$features

            # 1. Alignment
            for (v in expvars) {
                if (v %in% names(self$factor_levels)) {
                    complete_data[[v]] <- factor(complete_data[[v]], levels = self$factor_levels[[v]])
                }
            }

            # 2. Features
            prep_new <- ttah_prepare_features(
                data = complete_data,
                expvars = expvars,
                factor_levels = self$factor_levels,
                basis_specs = self$basis_specs
            )
            phi_new <- prep_new$phi

            latent_dim <- if (is.null(self$latent_projection)) 0 else ncol(self$latent_projection)
            phi_latent_new <- if (latent_dim > 0) {
                phi_new %*% self$latent_projection
            } else {
                NULL
            }

            # 3. Time Basis
            time_grid <- self$time_grid
            time_basis <- ttah_eval_time_basis(time_grid, self$model$time_basis_specs)
            if (!is.null(self$model$time_basis_names)) colnames(time_basis) <- self$model$time_basis_names

            n_obs <- nrow(phi_new)
            K <- length(time_grid)

            # 4. Long Design for Prediction (expand each obs into K rows)
            # This is expensive but matches multinom predict structure.
            # Alternatively, compute linear predictor manually.
            # Multinomial coeffs: (K_classes - 1) sets of coeffs.
            # Coeffs are relative to baseline? nnet::multinom baseline is first level.
            # model$cause_levels: "no_event", "cause_1"...
            # multinom output for K classes.

            # It's safer to use predict(multinom) on the design matrix.
            pred_desc <- ttah_build_long_design(
                phi = phi_new,
                phi_latent = phi_latent_new,
                time_basis = time_basis,
                interval_index = rep(K, n_obs), # expand fully
                event = rep(0, n_obs)
            )

            design_df <- as.data.frame(pred_desc$X)
            # Ensure cols match
            missing <- setdiff(self$model$feature_columns, colnames(design_df))
            if (length(missing) > 0) {
                zeros <- matrix(0, nrow = nrow(design_df), ncol = length(missing))
                colnames(zeros) <- missing
                design_df <- cbind(design_df, zeros)
            }
            design_df <- design_df[, self$model$feature_columns, drop = FALSE]

            probs <- predict(self$model$multinom, newdata = design_df, type = "probs")

            # Format probs [TotalRows x Classes]
            # Classes: no_event, cause_X...
            # Align:
            prob_cols <- setdiff(self$model$cause_levels, "no_event")
            if (is.null(dim(probs))) { # Single row or single class?
                probs <- matrix(probs, nrow = nrow(design_df))
                # If strictly 2 classes (no_event vs 1 cause), probs is vector of prob(Event).
                if (length(self$model$cause_levels) == 2) {
                    # vector implies P(Event).
                    probs <- cbind(1 - probs, probs)
                    colnames(probs) <- self$model$cause_levels
                } else {
                    # Should have dims
                }
            }
            # Assuming typical multinom matrix output

            # Map to (Time, Obs, Cause)
            # design_df rows correspond to: Obs 1 (t=1..K), Obs 2 (t=1..K)... fit order:
            # ttah_build_long_design iterates i=1..n, k=1..K.
            # So blocks of K rows per observation.

            # We need Hazard[k, i, c] = P(Event_c at k | Survival to k)
            # This is exactly what the multinom on (Person-Time) predicts.

            hazards_array <- array(0, dim = c(K, n_obs, length(self$cause_codes)))
            stay_matrix <- matrix(0, nrow = K, ncol = n_obs)

            cause_indices <- match(self$cause_codes, self$model$cause_levels)
            # But wait, predict gives columns by name.
            # Let's use column names if available.
            # If probs is matrix, colnames match levels.
            p_names <- colnames(probs)

            # We assume `pred_desc$obs` and `pred_desc$time` map correctly.
            obs_vec <- pred_desc$obs
            time_vec <- pred_desc$time

            # Vectorized filling?
            # hazards_array[time_vec, obs_vec, c] <- probs[, col_c]
            # Indices are flattened.

            for (idx in seq_along(self$cause_codes)) {
                c_code <- self$cause_codes[idx]
                c_label <- paste0("cause_", c_code)
                if (c_label %in% p_names) {
                    h_vals <- probs[, c_label]
                    # Map to array
                    # Matrix indexing: hazards_array[cbind(time, obs, idx)]
                    hazards_array[cbind(time_vec, obs_vec, rep(idx, length(time_vec)))] <- h_vals
                }
            }

            if ("no_event" %in% p_names) {
                stay_matrix[cbind(time_vec, obs_vec)] <- probs[, "no_event"]
            } else {
                # recalculate from 1 - sum(hazards)
                # stay_matrix = 1 - sum over causes
                sum_haz <- apply(hazards_array, c(1, 2), sum)
                stay_matrix <- 1 - sum_haz
            }

            # 5. Compute CIF
            # CIF_c(k+1) = CIF_c(k) + S_total(k) * hazard_c(k)
            # S_total(k+1) = S_total(k) * stay(k)

            Times <- c(0, time_grid)
            K_plus <- K + 1

            CIF_array <- array(0, dim = c(K_plus, n_obs, length(self$cause_codes)))
            S_total <- matrix(1, nrow = K_plus, ncol = n_obs)

            for (k in 1:K) {
                S_prev <- S_total[k, ]
                # Update CIFs
                for (idx in seq_along(self$cause_codes)) {
                    CIF_array[k + 1, , idx] <- CIF_array[k, , idx] + S_prev * hazards_array[k, , idx]
                }
                # Update S_total
                S_total[k + 1, ] <- S_prev * stay_matrix[k, ]
            }

            # 6. Interpolate to req_times
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))
            preds_list <- list()

            id_col <- complete_data[[self$task$id_col]]

            for (idx in seq_along(self$cause_codes)) {
                # Interpolate this cause's CIF matrix
                cif_mat <- CIF_array[, , idx]
                cif_interp <- cifMatInterpolator(cif_mat, Times, req_times)

                flat_cif <- as.vector(cif_interp)

                preds_list[[length(preds_list) + 1]] <- new_cif_prediction(
                    id = rep(id_col, each = length(req_times)),
                    time = rep(req_times, times = length(id_col)),
                    cause = rep(self$cause_codes[idx], length(flat_cif)),
                    cif = flat_cif,
                    model = rep("ttah", length(flat_cif)),
                    ensemble = FALSE,
                    set = set
                )
            }

            dplyr::bind_rows(preds_list)
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "TTAH Competing Risk"
            info
        },
        required_packages = function() {
            c("nnet")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
        }
    )
)

.register_time_to_event_model(
    engine = "cr_ttah",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        TtahCompetingRiskModel$new(spec = modifyList(list(engine = "cr_ttah"), spec))
    },
    packages = c("nnet"),
    tags = c("semiparametric", "discrete-time"),
    label = "TTAH Competing Risk"
)
