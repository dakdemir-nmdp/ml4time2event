#' Shallow Neural Network Competing Risk Model (R6)
#'
#' Fits a shallow neural network for competing risks using cause-specific
#' hazard models (Cox loss for each cause). Combines hazards via Aalen-Johansen
#' to estimate cumulative incidence functions.
#'
#' @keywords internal
#' @noRd
ShallowNNCompetingRiskModel <- R6::R6Class(
    classname = "ShallowNNCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL,
        time_grid = NULL,
        task = NULL,
        varprof = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            # Extract specifications
            spec <- self$spec
            size <- spec$size %||% 5
            decay <- spec$decay %||% 0.01
            maxit <- spec$maxit %||% 1000

            # Get data
            data <- as.data.frame(task$data)
            expvars <- task$features
            timevar <- task$time_col
            eventvar <- task$event_col

            # Variable profiling for diagnostics
            self$varprof <- VariableProfile(data, expvars)

            # Prepare data: handle factors, scaling, create model matrix
            prepared <- .nn_prepare_data(data, expvars)

            # Numeric event
            data[[eventvar]] <- as.numeric(data[[eventvar]])
            time <- data[[timevar]]

            # Train one model per cause (cause-specific hazard approach)
            cause_map <- task$metadata$cause_map
            causes <- as.character(cause_map$code)
            models_list <- list()

            for (cause_str in causes) {
                cause <- as.numeric(cause_str)

                # Binary event indicator: 1 if this cause, 0 otherwise (CSH model)
                event_binary <- as.numeric(data[[eventvar]] == cause)

                # Fit using shared Cox loss optimization
                fit_result <- .nn_fit_cox_loss(
                    X_train = prepared$X,
                    event = event_binary,
                    time = time,
                    n_features = prepared$n_features,
                    size = size,
                    decay = decay,
                    maxit = maxit
                )

                # Check convergence
                if (fit_result$convergence != 0) {
                    rlang::warn(glue::glue(
                        "Neural network optimization for cause {cause} did not converge ",
                        "fully (code {fit_result$convergence}). Consider increasing maxit."
                    ))
                }

                # Compute baseline hazard (vectorized for efficiency)
                res_final <- .nn_forward_pass(fit_result$sorted_data$X, fit_result$weights)
                risks <- exp(res_final$output)

                baseline_hazard <- .nn_compute_baseline_hazard_vectorized(
                    risks = risks,
                    time_sorted = fit_result$sorted_data$time,
                    event_sorted = fit_result$sorted_data$event
                )

                models_list[[cause_str]] <- list(
                    weights = fit_result$weights,
                    baseline_hazard = baseline_hazard,
                    scaling_params = prepared$scaling_params,
                    factor_levels = prepared$factor_levels
                )
            }

            # Store all cause-specific models
            self$model <- models_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()
            complete_data <- .ensure_prediction_data(newdata, self$task)

            # Use shared prediction function
            .nn_predict_cif(
                models_list = self$model,
                newdata = complete_data,
                times = times,
                task = self$task,
                expvars = self$task$features,
                time_grid = self$time_grid,
                set = set
            )
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "Shallow Neural Network (Competing Risks)"
            info$parameters <- list(
                hidden_units = self$model[[1]]$size %||% "unknown",
                n_causes = length(self$model)
            )
            info
        },
        required_packages = function() {
            character() # Base R only
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) {
                rlang::abort(
                    "Cannot predict: ShallowNN competing risks model has not been fitted yet.",
                    class = "unfitted_model_error"
                )
            }
            if (is.null(self$model) || length(self$model) == 0) {
                rlang::abort(
                    "Cannot predict: Model is empty or corrupted.",
                    class = "invalid_model_state_error"
                )
            }
        }
    )
)

# Register model
.register_time_to_event_model(
    engine = "cr_shallownn",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        ShallowNNCompetingRiskModel$new(spec = modifyList(list(engine = "cr_shallownn"), spec))
    },
    packages = character(),
    tags = c("neural-network", "cause-specific-hazard", "nonparametric"),
    label = "Shallow Neural Network (Competing Risks)"
)
