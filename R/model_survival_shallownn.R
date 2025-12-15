#' Shallow Neural Network Survival Model (R6)
#'
#' Fits a shallow neural network with a Cox partial likelihood loss.
#' Uses shared neural network core functions from nn_core.R for efficiency.
#'
#' @keywords internal
#' @noRd
ShallowNNSurvivalModel <- R6::R6Class(
  classname = "ShallowNNSurvivalModel",
  inherit = SurvivalModel,
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

      # Binary event indicator
      event <- as.numeric(data[[eventvar]] == 1)
      time <- data[[timevar]]

      # Fit using shared Cox loss optimization
      fit_result <- .nn_fit_cox_loss(
        X_train = prepared$X,
        event = event,
        time = time,
        n_features = prepared$n_features,
        size = size,
        decay = decay,
        maxit = maxit
      )

      # Check convergence
      if (fit_result$convergence != 0) {
        rlang::warn(glue::glue(
          "Neural network optimization did not converge fully (code {fit_result$convergence}). ",
          "Consider increasing maxit or adjusting decay parameter."
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

      # Store model
      self$model <- list(
        weights = fit_result$weights,
        baseline_hazard = baseline_hazard,
        scaling_params = prepared$scaling_params,
        factor_levels = prepared$factor_levels,
        size = size
      )

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)

      # Use shared prediction function
      .nn_predict_survival(
        model = self$model,
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
      info$label <- "Shallow Neural Network (Survival)"
      info$parameters <- list(
        hidden_units = self$model$size %||% "unknown",
        converged = TRUE # Would need to store this from fit
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
          "Cannot predict: ShallowNN model has not been fitted yet.",
          class = "unfitted_model_error"
        )
      }
      if (is.null(self$model$weights)) {
        rlang::abort(
          "Cannot predict: Model weights are missing.",
          class = "invalid_model_state_error"
        )
      }
    }
  )
)

# Register model
.register_time_to_event_model(
  engine = "shallownn",
  outcome = "survival",
  constructor = function(spec = list()) {
    ShallowNNSurvivalModel$new(spec = modifyList(list(engine = "shallownn"), spec))
  },
  packages = character(),
  tags = c("neural-network", "cox-loss", "nonparametric"),
  label = "Shallow Neural Network (Survival)"
)
