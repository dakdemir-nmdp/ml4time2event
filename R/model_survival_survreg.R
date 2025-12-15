#' Parametric Survival Model (R6)
#'
#' Fits a parametric survival model (survreg) using `survival` package.
#' Currently supports 'exponential' distribution.
#'
#' @keywords internal
#' @noRd
SurvRegSurvivalModel <- R6::R6Class(
  classname = "SurvRegSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    feature_names = NULL,
    varprof = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      self$varprof <- VariableProfile(data, task$features)

      spec <- self$spec
      dist <- spec$dist %||% "exponential"

      # Prepare Formula
      # We fit on all features. (Removing legacy forward selection)
      # Handle potential failure if features > N or collinearity?
      # Standard survreg behavior.

      timevar <- task$time_col
      eventvar <- task$event_col
      expvars <- task$features

      # Ensure event is numeric
      data[[eventvar]] <- as.numeric(data[[eventvar]] == 1)

      formula_str <- paste("survival::Surv(", timevar, ",", eventvar, ") ~", paste(expvars, collapse = "+"))
      formula <- stats::as.formula(formula_str)

      fit_obj <- tryCatch(
        {
          survival::survreg(formula, data = data, dist = dist, x = TRUE, y = TRUE)
        },
        error = function(e) {
          rlang::abort(glue::glue("survreg fit failed: {e$message}"))
        }
      )

      self$model <- fit_obj

      # Store feature names (from model matrix or terms)
      self$feature_names <- expvars

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()

      # Prepare newdata
      complete_data <- .ensure_prediction_data(newdata, self$task)

      # Predict LP
      # predict.survreg(obj, type="lp") returns X * beta + intercept + scale?
      # type="lp" is "linear predictor".
      # For survreg, lp usually refers to the location parameter mu?
      # mu = X*beta.
      # Scale is sigma.
      # dist="exponential": scale=1 (fixed).
      # T ~ Exp(lambda). S(t) = exp(-lambda t).
      # survreg parameterizes usually as log(T) = X beta + sigma W.
      # For exponential: W ~ ExtremeValue. sigma=1.
      # log(T) = X beta + W.
      # lambda = exp(-X beta).
      # S(t) = exp( - exp(-X beta) * t ).
      # Let's verify legacy code:
      # rate = exp(-lp). surv = exp(-rate * t).
      # So lp = X beta.
      # rate = exp(-lp) -> lambda = exp(-X beta).
      # Matches log(T) = -log(lambda) + W?
      # If mean log time matches lp.

      lp <- predict(self$model, newdata = complete_data, type = "lp")

      # Dist check
      dist <- self$model$dist
      if (dist != "exponential") {
        # Legacy only supported exponential
        rlang::warn("Only 'exponential' distribution is currently supported for prediction.")
        # Return NA?
      }

      n_obs <- nrow(complete_data)
      req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Calc probabilities
      # S(t) = exp( - exp(-lp) * t )

      rate <- exp(-lp)
      # Matrix: rows=obs, cols=times (outer)
      # outer(rate, req_times) -> rate[i] * t[j]
      hazard_cum <- outer(rate, req_times)
      surv_probs <- exp(-hazard_cum)

      # Flatten
      id_col <- complete_data[[self$task$id_col]]

      new_survival_prediction(
        id = rep(id_col, each = length(req_times)),
        time = rep(req_times, times = length(id_col)),
        surv = as.vector(t(surv_probs)), # Transpose because we want time cycling fast?
        # outer: [obs, times]. as.vector -> obs1_t1, obs2_t1... NO.
        # outer fills col-major. cols are times.
        # So [obs1_t1, obs1_t2... (col1), obs2_t1...(col1)...] -> NO.
        # Matrix is [obs, times].
        # as.vector flattens column by column.
        # So t1 for all obs, then t2 for all obs.
        # new_survival_prediction expects: id1_t1, id1_t2, ... id2_t1...
        # So we need to transpose to [times, obs] then flatten?
        # Or transpose to [obs, times] but we want row-major?
        # R is column major.
        # If we want id1_t1, id1_t2... we equivalent to row-major of [obs, times].
        # So transpose [obs, times] -> [times, obs]. Flatten [times, obs] -> t1_id1, t2_id1... Correct.

        model = rep("survreg", length(id_col) * length(req_times)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "Parametric Survival (Exponential)"
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
  engine = "survreg",
  outcome = "survival",
  constructor = function(spec = list()) {
    SurvRegSurvivalModel$new(spec = modifyList(list(engine = "survreg"), spec))
  },
  packages = "survival",
  tags = c("parametric"),
  label = "Parametric survival regression"
)
