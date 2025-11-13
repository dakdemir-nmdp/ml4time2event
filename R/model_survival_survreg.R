#' Parametric survival regression (R6 implementation)
#'
#' Leverages `SurvModel_SurvReg()` to provide a parametric engine that fits into
#' the unified `ml4time2event` interface.
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

    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      args <- c(
        list(
          data = as.data.frame(task$data),
          expvars = task$features,
          timevar = task$time_col,
          eventvar = task$event_col
        ),
        private$spec_args()
      )
      self$model <- do.call(SurvModel_SurvReg, args)
      self$time_grid <- time_grid
      self$task <- task
      invisible(self)
    },

    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)
      feature_frame <- as.data.frame(complete_data[, self$task$features, drop = FALSE])
      complete_idx <- stats::complete.cases(feature_frame)
      ids <- complete_data[[self$task$id_col]]
      model_name <- self$spec$engine %||% "survreg"

      if (is.null(times)) {
        target_times <- self$time_grid
      } else {
        target_times <- sort(unique(as.numeric(times)))
      }
      if (length(target_times) == 0) {
        rlang::abort("`times` must be numeric and non-empty.")
      }

      preds_complete <- NULL
      if (any(complete_idx)) {
        newdata_complete <- feature_frame[complete_idx, , drop = FALSE]
        id_complete <- ids[complete_idx]
        pred <- Predict_SurvModel_SurvReg(
          modelout = self$model,
          newdata = newdata_complete,
          new_times = target_times
        )
        preds_complete <- new_survival_prediction(
          id = rep(id_complete, each = length(target_times)),
          time = rep(target_times, times = length(id_complete)),
          surv = as.vector(pred$Probs),
          model = rep(model_name, length(id_complete) * length(target_times)),
          ensemble = FALSE,
          set = set
        )
      }

      preds_missing <- NULL
      if (!all(complete_idx)) {
        missing_ids <- ids[!complete_idx]
        rlang::warn(glue::glue(
          "Omitting {length(missing_ids)} rows with missing predictors for engine '{model_name}'."
        ))
        preds_missing <- new_survival_prediction(
          id = rep(missing_ids, each = length(target_times)),
          time = rep(target_times, times = length(missing_ids)),
          surv = rep(NA_real_, length(missing_ids) * length(target_times)),
          model = rep(model_name, length(missing_ids) * length(target_times)),
          ensemble = FALSE,
          set = set
        )
      }

      pieces <- list(preds_complete, preds_missing)
      pieces <- pieces[!vapply(pieces, is.null, logical(1))]
      if (length(pieces) == 0) {
        return(new_survival_prediction(
          id = integer(0),
          time = numeric(0),
          surv = numeric(0),
          model = character(0),
          ensemble = logical(0),
          set = character(0)
        ))
      }
      dplyr::bind_rows(pieces)
    },

    predict_risk = function(newdata, times = NULL, set = "test", ...) {
      private$ensure_fitted()
      if (is.null(times)) {
        times <- self$time_grid
      }
      surv_tbl <- self$predict_survival(newdata = newdata, times = times, set = set, ...)
      if (nrow(surv_tbl) == 0) {
        return(new_risk_prediction(
          id = integer(0),
          risk = numeric(0),
          model = character(0),
          time = numeric(0),
          ensemble = logical(0),
          set = character(0)
        ))
      }
      last_time <- max(times)
      last_slice <- surv_tbl[surv_tbl$time == last_time, , drop = FALSE]
      new_risk_prediction(
        id = last_slice$id,
        risk = 1 - last_slice$surv,
        model = last_slice$model,
        time = last_time,
        ensemble = last_slice$ensemble,
        set = last_slice$set
      )
    },

    model_info = function() {
      info <- super$model_info()
      info$label <- "Parametric survival regression"
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
    },

    spec_args = function() {
      spec <- self$spec %||% list()
      spec[c("engine", "package", "outcome")] <- NULL
      spec
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
