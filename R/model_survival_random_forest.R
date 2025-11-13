#' Random forest survival model (R6 implementation)
#'
#' Wraps `randomForestSRC::rfsrc()` to provide survival predictions on the shared
#' `TimeToEventModel` interface.
#'
#' @keywords internal
#' @noRd
RandomForestSurvivalModel <- R6::R6Class(
  classname = "RandomForestSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,

    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      data <- task$data
      required <- c(task$time_col, task$event_col, task$features)
      data <- data[, required, drop = FALSE]
      data <- stats::na.omit(data)
      data <- as.data.frame(data)

      formula <- stats::as.formula(paste0(
        "Surv(", task$time_col, ", ", task$event_col, ") ~ ",
        paste(task$features, collapse = " + ")
      ))
      spec <- self$spec
      ntree <- spec$ntree
      if (is.null(ntree)) {
        ntree <- 500L
      }
      nodesize <- spec$nodesize
      if (is.null(nodesize)) {
        nodesize <- 15
      }
      mtry <- spec$mtry
      if (is.null(mtry)) {
        mtry <- max(1, floor(sqrt(length(task$features))))
      }
      importance <- spec$importance
      if (is.null(importance)) {
        importance <- "permute"
      }
      splitrule <- spec$splitrule
      if (is.null(splitrule)) {
        splitrule <- "logrank"
      }

      self$model <- randomForestSRC::rfsrc(
        formula = formula,
        data = data,
        ntree = as.integer(ntree),
        nodesize = as.integer(nodesize),
        mtry = as.integer(mtry),
        importance = importance,
        splitrule = splitrule,
        forest = TRUE
      )
      self$time_grid <- time_grid
      self$task <- task
      invisible(self)
    },

    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      newdata <- .ensure_prediction_data(newdata, self$task)
      model_name <- self$spec$engine
      if (is.null(model_name)) {
        model_name <- "surv_rf"
      }

      complete_idx <- stats::complete.cases(newdata)
      id_values <- newdata[[self$task$id_col]]
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
        newdata_complete <- newdata[complete_idx, , drop = FALSE]
        id_complete <- newdata_complete[[self$task$id_col]]
        rf_pred <- randomForestSRC::predict.rfsrc(self$model, newdata = newdata_complete)
        base_times <- rf_pred$time.interest
        surv_matrix <- rf_pred$survival
        aligned <- .rf_align_survival(surv_matrix, base_times, target_times)

        preds_complete <- new_survival_prediction(
          id = rep(id_complete, each = length(target_times)),
          time = rep(target_times, times = length(id_complete)),
          surv = as.vector(aligned),
          model = rep(model_name, length(id_complete) * length(target_times)),
          ensemble = FALSE,
          set = set
        )
      }

      preds_missing <- NULL
      if (!all(complete_idx)) {
        missing_ids <- id_values[!complete_idx]
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
      target_times <- times
      if (is.null(target_times)) {
        target_times <- self$time_grid
      }
      survival_tbl <- self$predict_survival(newdata = newdata, times = target_times, set = set, ...)
      if (nrow(survival_tbl) == 0) {
        return(new_risk_prediction(
          id = integer(0),
          risk = numeric(0),
          model = character(0),
          time = numeric(0),
          ensemble = logical(0),
          set = character(0)
        ))
      }
      last_time <- max(target_times)
      last_slice <- survival_tbl[survival_tbl$time == last_time, , drop = FALSE]
      new_risk_prediction(
        id = last_slice$id,
        risk = 1 - last_slice$surv,
        model = last_slice$model,
        time = last_time,
        ensemble = FALSE,
        set = set
      )
    },

    model_info = function() {
      info <- super$model_info()
      info$label <- "Random forest survival"
      info
    },

    required_packages = function() {
      c("randomForestSRC")
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

.rf_align_survival <- function(surv_matrix, base_times, target_times) {
  n_obs <- nrow(surv_matrix)
  out <- matrix(NA_real_, nrow = length(target_times), ncol = n_obs)
  for (i in seq_len(n_obs)) {
    surv_vec <- surv_matrix[i, ]
    interpolator <- stats::approxfun(
      x = base_times,
      y = surv_vec,
      method = "linear",
      yleft = 1,
      yright = tail(surv_vec, 1)
    )
    out[, i] <- interpolator(target_times)
  }
  out
}

.register_time_to_event_model(
  engine = "random_forest",
  outcome = "survival",
  constructor = function(spec = list()) {
    defaults <- list(engine = "random_forest", package = "randomForestSRC")
    RandomForestSurvivalModel$new(spec = modifyList(defaults, spec))
  },
  packages = "randomForestSRC",
  tags = c("tree-based", "ensemble"),
  label = "Random forest survival"
)
