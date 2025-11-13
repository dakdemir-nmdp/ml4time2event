#' Cox proportional hazards model (R6 implementation)
#'
#' This class wraps `survival::coxph` and provides the interface required by the
#' new `ml4time2event` architecture.
#'
#' @keywords internal
#' @noRd
CoxSurvivalModel <- R6::R6Class(
  classname = "CoxSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    basehaz_function = NULL,
    time_grid = NULL,
    task = NULL,

    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      data <- task$data
      required <- c(task$time_col, task$event_col, task$features)
      data <- data[, required, drop = FALSE]
      data <- stats::na.omit(data)

      outcome <- survival::Surv(time = data[[task$time_col]], event = data[[task$event_col]])
      formula <- stats::as.formula(paste("outcome ~", paste(task$features, collapse = "+")))
      fit <- survival::coxph(formula, data = data, ties = "efron", x = TRUE)
      self$model <- fit
      self$time_grid <- time_grid
      self$task <- task

      basehaz <- survival::basehaz(fit, centered = FALSE)
      basehaz <- basehaz[order(basehaz$time), , drop = FALSE]
      basehaz_time <- basehaz$time
      basehaz_hazard <- basehaz$hazard
      if (length(basehaz_time) == 0) {
        rlang::abort("Baseline hazard could not be computed for the Cox model.")
      }
      if (isTRUE(all.equal(basehaz_time[1], 0))) {
        basehaz_hazard[1] <- 0
      } else {
        basehaz_time <- c(0, basehaz_time)
        basehaz_hazard <- c(0, basehaz_hazard)
      }
      self$basehaz_function <- stats::approxfun(
        x = basehaz_time,
        y = basehaz_hazard,
        method = "linear",
        rule = 2,
        yleft = 0
      )

      self$fitted <- TRUE
      invisible(self)
    },

    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      newdata <- .ensure_prediction_data(newdata, self$task)
      times <- sort(unique(times))
      H0 <- self$basehaz_function(times)

      model_name <- self$spec$engine
      if (is.null(model_name)) {
        model_name <- "cox"
      }

      complete_idx <- stats::complete.cases(newdata)
      if (!all(complete_idx)) {
        missing_ids <- newdata[[self$task$id_col]][!complete_idx]
        rlang::warn(glue::glue(
          "Omitting {length(missing_ids)} rows with missing predictors for engine '{model_name}'."
        ))
      }

      newdata_complete <- newdata[complete_idx, , drop = FALSE]
      id_values_complete <- newdata_complete[[self$task$id_col]]
      lp <- stats::predict(self$model, newdata = newdata_complete, type = "lp")
      risk_factor <- exp(lp)
      surv_matrix <- exp(-outer(H0, risk_factor))

      id_values <- newdata[[self$task$id_col]]
      preds_complete <- new_survival_prediction(
        id = rep(id_values_complete, each = length(times)),
        time = rep(times, times = length(id_values_complete)),
        surv = as.vector(surv_matrix),
        model = rep(model_name, length(id_values_complete) * length(times)),
        ensemble = FALSE,
        set = set
      )
      if (all(complete_idx)) {
        return(preds_complete)
      }

      missing_ids <- id_values[!complete_idx]
      preds_missing <- new_survival_prediction(
        id = rep(missing_ids, each = length(times)),
        time = rep(times, times = length(missing_ids)),
        surv = rep(NA_real_, length(missing_ids) * length(times)),
        model = rep(model_name, length(missing_ids) * length(times)),
        ensemble = FALSE,
        set = set
      )
      dplyr::bind_rows(preds_complete, preds_missing)
    },

    predict_risk = function(newdata, times = NULL, set = "test", ...) {
      private$ensure_fitted()
      newdata <- .ensure_prediction_data(newdata, self$task)
      model_name <- self$spec$engine
      if (is.null(model_name)) {
        model_name <- "cox"
      }
      complete_idx <- stats::complete.cases(newdata)
      if (!all(complete_idx)) {
        missing_ids <- newdata[[self$task$id_col]][!complete_idx]
        rlang::warn(glue::glue(
          "Omitting {length(missing_ids)} rows with missing predictors for engine '{model_name}'."
        ))
      }
      newdata_complete <- newdata[complete_idx, , drop = FALSE]
      id_values_complete <- newdata_complete[[self$task$id_col]]
      lp <- stats::predict(self$model, newdata = newdata_complete, type = "lp")
      preds_complete <- new_risk_prediction(
        id = id_values_complete,
        risk = exp(lp),
        model = model_name,
        time = if (is.null(times)) tail(self$time_grid, 1) else max(times),
        ensemble = FALSE,
        set = set
      )
      if (all(complete_idx)) {
        return(preds_complete)
      }
      missing_ids <- newdata[[self$task$id_col]][!complete_idx]
      preds_missing <- new_risk_prediction(
        id = missing_ids,
        risk = rep(NA_real_, length(missing_ids)),
        model = model_name,
        time = if (is.null(times)) tail(self$time_grid, 1) else max(times),
        ensemble = FALSE,
        set = set
      )
      dplyr::bind_rows(preds_complete, preds_missing)
    },

    model_info = function() {
      info <- super$model_info()
      info$label <- "Cox proportional hazards"
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

.ensure_prediction_data <- function(newdata, task) {
  if (!task$id_col %in% colnames(newdata)) {
    newdata[[task$id_col]] <- seq_len(nrow(newdata))
  }
  missing <- setdiff(task$features, colnames(newdata))
  if (length(missing) > 0) {
    rlang::abort(paste0("Missing feature columns: ", paste(missing, collapse = ", ")))
  }
  newdata[, c(task$id_col, task$features), drop = FALSE]
}

.register_time_to_event_model(
  engine = "cox",
  outcome = "survival",
  constructor = function(spec = list()) {
    defaults <- list(engine = "cox", package = "survival")
    CoxSurvivalModel$new(spec = modifyList(defaults, spec))
  },
  packages = "survival",
  tags = c("cox", "parametric"),
  label = "Cox proportional hazards"
)
