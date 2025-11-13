#' Base classes for ensembling time-to-event models
#'
#' Provide reusable infrastructure for combining predictions across multiple
#' fitted engines. Supports both survival and competing-risk outcomes with
#' pluggable strategies (currently simple averaging, with stacking hooks).
#'
#' @keywords internal
#' @noRd
EnsemblerBase <- R6::R6Class(
  classname = "EnsemblerBase",
  public = list(
    models = NULL,
    model_names = NULL,
    outcome_type = NULL,
    strategy = NULL,
    time_grid = NULL,
    cause_map = NULL,
    controls = NULL,

    initialize = function(models,
                          model_names,
                          outcome_type,
                          strategy,
                          time_grid,
                          controls = list(),
                          cause_map = NULL) {
      self$models <- models
      self$model_names <- model_names
      self$outcome_type <- outcome_type
      self$strategy <- tolower(strategy %||% "simple")
      self$time_grid <- time_grid
      self$cause_map <- cause_map
      self$controls <- controls %||% list()
    }
  ),
  private = list(
    ensure_strategy_available = function(allowed, context) {
      if (!self$strategy %in% allowed) {
        rlang::abort(glue::glue(
          "Ensemble strategy '{self$strategy}' is not available for {context} predictions."
        ))
      }
    },

    fallback_to_simple = function() {
      rlang::warn("Stacked ensembling is not yet implemented; using simple averaging instead.")
      self$strategy <- "simple"
    }
  )
)

#' @keywords internal
#' @noRd
SurvivalEnsembler <- R6::R6Class(
  classname = "SurvivalEnsembler",
  inherit = EnsemblerBase,
  public = list(
    predict_survival = function(newdata, times = NULL, set = "test") {
      target_times <- .ensemble_resolve_times(times, self$time_grid)
      pred_list <- lapply(self$models, function(model) {
        model$predict_survival(newdata = newdata, times = target_times, set = set)
      })
      names(pred_list) <- self$model_names

      if (identical(self$strategy, "stack")) {
        private$fallback_to_simple()
      }
      private$ensure_strategy_available("simple", "survival")

      .ensemble_average_survival(pred_list, set_label = set)
    },

    predict_risk = function(newdata, times = NULL, set = "test") {
      surv_tbl <- self$predict_survival(newdata = newdata, times = times, set = set)
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
      last_time <- max(surv_tbl$time)
      last_slice <- surv_tbl[surv_tbl$time == last_time, , drop = FALSE]
      new_risk_prediction(
        id = last_slice$id,
        risk = 1 - last_slice$surv,
        model = last_slice$model,
        time = last_time,
        ensemble = last_slice$ensemble,
        set = last_slice$set
      )
    }
  )
)

#' @keywords internal
#' @noRd
CompetingRiskEnsembler <- R6::R6Class(
  classname = "CompetingRiskEnsembler",
  inherit = EnsemblerBase,
  public = list(
    predict_cif = function(newdata, times = NULL, set = "test") {
      target_times <- .ensemble_resolve_times(times, self$time_grid)
      pred_list <- lapply(self$models, function(model) {
        model$predict_cif(newdata = newdata, times = target_times, set = set)
      })
      names(pred_list) <- self$model_names

      if (identical(self$strategy, "stack")) {
        private$fallback_to_simple()
      }
      private$ensure_strategy_available("simple", "competing-risk")

      .ensemble_average_cif(pred_list, set_label = set)
    }
  )
)

# -------------------------------------------------------------------------
# Helper functions shared by ensemble classes
# -------------------------------------------------------------------------

.ensemble_resolve_times <- function(times, default) {
  if (is.null(times)) {
    target_times <- default
  } else {
    target_times <- sort(unique(as.numeric(times)))
  }
  target_times
}

.ensemble_average_survival <- function(pred_list, set_label) {
  pred_list <- pred_list[!vapply(pred_list, is.null, logical(1))]
  if (length(pred_list) == 0) {
    return(new_survival_prediction(
      id = integer(0),
      time = numeric(0),
      surv = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  combined <- dplyr::bind_rows(pred_list)
  if (nrow(combined) == 0) {
    return(new_survival_prediction(
      id = integer(0),
      time = numeric(0),
      surv = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  aggregated <- combined %>%
    dplyr::group_by(id, time) %>%
    dplyr::summarise(surv = mean(surv, na.rm = TRUE), .groups = "drop")

  if (nrow(aggregated) == 0) {
    return(new_survival_prediction(
      id = integer(0),
      time = numeric(0),
      surv = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  aggregated$surv[is.nan(aggregated$surv)] <- NA_real_
  aggregated$model <- "ensemble"
  aggregated$ensemble <- TRUE
  aggregated$set <- set_label

  new_survival_prediction(
    id = aggregated$id,
    time = aggregated$time,
    surv = aggregated$surv,
    model = aggregated$model,
    ensemble = aggregated$ensemble,
    set = aggregated$set
  )
}

.ensemble_average_cif <- function(pred_list, set_label) {
  pred_list <- pred_list[!vapply(pred_list, is.null, logical(1))]
  if (length(pred_list) == 0) {
    return(new_cif_prediction(
      id = integer(0),
      time = numeric(0),
      cause = character(0),
      cif = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  combined <- dplyr::bind_rows(pred_list)
  if (nrow(combined) == 0) {
    return(new_cif_prediction(
      id = integer(0),
      time = numeric(0),
      cause = character(0),
      cif = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  aggregated <- combined %>%
    dplyr::group_by(id, cause, time) %>%
    dplyr::summarise(cif = mean(cif, na.rm = TRUE), .groups = "drop")

  if (nrow(aggregated) == 0) {
    return(new_cif_prediction(
      id = integer(0),
      time = numeric(0),
      cause = character(0),
      cif = numeric(0),
      model = character(0),
      ensemble = logical(0),
      set = character(0)
    ))
  }

  aggregated$cif[is.nan(aggregated$cif)] <- NA_real_
  aggregated$model <- "ensemble"
  aggregated$ensemble <- TRUE
  aggregated$set <- set_label

  new_cif_prediction(
    id = aggregated$id,
    time = aggregated$time,
    cause = aggregated$cause,
    cif = aggregated$cif,
    model = aggregated$model,
    ensemble = aggregated$ensemble,
    set = aggregated$set
  )
}
