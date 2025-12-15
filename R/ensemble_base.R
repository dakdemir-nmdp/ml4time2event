# Suppress R CMD check NOTEs about variables used in NSE contexts
utils::globalVariables(c("id", "time", "cause", "surv", "cif"))

#' Base classes for ensembling time-to-event models
#'
#' Provide reusable infrastructure for combining predictions across multiple
#' fitted engines. Supports both survival and competing-risk outcomes with
#' pluggable strategies (currently simple averaging, with stacking hooks).
#'
#' @keywords internal
#' @noRd
#' @importFrom magrittr %>%

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
    weights = NULL,
    initialize = function(models,
                          model_names,
                          outcome_type,
                          strategy,
                          time_grid,
                          controls = list(),
                          cause_map = NULL,
                          sl_predictions = NULL,
                          sl_actual = NULL,
                          sl_weights = NULL) {
      self$models <- models
      self$model_names <- model_names
      self$outcome_type <- outcome_type
      self$strategy <- tolower(strategy %||% "simple")
      self$time_grid <- time_grid
      self$cause_map <- cause_map
      self$controls <- controls %||% list()

      if (identical(self$strategy, "stack")) {
        if (!is.null(sl_predictions) && !is.null(sl_actual)) {
          # Use default MSE loss for now, can be parameterized via controls
          loss <- self$controls$sl_loss %||% "mse"
          self$weights <- optimizeSuperLearnerWeights(
            predictions_list = sl_predictions,
            actual_surv = sl_actual,
            loss_type = loss,
            weights_matrix = sl_weights
          )
        } else {
          rlang::warn("Stacking requested but no training predictions provided. Falling back to simple averaging.")
          self$strategy <- "simple"
        }
      }
    }
  ),
  private = list(
    ensure_strategy_available = function(allowed, context) {
      if (!self$strategy %in% allowed) {
        rlang::abort(glue::glue(
          "Ensemble strategy '{self$strategy}' is not available for {context} predictions."
        ))
      }
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
        # Convert tidy preds to matrices
        mat_list <- list()
        valid_models <- names(pred_list)
        for (m in valid_models) {
          p <- pred_list[[m]]
          if (nrow(p) > 0) {
            # Pivot to wide matrix
            wide <- ml4t2e_reshape_preds_to_matrix(p, newdata, target_times, "surv")
            mat_list[[m]] <- wide
          }
        }

        # Apply weighted averaging
        w_probs <- survprobMatWeightedAveraging(mat_list, self$weights)

        # Convert back to tidy
        # Reuse internal helper logic but for matrix
        # Create output dataframe manually
        # Row-major matching of wide matrix:
        # Rows = Obs, Cols = Times.
        # Vectorize: c(Mat) -> Time 1 (all obs), Time 2 (all obs)...
        # But tidy format is usually: Obs 1 (all times), Obs 2.
        # Let's align carefully.

        # wide matrix cols are times. rows are obs.
        # as.vector(t(w_probs)) -> Obs 1 (t1..tk), Obs 2...

        n_obs <- nrow(newdata)
        n_times <- length(target_times)

        # Ensure ID column exists
        ids <- if ("id" %in% names(newdata)) newdata$id else seq_len(n_obs)

        flat_surv <- as.vector(t(w_probs))

        return(new_survival_prediction(
          id = rep(ids, each = n_times),
          time = rep(target_times, times = n_obs),
          surv = flat_surv,
          model = rep("ensemble", length(flat_surv)),
          ensemble = TRUE,
          set = set
        ))
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
        # CR Stacking: Weighted average of CIFs
        # Need to do this PER CAUSE?
        # Current optimizeSuperLearnerWeights is single-outcome.
        # For CR, "sl_training_predictions" and "sl_actual" passed to initialize
        # should probably be a list per cause?
        # Or we learn one set of weights for the primary cause?
        # For simplicity in this iteration: Shared weights or Primary Cause weights.

        # However, `EnsembleBase` gets one `sl_predictions` list.
        # If CR, `sl_predictions` should probably be for the cause of interest?
        # Or list of lists?
        # Let's assume common weights for all causes for now to keep it feasible.

        # Logic:
        # For each cause, extract matrix, avg, combine.

        causes <- as.character(self$cause_map$code)
        res_list <- list()

        for (cause_val in causes) {
          # Extract matrix for this cause from each model
          mat_list <- list()
          for (m in self$model_names) {
            p <- pred_list[[m]]
            p_c <- p[p$cause == cause_val, ]
            if (nrow(p_c) > 0) {
              wide <- ml4t2e_reshape_preds_to_matrix(p_c, newdata, target_times, "cif")
              mat_list[[m]] <- wide
            }
          }

          # Average
          w_cif <- cifMatWeightedAveraging(mat_list, self$weights, type = "CumHaz") # Default

          # Convert to tidy
          n_obs <- nrow(newdata)
          n_times <- length(target_times)
          ids <- if ("id" %in% names(newdata)) newdata$id else seq_len(n_obs)

          flat_cif <- as.vector(t(w_cif))

          res_list[[cause_val]] <- new_cif_prediction(
            id = rep(ids, each = n_times),
            time = rep(target_times, times = n_obs),
            cause = rep(cause_val, length(flat_cif)),
            cif = flat_cif,
            model = rep("ensemble", length(flat_cif)),
            ensemble = TRUE,
            set = set
          )
        }
        return(dplyr::bind_rows(res_list))
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
