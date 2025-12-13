#' Internal helpers for registering competing-risk engines
#'
#' These wrappers adapt the existing `CRModel_*` fitters into the new R6-based
#' infrastructure without duplicating fitting logic.
#'
#' @keywords internal
.competing_risk_align_cif <- function(prediction, target_times) {
  cifs <- prediction$CIFs
  if (is.null(cifs)) cifs <- prediction$cifs
  if (is.null(cifs)) cifs <- prediction$Probs
  times <- prediction$Times
  if (is.null(times)) times <- prediction$times
  if (is.null(cifs) || is.null(times)) {
    rlang::abort("Competing-risk prediction must contain matrices under `CIFs` (or similar) and `Times`.")
  }
  times <- as.numeric(times)
  if (length(target_times) > 0) {
    target_times <- as.numeric(target_times)
    if (length(times) == length(target_times) && isTRUE(all.equal(times, target_times))) {
      matrix_out <- cifs
      times_out <- target_times
    } else {
      matrix_out <- cifMatInterpolator(
        probsMat = cifs,
        times = times,
        new_times = target_times
      )
      times_out <- target_times
    }
  } else {
    matrix_out <- cifs
    times_out <- times
  }
  list(matrix = matrix_out, times = times_out)
}

.register_competing_risk_engine <- function(engine,
                                            fit_fun_name,
                                            predict_fun_name,
                                            packages = character(),
                                            label = NULL,
                                            tags = character()) {
  fit_fun <- get(fit_fun_name, envir = asNamespace("ml4time2event"))
  predict_fun <- get(predict_fun_name, envir = asNamespace("ml4time2event"))
  classname <- paste0("CompetingRisk_", gsub("[^A-Za-z0-9]", "_", engine))

  CompetingRiskClass <- R6::R6Class(
    classname = classname,
    inherit = CompetingRiskModel,
    public = list(
      model = NULL,
      time_grid = NULL,
      task = NULL,
      cause_codes = NULL,
      model_per_cause = FALSE,

      fit = function(task, time_grid, ...) {
        super$fit(task = task, ...)

        # Ensure engine is set in spec (defensive check)
        if (is.null(self$spec$engine)) {
          self$spec$engine <- engine
        }

        spec_args <- private$spec_args()
        cause_map <- task$metadata$cause_map
        if (is.null(cause_map) || nrow(cause_map) == 0) {
          rlang::abort("Competing-risk tasks must contain a non-empty `cause_map` in metadata.")
        }
        data <- as.data.frame(task$data)
        cause_labels <- as.character(cause_map$cause)
        cause_codes <- cause_map$code
        cause_codes_named <- setNames(cause_codes, cause_labels)
        args <- c(
          list(
            data = data,
            expvars = task$features,
            timevar = task$time_col,
            eventvar = task$event_col
          ),
          spec_args
        )
        engine_name <- self$spec$engine %||% NA_character_
        event_param <- NULL
        fit_formals <- names(formals(fit_fun))
        if ("event_codes" %in% fit_formals) {
          event_param <- "event_codes"
        } else if ("event_of_interest" %in% fit_formals) {
          event_param <- "event_of_interest"
        }

        requires_single <- !is.na(engine_name) && engine_name %in% private$single_cause_engines
        multi_fit <- isTRUE(requires_single) && length(cause_codes) > 1 && !is.null(event_param)

        if (multi_fit) {
          models <- vector("list", length(cause_codes))
          names(models) <- cause_labels
          for (idx in seq_along(cause_codes)) {
            args_single <- args
            args_single[[event_param]] <- cause_codes[[idx]]
            models[[idx]] <- do.call(fit_fun, args_single)
          }
          self$model <- models
          self$model_per_cause <- TRUE
        } else {
          if (!is.null(event_param)) {
            args[[event_param]] <- cause_codes
          }
          model_obj <- do.call(fit_fun, args)
          self$model <- model_obj
          self$model_per_cause <- FALSE
        }

        self$cause_codes <- cause_codes_named
        self$time_grid <- time_grid
        self$task <- task
        invisible(self)
      },

      predict_cif = function(newdata, times, set = "test", ...) {
        private$ensure_fitted()
        model_name <- self$spec$engine
        complete_data <- .ensure_prediction_data(newdata, self$task)
        feature_frame <- as.data.frame(complete_data[, self$task$features, drop = FALSE])
        complete_idx <- stats::complete.cases(feature_frame)
        ids <- complete_data[[self$task$id_col]]

        if (is.null(times)) {
          target_times <- self$time_grid
        } else {
          target_times <- sort(unique(as.numeric(times)))
        }
        if (length(target_times) == 0) {
          rlang::abort("`times` must be numeric and non-empty.")
        }

        preds_complete <- NULL
        cause_labels <- names(self$cause_codes)
        if (is.null(cause_labels) || length(cause_labels) == 0) {
          rlang::abort("Competing-risk model does not contain any cause mappings.")
        }
        if (any(complete_idx)) {
          newdata_complete <- feature_frame[complete_idx, , drop = FALSE]
          id_complete <- ids[complete_idx]
          cause_predictions <- lapply(cause_labels, function(cause_label) {
            cause_code <- self$cause_codes[[cause_label]]
            model_obj <- private$get_model_for_cause(cause_label, cause_code)
            if (is.null(model_obj)) {
              rlang::warn(glue::glue(
                "No fitted model available for cause '{cause_label}' in engine '{model_name}'. Returning NA predictions."
              ))
            }
            predict_args <- list(
              modelout = model_obj,
              newdata = newdata_complete
            )
            if ("new_times" %in% names(formals(predict_fun))) {
              predict_args$new_times <- target_times
            }
            if ("event_of_interest" %in% names(formals(predict_fun))) {
              predict_args$event_of_interest <- cause_code
            } else if ("event_codes" %in% names(formals(predict_fun))) {
              predict_args$event_codes <- cause_code
            }
            if (is.null(model_obj)) {
              matrix <- matrix(NA_real_, nrow = length(target_times), ncol = sum(complete_idx))
              times_out <- target_times
            } else {
              pred <- do.call(predict_fun, predict_args)
              aligned <- .competing_risk_align_cif(pred, target_times)
              matrix <- aligned$matrix
              times_out <- aligned$times
            }
            new_cif_prediction(
              id = rep(id_complete, each = length(times_out)),
              time = rep(times_out, times = length(id_complete)),
              cause = rep(cause_label, length(id_complete) * length(times_out)),
              cif = as.vector(matrix),
              model = rep(model_name, length(id_complete) * length(times_out)),
              ensemble = FALSE,
              set = set
            )
          })
          preds_complete <- dplyr::bind_rows(cause_predictions)
        }

        preds_missing <- NULL
        if (!all(complete_idx)) {
          missing_ids <- ids[!complete_idx]
          rlang::warn(glue::glue(
            "Omitting {length(missing_ids)} rows with missing predictors for engine '{model_name}'."
          ))
          missing_list <- lapply(cause_labels, function(cause_label) {
            new_cif_prediction(
              id = rep(missing_ids, each = length(target_times)),
              time = rep(target_times, times = length(missing_ids)),
              cause = rep(cause_label, length(missing_ids) * length(target_times)),
              cif = rep(NA_real_, length(missing_ids) * length(target_times)),
              model = rep(model_name, length(missing_ids) * length(target_times)),
              ensemble = FALSE,
              set = set
            )
          })
          preds_missing <- dplyr::bind_rows(missing_list)
        }

        pieces <- list(preds_complete, preds_missing)
        pieces <- pieces[!vapply(pieces, is.null, logical(1))]
        if (length(pieces) == 0) {
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
        result <- dplyr::bind_rows(pieces)
        result <- .cif_enforce_bounds(result)
        result
      },

      model_info = function() {
        info <- super$model_info()
        if (is.null(label)) {
          info$label <- engine
        } else {
          info$label <- label
        }
        info
      },

      required_packages = function() {
        packages
      }
    ),
    private = list(
      single_cause_engines = c(
        "fine_gray",
        "cr_rulefit",
        "cr_shallownn",
        "cr_random_forest",
        "cr_bart"
      ),

      ensure_fitted = function() {
        if (!isTRUE(self$fitted)) {
          rlang::abort("Model must be fitted before predictions can be generated.")
        }
      },

      spec_args = function() {
        spec <- self$spec %||% list()
        spec[c("engine", "package", "outcome")] <- NULL
        spec
      },

      get_model_for_cause = function(cause_label, cause_code) {
        if (!isTRUE(self$model_per_cause)) {
          return(self$model)
        }
        models <- self$model
        if (is.null(models)) {
          return(NULL)
        }
        if (!is.null(models[[cause_label]])) {
          return(models[[cause_label]])
        }
        cause_key <- as.character(cause_code)
        if (!is.null(models[[cause_key]])) {
          return(models[[cause_key]])
        }
        NULL
      }
    )
  )

  .register_time_to_event_model(
    engine = engine,
    outcome = "competing_risk",
    constructor = function(spec = list()) {
      CompetingRiskClass$new(spec = modifyList(list(engine = engine), spec))
    },
    packages = packages,
    tags = tags,
    label = if (is.null(label)) engine else label
  )
}

.cif_enforce_bounds <- function(tbl) {
  if (nrow(tbl) == 0) {
    return(tbl)
  }
  tbl <- dplyr::arrange(tbl, id, time, cause)
  ids <- unique(tbl$id)
  out <- vector("list", length(ids))

  for (i in seq_along(ids)) {
    id_val <- ids[i]
    subset_tbl <- tbl[tbl$id == id_val, , drop = FALSE]
    cause_levels <- unique(subset_tbl$cause)
    prev_vals <- setNames(rep(0, length(cause_levels)), cause_levels)
    remaining <- 1
    res <- list()

    for (time_val in unique(subset_tbl$time)) {
      slice <- subset_tbl[subset_tbl$time == time_val, , drop = FALSE]
      if (all(is.na(slice$cif))) {
        res[[length(res) + 1]] <- slice
        next
      }
      deltas <- slice$cif - prev_vals[slice$cause]
      deltas[is.na(deltas)] <- 0
      deltas[deltas < 0] <- 0
      delta_sum <- sum(deltas)
      if (delta_sum > 0 && remaining > 0) {
        scale <- if (delta_sum > remaining) remaining / delta_sum else 1
        adj <- deltas * scale
        prev_vals[slice$cause] <- prev_vals[slice$cause] + adj
        remaining <- max(0, remaining - sum(adj))
      }
      slice$cif <- prev_vals[slice$cause]
      res[[length(res) + 1]] <- slice
    }
    out[[i]] <- dplyr::bind_rows(res)
  }

  dplyr::bind_rows(out)
}

.register_competing_risk_engine(
  engine = "cox",
  fit_fun_name = "CRModel_Cox",
  predict_fun_name = "Predict_CRModel_Cox",
  packages = "survival",
  label = "Cause-specific Cox",
  tags = c("cox", "competing-risk")
)

.register_competing_risk_engine(
  engine = "fine_gray",
  fit_fun_name = "CRModel_FineGray",
  predict_fun_name = "Predict_CRModel_FineGray",
  packages = "fastcmprsk",
  label = "Fine-Gray",
  tags = c("semi-parametric")
)

.register_competing_risk_engine(
  engine = "cr_bart",
  fit_fun_name = "CRModel_BART",
  predict_fun_name = "Predict_CRModel_BART",
  packages = "BART",
  label = "BART competing risk",
  tags = c("bart")
)

.register_competing_risk_engine(
  engine = "cr_gam",
  fit_fun_name = "CRModel_GAM",
  predict_fun_name = "Predict_CRModel_GAM",
  packages = "mgcv",
  label = "GAM competing risk",
  tags = c("smooth", "semiparametric")
)

.register_competing_risk_engine(
  engine = "cr_rulefit",
  fit_fun_name = "CRModel_rulefit",
  predict_fun_name = "Predict_CRModel_rulefit",
  packages = c("rpart", "partykit", "glmnet", "survival"),
  label = "RuleFit competing risk",
  tags = c("rules", "tree-based")
)

.register_competing_risk_engine(
  engine = "cr_random_forest",
  fit_fun_name = "CRModel_RF",
  predict_fun_name = "Predict_CRModel_RF",
  packages = "randomForestSRC",
  label = "Random forest competing risk",
  tags = c("tree-based", "ensemble")
)

.register_competing_risk_engine(
  engine = "cr_survreg",
  fit_fun_name = "CRModel_SurvReg",
  predict_fun_name = "Predict_CRModel_SurvReg",
  packages = "survival",
  label = "Parametric competing risk",
  tags = c("parametric")
)

.register_competing_risk_engine(
  engine = "cr_shallownn",
  fit_fun_name = "CRModel_ShallowNN",
  predict_fun_name = "Predict_CRModel_ShallowNN",
  packages = character(),
  label = "Shallow NN competing risk",
  tags = c("neural-network")
)

.register_competing_risk_engine(
  engine = "cr_ttah",
  fit_fun_name = "CRModel_TTAH",
  predict_fun_name = "Predict_CRModel_TTAH",
  packages = character(),
  label = "TTAH competing risk",
  tags = c("semiparametric")
)

.register_competing_risk_engine(
  engine = "cr_xgboost",
  fit_fun_name = "CRModel_xgboost",
  predict_fun_name = "Predict_CRModel_xgboost",
  packages = "xgboost",
  label = "XGBoost competing risk",
  tags = c("tree-based", "gradient-boosting")
)
