#' Fit time-to-event models via a unified interface
#'
#' @param task An object created with `ml4t2e_task_surv()` or `ml4t2e_task_cr()`.
#' @param models Character vector of engine names to fit. Use
#'   `ml4t2e_list_models()` to discover available options.
#' @param ensemble Controls ensembling strategy. Currently `"auto"` and `TRUE`
#'   are treated as aliases for stacking; `"simple"` performs an unweighted
#'   average. Set to `FALSE` to skip ensembling. (Only `FALSE` is implemented in
#'   this initial refactor.)
#' @param recipe Optional `recipes` recipe for preprocessing. Not yet supported
#'   in this refactor; placeholder for future work.
#' @param resampling Optional resampling specification from `rsample`.
#' @param seed Optional integer seed for reproducibility.
#' @param controls Named list of engine-specific parameter lists. Use entries
#'   like `controls = list(cox = list(penalized = TRUE))`.
#'
#' @return An object of class `t2e_fit`.
#' @export
ml4t2e_fit <- function(task,
                       models = "cox",
                       ensemble = c(FALSE, TRUE, "auto", "stack", "simple"),
                       recipe = NULL,
                       resampling = NULL,
                       seed = NULL,
                       controls = list()) {
  if (!inherits(task, "t2e_task")) {
    rlang::abort("`task` must be created with `ml4t2e_task_*()`.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      rlang::abort("`seed` must be a single numeric value.")
    }
    seed <- as.integer(seed)
    set.seed(seed)
  }

  ensemble_choice <- match.arg(as.character(ensemble[[1]]), c("FALSE", "TRUE", "auto", "stack", "simple"))
  ensemble_mode <- ensemble_choice
  if (ensemble_choice %in% c("TRUE", "auto")) {
    ensemble_mode <- "simple"
  }
  ensemble_mode <- tolower(ensemble_mode)
  ensemble_enabled <- !identical(ensemble_choice, "FALSE")

  outcome_type <- attr(task, "task_type")
  models <- unique(models)
  if (length(models) == 0) {
    rlang::abort("`models` must contain at least one engine name.")
  }

  time_grid <- .resolve_time_grid(task, controls)

  fitted_models <- list()
  for (engine in models) {
    spec <- controls[[engine]]
    if (is.null(spec)) {
      spec <- list()
    }
    instance <- .instantiate_model(engine = engine, outcome = outcome_type, spec = spec)
    trained <- instance$fit(task = task, time_grid = time_grid)
    if (inherits(trained, "R6")) {
      fitted_models[[engine]] <- trained
    } else {
      fitted_models[[engine]] <- instance
    }
  }

  ensemble_controls <- controls$ensemble
  ensemble_obj <- NULL
  if (ensemble_enabled) {
    ensemble_obj <- .build_ensembler(
      models = fitted_models,
      model_names = names(fitted_models),
      outcome_type = outcome_type,
      strategy = ensemble_mode,
      time_grid = time_grid,
      task = task,
      controls = ensemble_controls
    )
  }

  fit <- list(
    task = task,
    models = fitted_models,
    model_names = names(fitted_models),
    outcome_type = outcome_type,
    time_grid = time_grid,
    ensemble = ensemble_obj,
    ensemble_strategy = if (!is.null(ensemble_obj)) ensemble_obj$strategy else "none",
    seed = seed,
    controls = controls,
    recipe = recipe,
    resampling = resampling,
    created = Sys.time(),
    api_version = "0.1.0"
  )
  if (!is.null(ensemble_obj)) {
    fit$model_names <- c(fit$model_names, "ensemble")
  }
  class(fit) <- c("t2e_fit", "list")
  fit
}

.resolve_time_grid <- function(task, controls) {
  grid <- controls$times
  if (!is.null(grid)) {
    grid <- sort(unique(as.numeric(grid)))
    if (anyNA(grid)) {
      rlang::abort("`controls$times` must be numeric without NA values.")
    }
    if (!any(grid == 0)) {
      grid <- c(0, grid)
      grid <- sort(unique(grid))
    }
    return(grid)
  }
  range <- task$time_range
  if (!is.numeric(range) || length(range) != 2) {
    rlang::abort("Task `time_range` must be numeric of length 2.")
  }
  if (range[2] <= range[1]) {
    range[2] <- range[1] + 1
  }
  sort(unique(c(0, seq(range[1], range[2], length.out = 100))))
}

.build_ensembler <- function(models,
                             model_names,
                             outcome_type,
                             strategy,
                             time_grid,
                             task,
                             controls = list()) {
  strategy <- tolower(strategy %||% "simple")
  if (identical(strategy, "false")) {
    return(NULL)
  }

  if (identical(outcome_type, "survival")) {
    if (identical(strategy, "stack")) {
      rlang::warn("Stacked ensembling is not yet implemented; using simple averaging instead.")
      strategy <- "simple"
    }
    return(SurvivalEnsembler$new(
      models = models,
      model_names = model_names,
      outcome_type = outcome_type,
      strategy = strategy,
      time_grid = time_grid,
      controls = controls
    ))
  }

  if (identical(outcome_type, "competing_risk")) {
    if (identical(strategy, "stack")) {
      rlang::warn("Stacked ensembling is not yet implemented; using simple averaging instead.")
      strategy <- "simple"
    }
    cause_map <- task$metadata$cause_map
    if (is.null(cause_map) || nrow(cause_map) == 0) {
      rlang::warn("Competing risk ensemble could not determine a cause map; skipping ensemble construction.")
      return(NULL)
    }
    return(CompetingRiskEnsembler$new(
      models = models,
      model_names = model_names,
      outcome_type = outcome_type,
      strategy = strategy,
      time_grid = time_grid,
      cause_map = cause_map,
      controls = controls
    ))
  }

  NULL
}

#' @export
predict.t2e_fit <- function(object,
                            newdata = NULL,
                            times = NULL,
                            type = NULL,
                            include = "all",
                            ...) {
  if (!inherits(object, "t2e_fit")) {
    rlang::abort("`object` must be a `t2e_fit`.")
  }
  outcome_type <- object[["outcome_type"]]
  if (is.null(type)) {
    type <- if (identical(outcome_type, "survival")) "survival" else "cif"
  }
  type_choices <- if (identical(outcome_type, "survival")) {
    c("survival", "risk")
  } else {
    c("cif")
  }
  type <- match.arg(type, type_choices)
  data <- newdata
  set_label <- "test"
  if (is.null(data)) {
    task <- object[["task"]]
    data <- task$data
    set_label <- "train"
  } else {
    data <- dplyr::as_tibble(data)
    .verify_new_data_columns(data, object[["task"]])
  }
  time_grid <- object[["time_grid"]]
  times <- .resolve_prediction_times(times, time_grid)

  include <- .normalize_include(include, object[["model_names"]])
  model_store <- object[["models"]]
  predictions <- list()
  for (engine in include) {
    if (identical(engine, "ensemble")) {
      ensemble_obj <- object$ensemble
      if (is.null(ensemble_obj)) {
        rlang::warn("No ensemble available; skipping ensemble predictions.")
        next
      }
      if (identical(outcome_type, "survival")) {
        if (identical(type, "survival")) {
          pred <- ensemble_obj$predict_survival(newdata = data, times = times, set = set_label)
        } else {
          pred <- ensemble_obj$predict_risk(newdata = data, times = times, set = set_label)
        }
      } else {
        pred <- ensemble_obj$predict_cif(newdata = data, times = times, set = set_label)
      }
    } else {
      model_obj <- model_store[[engine]]
      if (is.null(model_obj)) {
        rlang::warn(glue::glue("Model '{engine}' not found; skipping."))
        next
      }
      if (identical(outcome_type, "survival")) {
        if (identical(type, "survival")) {
          pred <- model_obj$predict_survival(newdata = data, times = times, set = set_label)
        } else {
          pred <- model_obj$predict_risk(newdata = data, times = times, set = set_label)
        }
      } else {
        pred <- model_obj$predict_cif(newdata = data, times = times, set = set_label)
      }
    }
    predictions[[engine]] <- pred
  }

  if (length(predictions) == 0) {
    rlang::abort("No predictions were generated.")
  }

  result <- dplyr::bind_rows(predictions)
  attr(result, "pred_type") <- type
  attr(result, "time_grid") <- times
  class(result) <- unique(c("t2e_pred", class(result)))
  result
}

.normalize_include <- function(include, models) {
  if (identical(include, "all")) {
    return(models)
  }
  include <- unique(include)
  include <- include[include %in% models]
  if (length(include) == 0) {
    rlang::abort("`include` did not match any fitted models.")
  }
  include
}

.resolve_prediction_times <- function(times, default) {
  if (is.null(times)) {
    return(default)
  }
  times <- sort(unique(as.numeric(times)))
  if (anyNA(times)) {
    rlang::abort("`times` must be numeric without NA values.")
  }
  times
}

.verify_new_data_columns <- function(newdata, task) {
  required <- unique(task$features)
  missing <- setdiff(required, colnames(newdata))
  if (length(missing) > 0) {
    rlang::abort(paste0("`newdata` is missing columns: ", paste(missing, collapse = ", ")))
  }
}

#' @export
print.t2e_fit <- function(x, ...) {
  cat("<t2e_fit>\n")
  cat(" outcome:", x[["outcome_type"]], "\n", sep = " ")
  cat(" models :", paste(x[["model_names"]], collapse = ", "), "\n", sep = " ")
  cat(" time grid points:", length(x[["time_grid"]]), "\n")
  cat(" created:", format(x[["created"]], usetz = TRUE), "\n")
  invisible(x)
}

#' @export
summary.t2e_fit <- function(object, ...) {
  data.frame(
    model = object[["model_names"]],
    outcome = object[["outcome_type"]],
    n_times = length(object[["time_grid"]]),
    stringsAsFactors = FALSE
  )
}

#' Evaluate fitted models on metrics
#'
#' Currently supports the concordance index (`"c_index"`) and integrated Brier
#' score (`"ibs"`). The latter is computed using a simple trapezoidal integration
#' over the supplied time grid.
#'
#' @param fit_or_preds A `t2e_fit` object or predictions returned by
#'   `predict.t2e_fit()`.
#' @param task Optional task providing ground truth if `fit_or_preds` contains
#'   predictions on new data.
#' @param metrics Character vector of metrics to compute.
#' @param include Which models to include (same semantics as `predict()`).
#'
#' @return A tibble with columns `model`, `metric`, and `value`.
#' @export
ml4t2e_evaluate <- function(fit_or_preds,
                            task = NULL,
                            metrics = c("c_index", "ibs"),
                            include = "all") {
  predictions <- fit_or_preds
  truth_task <- task

  if (inherits(fit_or_preds, "t2e_fit")) {
    outcome_type <- fit_or_preds$outcome_type
    pred_type <- if (identical(outcome_type, "competing_risk")) "cif" else "survival"
    predictions <- predict(fit_or_preds, include = include, type = pred_type)
    truth_task <- fit_or_preds[["task"]]
  }

  if (!inherits(predictions, "t2e_pred")) {
    rlang::abort("`fit_or_preds` must be a `t2e_fit` or predictions from `predict()`.")
  }
  if (is.null(truth_task)) {
    rlang::abort("`task` must be supplied when evaluating detached predictions.")
  }
  task_type <- attr(truth_task, "task_type")
  if (identical(task_type, "competing_risk")) {
    return(.evaluate_competing_risk(predictions, truth_task, metrics))
  }

  metrics <- unique(metrics)
  allowed_metrics <- c("c_index", "ibs")
  bad <- setdiff(metrics, allowed_metrics)
  if (length(bad) > 0) {
    rlang::abort(paste0("Unsupported metrics: ", paste(bad, collapse = ", ")))
  }

  task_df <- truth_task[["data"]]
  id_col <- truth_task[["id_col"]]
  time_col <- truth_task[["time_col"]]
  event_col <- truth_task[["event_col"]]

  results <- list()
  for (model_name in unique(predictions[["model"]])) {
    pred_subset <- predictions[predictions[["model"]] == model_name, , drop = FALSE]
    metric_values <- list()
    if ("c_index" %in% metrics) {
      # Always use Expected Time Lost (ETL) for concordance
      lp <- .extract_etl_from_survival(pred_subset, task_df, id_col)
      surv_obj <- survival::Surv(task_df[[time_col]], task_df[[event_col]])
      valid_idx <- !is.na(lp)
      if (sum(valid_idx) == 0) {
        metric_values[["c_index"]] <- NA_real_
      } else {
        valid_df <- task_df[valid_idx, , drop = FALSE]
        surv_obj_valid <- survival::Surv(valid_df[[time_col]], valid_df[[event_col]])
        c_val <- tryCatch({
          concordance <- survival::concordance(surv_obj_valid ~ -lp[valid_idx])
          unname(concordance$concordance)
        }, error = function(e) {
          concordance <- survival::concordance(surv_obj_valid ~ lp[valid_idx])
          1 - unname(concordance$concordance)
        })
        metric_values[["c_index"]] <- c_val
      }
    }
    if ("ibs" %in% metrics) {
      ibs <- .integrated_brier(pred_subset, task_df, id_col, time_col, event_col)
      metric_values[["ibs"]] <- ibs
    }
    for (metric_name in names(metric_values)) {
      results[[length(results) + 1]] <- data.frame(
        model = model_name,
        metric = metric_name,
        value = metric_values[[metric_name]],
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(results)
}

.extract_risk_from_predictions <- function(risk_predictions, task_df, id_col) {
  # Extract risk scores from risk predictions
  # Match by id to ensure correct order
  risk_scores <- numeric(nrow(task_df))
  risk_scores[] <- NA_real_
  
  id_map <- match(task_df[[id_col]], risk_predictions$id)
  valid_idx <- !is.na(id_map)
  risk_scores[valid_idx] <- risk_predictions$risk[id_map[valid_idx]]
  
  risk_scores
}

.extract_etl_from_survival <- function(predictions, task_df, id_col) {
  # Calculate Expected Time Lost (ETL) from survival predictions
  # ETL = ∫₀^T (1 - S(t)) dt where T is the maximum prediction time
  # Higher ETL = higher risk (more expected time lost due to event)
  
  base_preds <- predictions[, c("id", "time", "surv")]
  unique_times <- sort(unique(base_preds$time))
  max_time <- max(unique_times, na.rm = TRUE)
  
  # Get unique IDs in task_df order
  ids <- task_df[[id_col]]
  unique_ids <- unique(ids)
  
  # Calculate ETL for each observation
  etl_scores <- numeric(length(unique_ids))
  names(etl_scores) <- unique_ids
  
  for (id_val in unique_ids) {
    id_preds <- base_preds[base_preds$id == id_val, ]
    if (nrow(id_preds) == 0) {
      etl_scores[as.character(id_val)] <- NA_real_
      next
    }
    
    # Sort by time
    id_preds <- id_preds[order(id_preds$time), ]
    surv_curve <- id_preds$surv
    time_curve <- id_preds$time
    
    # Calculate event probability: 1 - S(t)
    event_probs <- 1 - surv_curve
    
    # Use unified ETL calculation
    etl <- calculate_expected_time_lost(
      times = time_curve,
      event_probs = event_probs,
      upper_limit = max_time,
      lower_limit = 0
    )
    
    etl_scores[as.character(id_val)] <- etl
  }
  
  # Match to task_df order
  risk_scores <- numeric(nrow(task_df))
  risk_scores[] <- NA_real_
  
  id_map <- match(as.character(ids), names(etl_scores))
  valid_idx <- !is.na(id_map)
  risk_scores[valid_idx] <- etl_scores[id_map[valid_idx]]
  
  risk_scores
}

.extract_etl_from_cif <- function(predictions, task_df, id_col, cause_label = NULL) {
  # Calculate Expected Time Lost (ETL) from CIF predictions
  # ETL = ∫₀^T CIF(t) dt where T is the maximum prediction time
  # Higher ETL = higher risk (more expected time lost due to event)
  
  # Filter by cause if specified
  if (!is.null(cause_label)) {
    base_preds <- predictions[predictions$cause == cause_label, c("id", "time", "cif")]
  } else {
    base_preds <- predictions[, c("id", "time", "cif")]
  }
  
  if (nrow(base_preds) == 0) {
    return(rep(NA_real_, nrow(task_df)))
  }
  
  unique_times <- sort(unique(base_preds$time))
  max_time <- max(unique_times, na.rm = TRUE)
  
  # Get unique IDs in task_df order
  ids <- task_df[[id_col]]
  unique_ids <- unique(ids)
  
  # Calculate ETL for each observation
  etl_scores <- numeric(length(unique_ids))
  names(etl_scores) <- unique_ids
  
  for (id_val in unique_ids) {
    id_preds <- base_preds[base_preds$id == id_val, ]
    if (nrow(id_preds) == 0) {
      etl_scores[as.character(id_val)] <- NA_real_
      next
    }
    
    # Sort by time
    id_preds <- id_preds[order(id_preds$time), ]
    cif_curve <- id_preds$cif
    time_curve <- id_preds$time
    
    # Ensure CIF is bounded [0, 1]
    cif_curve <- pmax(0, pmin(1, cif_curve))
    
    # Use unified ETL calculation
    etl <- calculate_expected_time_lost(
      times = time_curve,
      event_probs = cif_curve,
      upper_limit = max_time,
      lower_limit = 0
    )
    
    etl_scores[as.character(id_val)] <- etl
  }
  
  # Match to task_df order
  risk_scores <- numeric(nrow(task_df))
  risk_scores[] <- NA_real_
  
  id_map <- match(as.character(ids), names(etl_scores))
  valid_idx <- !is.na(id_map)
  risk_scores[valid_idx] <- etl_scores[id_map[valid_idx]]
  
  risk_scores
}

.extract_linear_predictor <- function(predictions, task_df, id_col, time_col, event_col = NULL) {
  if ("risk" %in% colnames(predictions)) {
    # If risk column exists, use it directly
    return(.extract_risk_from_predictions(predictions, task_df, id_col))
  }
  
  # Use Expected Time Lost (ETL) as risk score
  # For survival: ETL = ∫₀^T (1 - S(t)) dt
  # Higher ETL = higher risk (more expected time lost)
  if ("surv" %in% colnames(predictions)) {
    return(.extract_etl_from_survival(predictions, task_df, id_col))
  }
  
  # Should not reach here if predictions are valid
  rlang::abort("Unable to extract risk scores from predictions.")
}

.integrated_brier <- function(predictions, task_df, id_col, time_col, event_col) {
  times <- sort(unique(predictions$time))
  base_preds <- as.data.frame(predictions[, c("id", "time", "surv")])
  surv_mat <- stats::reshape(
    base_preds,
    idvar = "id",
    timevar = "time",
    direction = "wide"
  )
  colnames(surv_mat) <- sub("^surv\\.", "", colnames(surv_mat))
  surv_mat <- surv_mat[match(task_df[[id_col]], surv_mat$id), , drop = FALSE]
  surv_mat <- surv_mat[, -1, drop = FALSE]
  surv_mat <- as.matrix(surv_mat)

  status <- task_df[[event_col]]
  time <- task_df[[time_col]]

  brier_values <- numeric(length(times))
  for (i in seq_along(times)) {
    t <- times[i]
    surv_t <- surv_mat[, i]
    # Only include observations at risk at time t
    # At risk if: time > t (still at risk) OR (time <= t AND status == 1) (event occurred)
    at_risk <- (time > t) | (time <= t & status == 1)
    valid_idx <- at_risk & !is.na(surv_t)
    
    if (!any(valid_idx)) {
      brier_values[i] <- NA_real_
      next
    }
    
    # Brier score: (predicted - observed)^2
    # For events by time t: observed = 0 (no survival)
    # For those still at risk: observed = 1 (still alive)
    observed <- ifelse(time[valid_idx] <= t & status[valid_idx] == 1, 0, 1)
    brier_values[i] <- mean((surv_t[valid_idx] - observed)^2, na.rm = TRUE)
  }
  # Remove NA values before integration
  valid_brier <- !is.na(brier_values)
  if (sum(valid_brier) < 2) {
    return(if (any(valid_brier)) brier_values[valid_brier][1] else NA_real_)
  }
  pracma::trapz(times[valid_brier], brier_values[valid_brier])
}

.evaluate_competing_risk <- function(predictions, task, metrics) {
  metrics <- unique(metrics)
  allowed_metrics <- c("c_index", "ibs")
  bad <- setdiff(metrics, allowed_metrics)
  if (length(bad) > 0) {
    rlang::abort(paste0("Unsupported metrics: ", paste(bad, collapse = ", ")))
  }

  task_df <- task[["data"]]
  id_col <- task[["id_col"]]
  time_col <- task[["time_col"]]
  event_col <- task[["event_col"]]
  cause_map <- task[["metadata"]][["cause_map"]]
  if (is.null(cause_map) || nrow(cause_map) == 0) {
    rlang::abort("Competing risk task metadata must include a cause map.")
  }

  ids <- task_df[[id_col]]
  time_grid <- attr(predictions, "time_grid")
  if (is.null(time_grid)) {
    time_grid <- sort(unique(predictions$time))
  }
  surv_obj <- survival::Surv(task_df[[time_col]], task_df[[event_col]], type = "mstate")

  results <- list()
  model_names <- unique(predictions$model)
  for (model_name in model_names) {
    model_pred <- predictions[predictions$model == model_name, , drop = FALSE]
    for (i in seq_len(nrow(cause_map))) {
      cause_label <- cause_map$cause[i]
      cause_code <- cause_map$code[i]
      cause_pred <- model_pred[model_pred$cause == cause_label, , drop = FALSE]
      if (nrow(cause_pred) == 0) {
        next
      }
      matrix <- .cif_prediction_matrix(cause_pred, ids, time_grid)

      if ("ibs" %in% metrics) {
        ibs_val <- integratedBrierCR(
          SurvObj = surv_obj,
          Predictions = matrix,
          eval_times = time_grid,
          cause = cause_code,
          pred_times = time_grid
        )
        results[[length(results) + 1]] <- data.frame(
          model = model_name,
          cause = cause_label,
          metric = "ibs",
          value = ibs_val,
          stringsAsFactors = FALSE
        )
      }
      if ("c_index" %in% metrics) {
        # Use ETL (Expected Time Lost) as risk score for C-index
        # ETL = ∫₀^T CIF(t) dt, which is more principled than time-dependent C-index
        cause_pred_subset <- cause_pred[cause_pred$cause == cause_label, ]
        if (nrow(cause_pred_subset) > 0) {
          etl_scores <- .extract_etl_from_cif(
            predictions = cause_pred_subset,
            task_df = task_df,
            id_col = id_col,
            cause_label = cause_label
          )
          # Use survival::concordance with ETL as risk score
          # For competing risks, check if there's an event_code column, otherwise use event_col
          event_code_col <- ".event_code"
          if (event_code_col %in% colnames(task_df)) {
            # Use numeric event codes from task metadata
            event_indicator <- ifelse(task_df[[event_code_col]] == cause_code, 1L, 0L)
          } else {
            # Fallback: try to match cause labels/values
            # The event_col should contain numeric codes where 0 = censored, >0 = cause
            event_indicator <- ifelse(task_df[[event_col]] == cause_code, 1L, 0L)
          }
          surv_obj_cause <- survival::Surv(task_df[[time_col]], event_indicator)
          
          # Calculate concordance using -ETL (survival::concordance expects lower = worse outcome)
          valid_idx <- !is.na(etl_scores) & is.finite(etl_scores)
          if (sum(valid_idx) > 0 && sum(event_indicator[valid_idx]) > 0) {
            c_val <- tryCatch({
              concordance <- survival::concordance(surv_obj_cause[valid_idx] ~ -etl_scores[valid_idx])
              unname(concordance$concordance)
            }, error = function(e) {
              concordance <- survival::concordance(surv_obj_cause[valid_idx] ~ etl_scores[valid_idx])
              1 - unname(concordance$concordance)
            })
          } else {
            c_val <- NA_real_
          }
        } else {
          c_val <- NA_real_
        }
        results[[length(results) + 1]] <- data.frame(
          model = model_name,
          cause = cause_label,
          metric = "c_index",
          value = c_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    rlang::abort("No competing risk predictions available for evaluation.")
  }
  dplyr::bind_rows(results)
}

.cif_prediction_matrix <- function(pred_subset, ids, time_grid) {
  pred_clean <- pred_subset %>%
    dplyr::group_by(id, time) %>%
    dplyr::summarise(cif = mean(cif, na.rm = TRUE), .groups = "drop")

  mat <- matrix(NA_real_, nrow = length(time_grid), ncol = length(ids))
  colnames(mat) <- ids
  rownames(mat) <- time_grid

  id_index <- match(pred_clean$id, ids)
  time_index <- match(pred_clean$time, time_grid)
  valid <- !is.na(id_index) & !is.na(time_index)
  mat[cbind(time_index[valid], id_index[valid])] <- pred_clean$cif[valid]
  mat
}

#' Persist fitted objects to disk
#'
#' @param fit A `t2e_fit` object.
#' @param path Path to the file to create.
#'
#' @export
ml4t2e_save <- function(fit, path) {
  if (!inherits(fit, "t2e_fit")) {
    rlang::abort("`fit` must be a `t2e_fit` object.")
  }
  saveRDS(fit, file = path)
  invisible(path)
}

#' Load a persisted fit
#'
#' @param path Path to a saved fit (created by `ml4t2e_save()`).
#'
#' @return A `t2e_fit` object.
#' @export
ml4t2e_load <- function(path) {
  fit <- readRDS(path)
  if (!inherits(fit, "t2e_fit")) {
    rlang::abort("The file does not contain a `t2e_fit` object.")
  }
  fit
}
