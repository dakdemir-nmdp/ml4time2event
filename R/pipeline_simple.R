#' Flexible ml4time2event pipeline
#'
#' Provides a lightweight orchestration layer around `ml4t2e_fit()` with optional
#' preprocessing (via `recipes`) and resampling (via `rsample`). The pipeline
#' stores the fitted model, preprocessing artefacts, and resampling summaries so
#' that predictions and evaluation can be reproduced consistently.
#'
#' @keywords internal
#' @noRd
T2EPipeline <- R6::R6Class(
  classname = "T2EPipeline",
  public = list(
    outcome = NULL,
    models = NULL,
    ensemble = NULL,
    controls = NULL,
    recipe = NULL,
    resampling = NULL,
    metrics = NULL,
    fit_object = NULL,
    task = NULL,
    prepped_recipe = NULL,
    features = NULL,
    training_metrics = NULL,
    resampling_results = NULL,

    initialize = function(outcome,
                          models = c("cox", "gbm"),
                          ensemble = "auto",
                          controls = list(),
                          recipe = NULL,
                          resampling = NULL,
                          metrics = NULL) {
      self$outcome <- .pipeline_validate_outcome(outcome)
      self$models <- unique(models)
      self$ensemble <- ensemble
      self$controls <- controls %||% list()
      self$recipe <- recipe
      self$resampling <- resampling
      self$metrics <- metrics
    },

    fit = function(data) {
      if (inherits(data, "t2e_task")) {
        processed <- .pipeline_prepare_from_task(data, self$recipe)
        processed_task <- processed$task
        self$prepped_recipe <- processed$prepped_recipe
        self$features <- processed_task$features
        self$task <- processed_task
      } else {
        prepared <- .pipeline_prepare_data(
          data = data,
          outcome = self$outcome,
          recipe = self$recipe
        )
        self$prepped_recipe <- prepared$prepped_recipe
        self$features <- prepared$features
        self$task <- prepared$task
      }

      self$fit_object <- ml4t2e_fit(
        task = self$task,
        models = self$models,
        ensemble = self$ensemble,
        controls = self$controls
      )

      metric_set <- self$metrics
      if (is.null(metric_set)) {
        metric_set <- .pipeline_default_metrics(self$task)
      }
      if (length(metric_set) > 0) {
        self$training_metrics <- ml4t2e_evaluate(self$fit_object, metrics = metric_set)
      } else {
        self$training_metrics <- NULL
      }

      if (!is.null(self$resampling)) {
        self$resampling_results <- .pipeline_resample(
          data = if (inherits(data, "t2e_task")) data$data else data,
          outcome = self$outcome,
          recipe = self$recipe,
          models = self$models,
          ensemble = self$ensemble,
          controls = self$controls,
          resampling = self$resampling,
          metrics = metric_set
        )
      }

      invisible(self)
    },

    predict = function(newdata = NULL, ...) {
      if (is.null(self$fit_object)) {
        rlang::abort("Pipeline must be fitted before calling `predict()`.")
      }
      if (is.null(newdata)) {
        return(predict(self$fit_object, ...))
      }
      processed <- .pipeline_process_new_data(
        pipeline = self,
        newdata = newdata
      )
      predict(self$fit_object, newdata = processed, ...)
    },

    evaluate = function(newdata = NULL,
                        metrics = NULL,
                        include = "all") {
      if (is.null(self$fit_object)) {
        rlang::abort("Pipeline must be fitted before calling `evaluate()`.")
      }
      metric_set <- metrics %||% self$metrics
      if (is.null(metric_set)) {
        metric_set <- .pipeline_default_metrics(self$task)
      }
      if (length(metric_set) == 0) {
        if (is.null(newdata)) {
          return(NULL)
        }
      }
      if (is.null(newdata)) {
        return(ml4t2e_evaluate(self$fit_object, metrics = metric_set, include = include))
      }
      processed <- .pipeline_process_new_data(
        pipeline = self,
        newdata = newdata,
        require_outcomes = TRUE
      )
      eval_task <- .pipeline_build_task(
        data = processed,
        outcome = self$outcome,
        features = self$features
      )
      preds <- predict(self$fit_object, newdata = processed, include = include)
      ml4t2e_evaluate(preds, task = eval_task, metrics = metric_set, include = include)
    },

    summary = function() {
      list(
        outcome = self$outcome,
        models = self$models,
        ensemble = self$ensemble,
        features = self$features,
        training_metrics = self$training_metrics,
        resampling_results = self$resampling_results
      )
    },

    print = function(...) {
      cat("ml4time2event Pipeline\n")
      cat("======================\n")
      cat("Outcome type   :", self$outcome$type, "\n")
      cat("Models         :", paste(self$models, collapse = ", "), "\n")
      ens_strategy <- if (!is.null(self$fit_object)) {
        self$fit_object$ensemble_strategy %||% "none"
      } else {
        "none"
      }
      cat("Ensemble       :", ens_strategy, "\n")
      if (!is.null(self$features)) {
        cat("Predictors     :", length(self$features), "\n")
      }
      if (!is.null(self$training_metrics)) {
        cat("Training metrics:\n")
        print(self$training_metrics)
      }
      if (!is.null(self$resampling_results)) {
        cat("Resampling metrics:\n")
        print(self$resampling_results)
      }
      invisible(self)
    }
  )
)

#' Create a flexible pipeline builder
#'
#' @param outcome Named list describing the analysis outcome:
#'   - Survival: `list(type = "survival", time = "time", event = "status", id = NULL, features = NULL, time_units = NULL)`
#'   - Competing risk: `list(type = "competing_risk", time = "time", status = "status", cause = "cause", id = NULL, features = NULL, time_units = NULL)`
#' @param models Character vector of model engines to fit.
#' @param ensemble Ensemble strategy passed to `ml4t2e_fit()`.
#' @param controls Optional named list of per-engine control settings.
#' @param recipe Optional `recipes::recipe` (or function returning one) used for preprocessing.
#' @param resampling Optional `rsample` resampling object (e.g., `vfold_cv`).
#' @param metrics Optional character vector of metric names for evaluation.
#'
#' @return An unfitted pipeline object. Call `$fit(data)` to train.
#' @export
ml4t2e_pipeline <- function(outcome,
                            models = c("cox", "gbm"),
                            ensemble = "auto",
                            controls = list(),
                            recipe = NULL,
                            resampling = NULL,
                            metrics = NULL) {
  pipeline <- T2EPipeline$new(
    outcome = outcome,
    models = models,
    ensemble = ensemble,
    controls = controls,
    recipe = recipe,
    resampling = resampling,
    metrics = metrics
  )
  class(pipeline) <- unique(c("ml4t2e_pipeline", class(pipeline)))
  pipeline
}

# -------------------------------------------------------------------------
# Helper functions (internal)
# -------------------------------------------------------------------------

.pipeline_validate_outcome <- function(outcome) {
  if (!is.list(outcome)) {
    rlang::abort("`outcome` must be a list describing the analysis outcome.")
  }
  outcome$type <- tolower(outcome$type %||% rlang::abort("`outcome$type` must be provided (survival or competing_risk)."))
  if (outcome$type == "survival") {
    required <- c("time", "event")
  } else if (outcome$type == "competing_risk") {
    required <- c("time", "status", "cause")
  } else {
    rlang::abort("`outcome$type` must be 'survival' or 'competing_risk'.")
  }
  missing <- setdiff(required, names(outcome))
  if (length(missing) > 0) {
    rlang::abort(paste0("Outcome definition missing fields: ", paste(missing, collapse = ", ")))
  }
  outcome
}

.pipeline_prepare_from_task <- function(task, recipe) {
  data <- task$data
  outcome <- list(
    type = attr(task, "task_type"),
    time = task$time_col,
    event = if (attr(task, "task_type") == "survival") task$event_col else NULL,
    status = if (attr(task, "task_type") == "competing_risk") task$metadata$status_col else NULL,
    cause = if (attr(task, "task_type") == "competing_risk") task$cause_col else NULL,
    id = task$id_col,
    features = task$features,
    time_units = task$time_units
  )
  prepared <- .pipeline_prepare_data(
    data = data,
    outcome = outcome,
    recipe = recipe
  )
  list(
    task = prepared$task,
    prepped_recipe = prepared$prepped_recipe
  )
}

.pipeline_prepare_data <- function(data, outcome, recipe) {
  if (!is.data.frame(data)) {
    rlang::abort("Training data must be a data.frame.")
  }
  required <- .pipeline_required_columns(outcome)
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0) {
    rlang::abort(paste0("Training data missing columns: ", paste(missing, collapse = ", ")))
  }

  if (!is.null(recipe)) {
    if (!requireNamespace("recipes", quietly = TRUE)) {
      rlang::abort("Package 'recipes' must be installed to use preprocessing.")
    }
    recipe_obj <- recipe
    if (is.function(recipe_obj)) {
      recipe_obj <- recipe_obj(data)
    }
    if (!inherits(recipe_obj, "recipe")) {
      rlang::abort("`recipe` must be a recipes::recipe or a function returning one.")
    }
    prepped <- recipes::prep(recipe_obj, training = data, retain = TRUE)
    processed <- recipes::bake(prepped, new_data = NULL)
    if (!is.null(outcome$cause) && !outcome$cause %in% colnames(processed)) {
      processed[[outcome$cause]] <- data[[outcome$cause]]
    }
  } else {
    prepped <- NULL
    processed <- data
  }

  features <- outcome$features
  if (is.null(features)) {
    features <- .pipeline_default_features(processed, outcome)
  } else {
    features <- intersect(features, colnames(processed))
  }

  task <- .pipeline_build_task(
    data = processed,
    outcome = outcome,
    features = features
  )

  list(
    task = task,
    prepped_recipe = prepped,
    features = features
  )
}

.pipeline_build_task <- function(data, outcome, features) {
  if (outcome$type == "survival") {
    ml4t2e_task_surv(
      data = data,
      time = outcome$time,
      event = outcome$event,
      features = features,
      id = outcome$id,
      time_units = outcome$time_units
    )
  } else {
    ml4t2e_task_cr(
      data = data,
      time = outcome$time,
      status = outcome$status,
      cause = outcome$cause,
      features = features,
      id = outcome$id,
      time_units = outcome$time_units
    )
  }
}

.pipeline_default_features <- function(data, outcome) {
  outcome_cols <- .pipeline_required_columns(outcome)
  setdiff(colnames(data), outcome_cols)
}

.pipeline_required_columns <- function(outcome) {
  base <- c(outcome$time, outcome$id)
  if (outcome$type == "survival") {
    c(base, outcome$event)
  } else {
    c(base, outcome$status, outcome$cause)
  }
}

.pipeline_process_new_data <- function(pipeline,
                                       newdata,
                                       require_outcomes = FALSE) {
  if (!is.data.frame(newdata)) {
    rlang::abort("New data must be a data.frame.")
  }
  required <- .pipeline_required_columns(pipeline$outcome)
  if (!require_outcomes) {
    missing <- setdiff(required, colnames(newdata))
    if (length(missing) > 0) {
      for (col in missing) {
        newdata[[col]] <- NA
      }
    }
  } else {
    missing <- setdiff(required, colnames(newdata))
    if (length(missing) > 0) {
      rlang::abort(paste0("Evaluation data missing columns: ", paste(missing, collapse = ", ")))
    }
  }

  if (!is.null(pipeline$prepped_recipe)) {
    recipes::bake(pipeline$prepped_recipe, new_data = newdata)
  } else {
    newdata
  }
}

.pipeline_default_metrics <- function(task) {
  if (attr(task, "task_type") == "survival") {
    c("c_index", "ibs")
  } else {
    c("c_index", "ibs")
  }
}

.pipeline_resample <- function(data,
                               outcome,
                               recipe,
                               models,
                               ensemble,
                               controls,
                               resampling,
                               metrics) {
  if (!requireNamespace("rsample", quietly = TRUE)) {
    rlang::abort("Package 'rsample' must be installed to use resampling.")
  }
  splits <- resampling$splits %||% resampling
  results <- list()
  for (i in seq_along(splits)) {
    split <- splits[[i]]
    analysis_data <- rsample::analysis(split)
    assessment_data <- rsample::assessment(split)

    train_prepared <- .pipeline_prepare_data(
      data = analysis_data,
      outcome = outcome,
      recipe = recipe
    )
    train_task <- train_prepared$task
    prepped <- train_prepared$prepped_recipe
    features <- train_prepared$features

    assessment_processed <- if (!is.null(prepped)) {
      recipes::bake(prepped, new_data = assessment_data)
    } else {
      assessment_data
    }
    assess_task <- .pipeline_build_task(
      data = assessment_processed,
      outcome = outcome,
      features = features
    )

    fit_obj <- ml4t2e_fit(
      task = train_task,
      models = models,
      ensemble = ensemble,
      controls = controls
    )

    preds <- predict(fit_obj, newdata = assessment_processed)
    if (length(metrics) > 0) {
      metrics_df <- ml4t2e_evaluate(preds, task = assess_task, metrics = metrics)
      metrics_df$split <- i
      results[[length(results) + 1]] <- metrics_df
    }
  }

  if (length(results) == 0) {
    return(NULL)
  }

  dplyr::bind_rows(results)
}
