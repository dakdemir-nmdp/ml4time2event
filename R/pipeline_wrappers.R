#' @title Fit an ml4time2event modeling pipeline
#'
#' @description Convenience wrapper that combines preprocessing, model fitting,
#' persistence, and prediction for survival or competing-risks analyses.
#'
#' @param data Training data frame.
#' @param analysis_type Either `"survival"` or `"competing_risks"`.
#' @param timevar Name of time-to-event column.
#' @param eventvar Name of event indicator column (0/1 for survival, 0/1/2... for CR).
#' @param expvars Optional character vector of predictors. Defaults to all columns
#'   except outcomes and IDs.
#' @param idvars Optional character vector of identifier columns that should pass
#'   through preprocessing unchanged.
#' @param recipe_fn Preprocessing function accepting `model_recipe` as first argument
#'   (defaults to [minimal_data_recipe()]).
#' @param recipe_args Named list of additional arguments passed to `recipe_fn`.
#' @param prediction_times Optional numeric vector specifying the default time grid
#'   for downstream predictions. When `NULL`, the grid spans the observed follow-up.
#' @param models Optional character vector passed to [RunSurvModels()] or
#'   [RunCRModels()]. When `NULL`, sensible defaults are used.
#' @param ntreeRF Number of trees for random forest models (forwarded to Run* helpers).
#' @param varsel Logical flag forwarded to [RunCRModels()] indicating variable selection.
#' @param include_rf Logical; if `FALSE`, skips fitting the baseline random forest
#'   models within [RunSurvModels()] and [RunCRModels()].
#' @param ... Additional arguments forwarded to [RunSurvModels()] or [RunCRModels()].
#'
#' @return An object of class `ml4t2e_pipeline` containing the fitted preprocessing
#'   recipe, ensemble model, and metadata.
#' @export
ml4t2e_fit_pipeline <- function(
    data,
    analysis_type = c("survival", "competing_risks"),
    timevar,
    eventvar,
    expvars = NULL,
    idvars = character(0),
    recipe_fn = minimal_data_recipe,
    recipe_args = list(pmiss = 0.25, pother = 0.05, dummy = FALSE),
    prediction_times = NULL,
    models = NULL,
    ntreeRF = 300,
    varsel = FALSE,
    include_rf = TRUE,
    ...
) {
  analysis_type <- match.arg(analysis_type)

  if (!requireNamespace("recipes", quietly = TRUE)) {
    stop("Package 'recipes' must be installed to use ml4t2e_fit_pipeline().", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  if (missing(timevar) || !timevar %in% colnames(data)) {
    stop("'timevar' must be present in 'data'.")
  }
  if (missing(eventvar) || !eventvar %in% colnames(data)) {
    stop("'eventvar' must be present in 'data'.")
  }
  if (!is.function(recipe_fn)) {
    stop("'recipe_fn' must be a function (e.g., minimal_data_recipe).")
  }
  if (is.null(expvars)) {
    expvars <- setdiff(colnames(data), c(timevar, eventvar, idvars))
  }
  if (length(expvars) == 0) {
    stop("No predictor variables identified. Provide 'expvars' or check input data.")
  }

  base_recipe <- t2emodel_data_recipe_init(
    timevar = timevar,
    eventvar = eventvar,
    expvar = expvars,
    idvars = idvars,
    traindata = data
  )

  recipe_call <- c(list(model_recipe = base_recipe), recipe_args)
  processed_recipe <- do.call(recipe_fn, recipe_call)
  prepped_recipe <- prep_data_recipe(processed_recipe, training = data)

  baked_train <- bake_data_recipe(prepped_recipe, data = data)
  processed_expvars <- setdiff(colnames(baked_train), c(timevar, eventvar, idvars))

  if (length(processed_expvars) == 0) {
    stop("Preprocessing removed all predictors. Adjust recipe settings.")
  }

  default_time_grid <- .ml4t2e_compute_time_grid(
    baked_train[[timevar]],
    prediction_times
  )

  model_defaults <- .ml4t2e_default_models(analysis_type, models)

  if (identical(analysis_type, "survival")) {
    fitted_models <- RunSurvModels(
      datatrain = baked_train,
      ExpVars = processed_expvars,
      timevar = timevar,
      eventvar = eventvar,
      models = model_defaults,
      ntreeRF = ntreeRF,
      run_rf = include_rf,
      ...
    )
    ensemble_object <- SurvEnsemble(fitted_models)
  } else {
    fitted_models <- RunCRModels(
      datatrain = baked_train,
      ExpVars = processed_expvars,
      timevar = timevar,
      eventvar = eventvar,
      models = model_defaults,
      ntreeRF = ntreeRF,
      varsel = varsel,
      run_rf = include_rf,
      ...
    )
    ensemble_object <- CREnsemble(fitted_models)
  }

  pipeline <- list(
    analysis_type = analysis_type,
    timevar = timevar,
    eventvar = eventvar,
    idvars = idvars,
    original_expvars = expvars,
    processed_expvars = processed_expvars,
    recipe = prepped_recipe,
    recipe_fn = .ml4t2e_fn_name(recipe_fn, substitute(recipe_fn)),
    recipe_args = recipe_args,
    model = ensemble_object,
    prediction_grid = default_time_grid,
    training_summary = .ml4t2e_training_summary(baked_train, timevar, eventvar, analysis_type),
    include_rf = include_rf
  )

  class(pipeline) <- c("ml4t2e_pipeline", "list")
  pipeline
}

#' @export
print.ml4t2e_pipeline <- function(x, ...) {
  cat("ml4time2event Pipeline\n")
  cat("======================\n")
  cat("Analysis type :", x$analysis_type, "\n")
  cat("Time variable :", x$timevar, "\n")
  cat("Event variable:", x$eventvar, "\n")
  cat(sprintf("Predictors    : %d processed (%d original)\n",
              length(x$processed_expvars), length(x$original_expvars)))
  cat("Recipe        :", x$recipe_fn, "\n")
  if (!is.null(x$training_summary$rows)) {
    cat("Observations  :", x$training_summary$rows, "\n")
  }
  invisible(x)
}

#' @title Save an ml4time2event pipeline
#' @param pipeline Object returned by [ml4t2e_fit_pipeline()].
#' @param file File path (an `.rds` extension is added when missing).
#' @param compress Compression setting passed to [saveRDS()].
#' @return Invisibly returns the file path.
#' @export
ml4t2e_save_pipeline <- function(pipeline, file, compress = TRUE) {
  if (!inherits(pipeline, "ml4t2e_pipeline")) {
    stop("'pipeline' must be created with ml4t2e_fit_pipeline().")
  }
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".rds")
  }
  pipeline$.__metadata__ <- list(
    saved_date = Sys.time(),
    r_version = R.version.string,
    package_version = utils::packageVersion("ml4time2event")
  )
  saveRDS(pipeline, file = file, compress = compress)
  invisible(file)
}

#' @title Load an ml4time2event pipeline
#' @param file Path to an `.rds` file created by [ml4t2e_save_pipeline()].
#' @return A reconstructed `ml4t2e_pipeline` object.
#' @export
ml4t2e_load_pipeline <- function(file) {
  if (!file.exists(file)) {
    stop("Pipeline file not found: ", file)
  }
  pipeline <- readRDS(file)
  if (!inherits(pipeline, "ml4t2e_pipeline")) {
    stop("File does not contain an ml4t2e_pipeline object.")
  }
  pipeline
}

#' @export
predict.ml4t2e_pipeline <- function(object, newdata, new_times = NULL,
                                    ensemble_method = "average", ...) {
  if (!inherits(object, "ml4t2e_pipeline")) {
    stop("object must be created with ml4t2e_fit_pipeline().")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data.frame.")
  }

  completed_newdata <- .ml4t2e_prepare_newdata(
    newdata,
    required = unique(c(object$timevar, object$eventvar,
                        object$original_expvars, object$idvars))
  )

  baked_new <- bake_data_recipe(object$recipe, data = completed_newdata)

  prediction_grid <- .ml4t2e_compute_time_grid(
    baked_new[[object$timevar]],
    new_times,
    fallback = object$prediction_grid
  )

  if (identical(object$analysis_type, "survival")) {
    preds <- PredictSurvModels(
      models = object$model,
      newdata = baked_new,
      new_times = prediction_grid,
      ensemble_method = ensemble_method,
      ...
    )
  } else {
    preds <- PredictCRModels(
      models = object$model,
      newdata = baked_new,
      new_times = prediction_grid,
      ensemble_method = ensemble_method,
      ...
    )
  }

  list(
    predictions = preds,
    baked_data = baked_new,
    times = prediction_grid,
    analysis_type = object$analysis_type
  )
}

.ml4t2e_compute_time_grid <- function(observed_times, override = NULL, fallback = NULL) {
  if (!is.null(override)) {
    override <- sort(unique(as.numeric(override)))
    override <- override[!is.na(override)]
    if (length(override) > 0) {
      return(unique(c(0, override)))
    }
  }

  if (!is.null(fallback)) {
    return(fallback)
  }

  observed_times <- as.numeric(observed_times)
  observed_times <- observed_times[is.finite(observed_times)]
  if (length(observed_times) == 0) {
    return(seq(0, 1, length.out = 50))
  }
  max_time <- max(observed_times, na.rm = TRUE)
  seq(0, max_time, length.out = 50)
}

.ml4t2e_training_summary <- function(data, timevar, eventvar, analysis_type) {
  out <- list(rows = nrow(data))
  if (analysis_type == "survival") {
    events <- data[[eventvar]]
    out$event_rate <- mean(events == 1, na.rm = TRUE)
    out$time_range <- range(data[[timevar]], na.rm = TRUE)
  } else {
    out$event_counts <- as.list(table(data[[eventvar]], useNA = "no"))
    out$time_range <- range(data[[timevar]], na.rm = TRUE)
  }
  out
}

.ml4t2e_default_models <- function(analysis_type, models) {
  if (!is.null(models)) {
    return(models)
  }
  if (identical(analysis_type, "survival")) {
    return(c("glmnet", "coxph", "xgboost", "gam"))
  }
  c("FG", "cox", "xgboost")
}

.ml4t2e_prepare_newdata <- function(newdata, required) {
  missing_cols <- setdiff(required, colnames(newdata))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      newdata[[col]] <- NA
    }
  }
  newdata[required]
}

.ml4t2e_fn_name <- function(fn, expr = NULL) {
  if (is.character(fn)) {
    return(fn)
  }
  if (!is.null(expr)) {
    expr_char <- tryCatch(deparse(expr), error = function(e) character())
    if (length(expr_char) > 0) {
      candidate <- trimws(expr_char[1])
      if (nzchar(candidate) && !candidate %in% c("fn")) {
        return(candidate)
      }
    }
  }
  fn_name <- tryCatch({
    attr(fn, "name")
  }, error = function(e) NULL)
  if (!is.null(fn_name) && nzchar(fn_name)) {
    return(fn_name)
  }
  "function"
}
