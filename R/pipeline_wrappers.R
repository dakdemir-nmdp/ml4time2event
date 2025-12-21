#' Fit a time-to-event pipeline (backwards-compatible helper)
#'
#' This convenience wrapper bridges the legacy `ml4t2e_fit_pipeline()` workflow
#' onto the new `ml4t2e_pipeline()` + `T2EPipeline` infrastructure. It accepts
#' the original argument names but internally builds a modern pipeline,
#' returning an object with classes `c("ml4t2e_pipeline", "T2EPipeline", "R6")`.
#'
#' @param data Training data frame.
#' @param analysis_type `"survival"` or `"competing_risks"`.
#' @param timevar Name of the time column.
#' @param eventvar Name of the event/status column.
#' @param causevar Optional cause column for competing-risk tasks. When
#'   `NULL`, event codes from `eventvar` are reused for causes.
#' @param expvars Optional predictor set (defaults to all non-outcome columns).
#' @param idvars Optional identifier column (first entry is used).
#' @param recipe_fn Optional preprocessing function (defaults to
#'   `minimal_data_recipe()`); set to `NULL` to skip recipes.
#' @param recipe_args Named list of additional arguments for `recipe_fn`.
#' @param prediction_times Optional numeric time grid forwarded to `controls$times`.
#' @param models Character vector of engine names (legacy aliases are supported).
#' @param ensemble Ensemble strategy passed to `ml4t2e_fit()`.
#' @param controls Optional named list of per-engine controls.
#' @param resampling Optional `rsample` resampling object.
#' @param metrics Optional metric identifiers for `ml4t2e_evaluate()`.
#' @param include_rf,ntreeRF,varsel Legacy arguments retained for compatibility.
#'   They are ignored (with a warning when changed from defaults).
#' @param ... Unused; retained for compatibility.
#'
#' @return A fitted `T2EPipeline` object (also inheriting the
#'   `ml4t2e_pipeline` class for S3 helpers).
#' @keywords internal
#' @export
ml4t2e_fit_pipeline <- function(data,
                                analysis_type = c("survival", "competing_risks"),
                                timevar,
                                eventvar,
                                causevar = NULL,
                                expvars = NULL,
                                idvars = character(0),
                                recipe_fn = minimal_data_recipe,
                                recipe_args = list(pmiss = 0.25, pother = 0.05, dummy = FALSE),
                                prediction_times = NULL,
                                models = NULL,
                                ensemble = "auto",
                                controls = list(),
                                resampling = NULL,
                                metrics = NULL,
                                include_rf = TRUE,
                                ntreeRF = NULL,
                                varsel = NULL,
                                ...) {
  analysis_type <- match.arg(analysis_type)
  data <- as.data.frame(data)

  # Warn about deprecated knobs
  if (!isTRUE(include_rf)) {
    rlang::warn("`include_rf` is deprecated; specify engines explicitly via `models`.")
  }
  if (!is.null(ntreeRF)) {
    rlang::warn("`ntreeRF` is deprecated and ignored; use `controls$random_forest` if needed.")
  }
  if (!is.null(varsel)) {
    rlang::warn("`varsel` is deprecated and ignored.")
  }
  if (length(idvars) > 1) {
    rlang::warn("Only one id column is supported; using the first entry in `idvars`.")
    idvars <- idvars[1]
  }

  required_cols <- c(timevar, eventvar)
  if (!all(required_cols %in% colnames(data))) {
    missing <- setdiff(required_cols, colnames(data))
    rlang::abort(paste0("Training data missing required columns: ", paste(missing, collapse = ", ")))
  }

  outcome_type <- if (identical(analysis_type, "survival")) "survival" else "competing_risk"
  data_prepped <- data
  id_col <- if (length(idvars) == 1) idvars else NULL
  cause_col <- NULL

  if (identical(outcome_type, "competing_risk")) {
    status_raw <- data_prepped[[eventvar]]
    if (is.null(causevar)) {
      cause_col <- ".ml4t2e_cause"
      data_prepped[[cause_col]] <- status_raw
    } else {
      if (!causevar %in% colnames(data_prepped)) {
        rlang::abort(paste0("Cause column '", causevar, "' not found in data."))
      }
      cause_col <- causevar
    }
    data_prepped[[eventvar]] <- ifelse(as.numeric(status_raw) == 0, 0L, 1L)
    data_prepped[[cause_col]][status_raw == 0] <- NA
  }

  drop_cols <- unique(c(timevar, eventvar, id_col, cause_col))
  if (is.null(expvars)) {
    expvars <- setdiff(colnames(data_prepped), drop_cols)
  }
  if (length(expvars) == 0) {
    rlang::abort("No predictor variables identified. Provide `expvars` or adjust input data.")
  }

  outcome <- list(
    type = outcome_type,
    time = timevar,
    id = id_col,
    features = expvars
  )
  if (identical(outcome_type, "survival")) {
    outcome$event <- eventvar
  } else {
    outcome$status <- eventvar
    outcome$cause <- cause_col
  }

  recipe_obj <- NULL
  if (!is.null(recipe_fn)) {
    if (!requireNamespace("recipes", quietly = TRUE)) {
      rlang::abort("Package 'recipes' must be installed to use preprocessing.")
    }
    base_recipe <- t2emodel_data_recipe_init(
      timevar = timevar,
      eventvar = eventvar,
      expvar = expvars,
      idvars = id_col %||% character(0),
      traindata = data_prepped
    )
    recipe_obj <- do.call(recipe_fn, c(list(model_recipe = base_recipe), recipe_args))
  }

  translated_models <- .ml4t2e_translate_models(models, outcome_type)
  if (!is.null(prediction_times)) {
    controls$times <- prediction_times
  }

  pipeline <- ml4t2e_pipeline(
    outcome = outcome,
    models = translated_models,
    ensemble = ensemble,
    controls = controls,
    recipe = recipe_obj,
    resampling = resampling,
    metrics = metrics
  )
  pipeline$fit(data_prepped)
  class(pipeline) <- unique(c("ml4t2e_pipeline", class(pipeline)))
  pipeline
}

#' @export
print.ml4t2e_pipeline <- function(x, ...) {
  if (inherits(x, "T2EPipeline")) {
    x$print(...)
  } else {
    NextMethod()
  }
}

#' @export
predict.ml4t2e_pipeline <- function(object, ...) {
  if (!inherits(object, "T2EPipeline")) {
    rlang::abort("`object` must be created with `ml4t2e_pipeline()` or `ml4t2e_fit_pipeline()`.")
  }
  object$predict(...)
}

#' @keywords internal
#' @noRd
ml4t2e_save_pipeline <- function(pipeline, file, compress = TRUE) {
  if (!inherits(pipeline, "T2EPipeline")) {
    rlang::abort("'pipeline' must be created with `ml4t2e_pipeline()` or `ml4t2e_fit_pipeline()`.")
  }
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".rds")
  }
  saveRDS(pipeline, file = file, compress = compress)
  invisible(file)
}

#' @keywords internal
#' @noRd
ml4t2e_load_pipeline <- function(file) {
  if (!file.exists(file)) {
    rlang::abort("Pipeline file not found: ", file)
  }
  pipeline <- readRDS(file)
  if (!inherits(pipeline, "T2EPipeline")) {
    rlang::abort("File does not contain a saved `ml4t2e_pipeline` object.")
  }
  class(pipeline) <- unique(c("ml4t2e_pipeline", class(pipeline)))
  pipeline
}

.ml4t2e_translate_models <- function(models, outcome_type) {
  if (is.null(models)) {
    return(if (identical(outcome_type, "survival")) c("cox", "random_forest") else c("cox", "fine_gray"))
  }
  map <- c(
    coxph = "cox",
    cox = "cox",
    fg = "fine_gray",
    finegray = "fine_gray",
    fine_gray = "fine_gray",
    randomforest = "random_forest",
    random_forest = "random_forest",
    rf = "random_forest"
  )
  normalized <- tolower(models)
  translated <- map[normalized]
  translated[is.na(translated)] <- normalized[is.na(translated)]
  unique(translated)
}
