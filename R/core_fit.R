#' @import stats
#' @import utils
#' @import grDevices
#' @importFrom stats AIC as.formula coef complete.cases cor median model.matrix optim predict reorder runif sd setNames terms var
#' @importFrom utils capture.output head modifyList tail
# Suppress R CMD check NOTEs about variables used in NSE contexts and data.table operators
utils::globalVariables(c(
  "id", "time", "cif", ".data", "self", "private", "super",
  "age", "sex", "status", "ftime", "ph.ecog", "ph.karno",
  "pat.karno", "meal.cal", "wt.loss", "inst", "phase", "d",
  "PredictSurvModels", "PredictCRModels", ":="
))
#'
#' Fit time-to-event models via a unified interface
#'
#' @param task An object created with `ml4t2e_task_surv()` or `ml4t2e_task_cr()`.
#' @param models Character vector of engine names to fit. Use
#'   `ml4t2e_list_models()` to discover available options.
#'   Common engines include:
#'   - Survival: `"cox"`, `"random_forest"`, `"gbm"`, `"xgboost"`, `"shallownn"`, `"ttah"`, `"survreg"`.
#'   - Competing Risk: `"cox"`, `"fine_gray"`, `"cr_random_forest"`, `"cr_xgboost"`, `"cr_shallownn"`, `"cr_ttah"`.
#' @param ensemble Character or logical; controls ensembling strategy:
#'   - `FALSE` or `"none"`: No ensembling; returns individual model predictions only.
#'   - `"simple"`: Unweighted average of predictions across models.
#'   - `"auto"`, `TRUE`: Alias for `"simple"`.
#'   - `"stack"`: Uses blending (data splitting) to learn ensemble weights via Super Learner optimization.
#'
#'   Default: `FALSE`.
#'
#'   **Note**: When `ensemble = "stack"`, a portion of the data (default 20%) is held out
#'   for meta-learner training.
#' @param recipe Optional `recipes` recipe for preprocessing. Not yet supported
#'   in this refactor; placeholder for future work.
#' @param resampling Optional resampling specification from `rsample`.
#' @param seed Optional integer seed for reproducibility.
#' @param keep_data Logical; if `TRUE`, retains the training data in the returned
#'   object. Default is `FALSE` to save memory.
#'   **Warning**: If `FALSE`, you must provide `newdata` to `predict()` and
#'   `task` to `ml4t2e_evaluate()`.
#' @param controls Named list of engine-specific parameter lists. Use entries
#'   like `controls = list(cox = list(penalized = TRUE))`. Can also include
#'   global controls like `time_limit` (in seconds) to stop training after a
#'   certain duration.
#' @param conformal_calibration Optional numeric value between 0 and 1. If provided,
#'   specifies the fraction of data to use for conformal calibration (and stacking).
#'
#' @return An object of class `t2e_fit`.
#' @examples
#' \dontrun{
#' lung_df <- survival::lung
#' lung_df$status <- lung_df$status - 1 # 0=censored, 1=event
#' task <- ml4t2e_task_surv(lung_df, time = "time", event = "status", features = c("age", "sex"))
#' fit <- ml4t2e_fit(task, models = c("cox", "random_forest"), keep_data = TRUE)
#' predict(fit, newdata = lung_df[1:5, ], times = c(100, 200, 300))
#' }
#' @export
ml4t2e_fit <- function(task,
                       models = "cox",
                       ensemble = c(FALSE, TRUE, "auto", "stack", "simple"),
                       recipe = NULL,
                       resampling = NULL,
                       seed = NULL,
                       keep_data = TRUE,
                       controls = list(),
                       conformal_calibration = NULL) {
  # Validate task structure
  .validate_task(task, context = "ml4t2e_fit()")

  # Validate seed if provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      rlang::abort("`seed` must be a single numeric value.")
    }
    seed <- as.integer(seed)
    set.seed(seed)
  }

  ensemble_choice <- match.arg(as.character(ensemble[[1]]), c("FALSE", "TRUE", "auto", "stack", "simple", "none"))
  ensemble_mode <- ensemble_choice
  if (ensemble_choice %in% c("TRUE", "auto")) {
    ensemble_mode <- "simple"
  }
  ensemble_mode <- tolower(ensemble_mode)
  ensemble_enabled <- !ensemble_choice %in% c("FALSE", "none")

  # Validate ensemble strategy
  if (ensemble_enabled) {
    .validate_ensemble_strategy(ensemble_mode, context = "ml4t2e_fit()")
  }

  outcome_type <- attr(task, "task_type")
  models <- unique(models)

  # Validate models against registry
  .validate_models(models, outcome_type, context = "ml4t2e_fit()")

  train_task <- task
  stack_data <- NULL
  conf_data <- NULL

  # Data splitting for stacking and/or conformal calibration
  is_stacking <- (ensemble_mode == "stack")
  is_conformal <- (!is.null(conformal_calibration))

  if (is_stacking || is_conformal) {
    # Validate feasibility first
    .validate_split_feasibility(
      n_total = nrow(task$data),
      stacking_requested = is_stacking,
      conformal_requested = is_conformal,
      conformal_ratio = conformal_calibration %||% 0.2
    )

    # Perform split
    split_result <- .split_data_for_training(
      task = task,
      is_stacking = is_stacking,
      is_conformal = is_conformal,
      conformal_ratio = conformal_calibration %||% 0.2
    )

    # Unpack results
    train_task <- task
    train_task$data <- split_result$train_data
    stack_data <- split_result$stack_data
    conf_data <- split_result$conf_data

    # Inform user
    rlang::inform(split_result$split_info)
  }

  time_grid <- .resolve_time_grid(train_task, controls)

  # Time limit check setup
  time_limit <- controls$time_limit
  start_time <- Sys.time()

  fitted_models <- list()
  for (engine in models) {
    if (!is.null(time_limit) && as.numeric(difftime(Sys.time(), start_time, units = "secs")) > time_limit) {
      rlang::warn("Time limit exceeded. Skipping remaining models.")
      break
    }

    spec <- controls[[engine]]
    if (is.null(spec)) {
      spec <- list()
    }
    instance <- .instantiate_model(engine = engine, outcome = outcome_type, spec = spec)
    trained <- instance$fit(task = train_task, time_grid = time_grid)
    if (inherits(trained, "R6")) {
      fitted_models[[engine]] <- trained
    } else {
      fitted_models[[engine]] <- instance
    }
  }

  ensemble_controls <- controls$ensemble
  ensemble_obj <- NULL

  if (ensemble_enabled) {
    sl_data <- list(preds = NULL, actual = NULL)

    # Generate Stacking Data if needed
    if (ensemble_mode == "stack" && !is.null(stack_data)) {
      sl_data <- prepare_stacking_data(
        task = task,
        fitted_models = fitted_models,
        time_grid = time_grid,
        outcome_type = outcome_type,
        cal_data = stack_data
      )
    }

    ensemble_obj <- .build_ensembler(
      models = fitted_models,
      model_names = names(fitted_models),
      outcome_type = outcome_type,
      strategy = ensemble_mode,
      time_grid = time_grid,
      task = train_task,
      controls = ensemble_controls,
      sl_predictions = sl_data$preds,
      sl_actual = sl_data$actual,
      sl_weights = sl_data$weights
    )
  }

  # Prepare task for storage: strip data if requested
  stored_task <- train_task
  if (!isTRUE(keep_data)) {
    stored_task$data <- NULL
  }

  fit <- list(
    task = stored_task,
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

  # Conformal calibration if requested
  if (!is.null(conformal_calibration)) {
    if (is.null(conf_data)) {
      rlang::warn("Calibration data invalid; skipping conformal calibration.")
    } else {
      fit <- ml4t2e_calibrate(fit, conf_data)
    }
  }

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
                             controls = list(),
                             sl_predictions = NULL,
                             sl_actual = NULL,
                             sl_weights = NULL) {
  strategy <- tolower(strategy %||% "simple")
  if (identical(strategy, "false")) {
    return(NULL)
  }

  if (identical(outcome_type, "survival")) {
    return(SurvivalEnsembler$new(
      models = models,
      model_names = model_names,
      outcome_type = outcome_type,
      strategy = strategy,
      time_grid = time_grid,
      controls = controls,
      sl_predictions = sl_predictions,
      sl_actual = sl_actual
    ))
  }

  if (identical(outcome_type, "competing_risk")) {
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
      controls = controls,
      sl_predictions = sl_predictions,
      sl_actual = sl_actual
    ))
  }

  NULL
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
