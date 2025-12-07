#' Base R6 classes for time-to-event models
#'
#' These classes define the contract that concrete survival and competing risk
#' engines need to satisfy. Each engine should inherit from one of the concrete
#' subclasses and implement the abstract methods documented below.
#'
#' @keywords internal
#' @noRd
TimeToEventModel <- R6::R6Class(
  classname = "TimeToEventModel",
  public = list(
    #' @field spec Engine specification (list of hyperparameters / metadata)
    spec = NULL,
    #' @field outcome_type Either "survival" or "competing_risk"
    outcome_type = NULL,
    #' @field fitted Logical flag indicating whether `fit()` has been called
    fitted = FALSE,

    #' @description
    #' Create a new model instance
    #' @param spec List containing engine configuration parameters
    initialize = function(spec = list()) {
      if (is.null(spec)) {
        spec <- list()
      }
      self$spec <- spec
      if (!is.list(self$spec)) {
        rlang::abort("`spec` must be a list.")
      }
      private$validate_spec()
    },

    #' @description
    #' Fit the underlying model on a task object.
    #' @param task An object created by `ml4t2e_task_surv()` or `ml4t2e_task_cr()`
    #' @param ... Currently ignored but kept for future extensibility
    fit = function(task, ...) {
      rlang::abort(glue::glue("`fit()` is not implemented for {class(self)[1]}"))
    },

    #' @description
    #' Predict survival probabilities over a time grid.
    #' Subclasses must return a tibble produced by `new_survival_prediction()`.
    predict_survival = function(newdata, times, ...) {
      rlang::abort(glue::glue(
        "`predict_survival()` is not implemented for {class(self)[1]}"
      ))
    },

    #' @description
    #' Predict risk scores (e.g., linear predictor or hazard ratios).
    predict_risk = function(newdata, times = NULL, ...) {
      rlang::abort(glue::glue(
        "`predict_risk()` is not implemented for {class(self)[1]}"
      ))
    },

    #' @description
    #' Predict cumulative incidence functions. Applicable to competing risk models.
    predict_cif = function(newdata, times, ...) {
      rlang::abort(glue::glue(
        "`predict_cif()` is not implemented for {class(self)[1]}"
      ))
    },

    #' @description
    #' Predict event time summaries (e.g., median survival). Optional.
    predict_time = function(newdata, q = 0.5, ...) {
      rlang::abort(glue::glue(
        "`predict_time()` is not implemented for {class(self)[1]}"
      ))
    },

    #' @description
    #' Return a list describing the model for registries and summaries.
    model_info = function() {
      engine <- self$spec$engine
      if (is.null(engine)) {
        engine <- class(self)[1]
      }
      package <- self$spec$package
      if (is.null(package)) {
        package <- NA_character_
      }
      list(engine = engine, outcome = self$outcome_type, package = package, mode = "time_to_event")
    },

    #' @description
    #' Return a character vector of packages required by the engine.
    required_packages = function() {
      character()
    }
  ),
  private = list(
    validate_spec = function() {
      # Intentionally left blank; subclasses can override
    }
  ),
  lock_objects = FALSE
)

#' @keywords internal
#' @noRd
SurvivalModel <- R6::R6Class(
  classname = "SurvivalModel",
  inherit = TimeToEventModel,
  public = list(
    initialize = function(spec = list()) {
      super$initialize(spec = modifyList(spec, list(outcome = "survival")))
      self$outcome_type <- "survival"
    },

    #' @description Ensure the supplied task is a survival task before fitting
    fit = function(task, ...) {
      private$check_task(task, expected = "survival")
      self$fitted <- TRUE
      invisible(self)
    }
  ),
  private = list(
    check_task = function(task, expected) {
      if (!inherits(task, "t2e_task")) {
        rlang::abort("`task` must be created with `ml4t2e_task_*()`.")
      }
      task_type <- attr(task, "task_type")
      if (!identical(task_type, expected)) {
        rlang::abort(glue::glue("Expected a {expected} task, got `{task_type}`."))
      }
    }
  )
)

#' @keywords internal
#' @noRd
CompetingRiskModel <- R6::R6Class(
  classname = "CompetingRiskModel",
  inherit = TimeToEventModel,
  public = list(
    initialize = function(spec = list()) {
      super$initialize(spec = modifyList(spec, list(outcome = "competing_risk")))
      self$outcome_type <- "competing_risk"
    },

    fit = function(task, ...) {
      private$check_task(task, expected = "competing_risk")
      self$fitted <- TRUE
      invisible(self)
    }
  ),
  private = list(
    check_task = function(task, expected) {
      if (!inherits(task, "t2e_task")) {
        rlang::abort("`task` must be created with `ml4t2e_task_*()`.")
      }
      task_type <- attr(task, "task_type")
      if (!identical(task_type, expected)) {
        rlang::abort(glue::glue("Expected a {expected} task, got `{task_type}`."))
      }
    }
  )
)
