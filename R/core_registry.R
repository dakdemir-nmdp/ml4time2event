#' Internal registry for time-to-event model engines
#'
#' Engines register themselves at package load so that the public API can look up
#' constructors, metadata, and dependencies.
#'
#' @keywords internal
.t2e_model_registry <- new.env(parent = emptyenv())
.t2e_model_registry$models <- list()

.register_time_to_event_model <- function(engine,
                                          outcome,
                                          constructor,
                                          packages = character(),
                                          tags = character(),
                                          label = NULL) {
  outcome <- match.arg(outcome, c("survival", "competing_risk"))
  if (!is.function(constructor)) {
    rlang::abort("`constructor` must be a function that returns an R6 model instance.")
  }
  key <- paste0(outcome, "::", engine)
  label_value <- label
  if (is.null(label_value)) {
    label_value <- engine
  }
  .t2e_model_registry$models[[key]] <- list(
    engine = engine,
    outcome = outcome,
    constructor = constructor,
    packages = packages,
    tags = tags,
    label = label_value
  )
  invisible(TRUE)
}

.t2e_model_registry_get <- function(engine, outcome) {
  key <- paste0(outcome, "::", engine)
  model <- .t2e_model_registry$models[[key]]
  if (is.null(model)) {
    rlang::abort(glue::glue("No registered {outcome} engine named '{engine}'."))
  }
  model
}

#' List available model engines
#'
#' @param outcome Filter by outcome type. `"survival"` and `"competing_risk"` are
#'   accepted. Use `"all"` to list everything.
#'
#' @return A data frame with columns `engine`, `outcome`, `label`, `packages`,
#'   and `tags`.
#' @export
ml4t2e_list_models <- function(outcome = c("survival", "competing_risk", "all")) {
  outcome <- match.arg(outcome)
  entries <- .t2e_model_registry$models
  if (length(entries) == 0) {
    return(data.frame(engine = character(), outcome = character(), label = character(),
                      packages = I(list()), tags = I(list()), stringsAsFactors = FALSE))
  }
  df <- lapply(entries, function(entry) {
    data.frame(
      engine = entry$engine,
      outcome = entry$outcome,
      label = entry$label,
      packages = I(list(entry$packages)),
      tags = I(list(entry$tags)),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df)
  if (outcome != "all") {
    df <- df[df$outcome == outcome, , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

.instantiate_model <- function(engine, outcome, spec = list()) {
  entry <- .t2e_model_registry_get(engine, outcome)
  needed <- entry$packages
  if (length(needed) > 0) {
    missing_pkgs <- needed[!vlapply(needed, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      rlang::abort(glue::glue(
        "The following packages are required for engine '{engine}': {paste(missing_pkgs, collapse = ', ')}"
      ))
    }
  }
  entry$constructor(spec = modifyList(spec, list(engine = entry$engine)))
}

vlapply <- function(x, FUN, ...) {
  vapply(x, FUN, logical(1), ...)
}
