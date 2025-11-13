#' Create a survival analysis task
#'
#' This constructor standardises the inputs needed to fit survival models and
#' ensures that downstream components can rely on a consistent structure.
#'
#' @param data A data frame containing the outcome and predictor variables.
#' @param time Name of the time-to-event column.
#' @param event Name of the event indicator (0/1) column.
#' @param features Optional character vector of predictor columns. Defaults to
#'   all non-outcome columns.
#' @param id Optional column containing subject identifiers. If `NULL`, a
#'   sequential identifier is created.
#' @param time_units Optional description of the time scale (e.g., "days").
#'
#' @return An object of class `t2e_task_surv` (inherits from `t2e_task`).
#' @export
ml4t2e_task_surv <- function(data,
                             time,
                             event,
                             features = NULL,
                             id = NULL,
                             time_units = NULL) {
  .task_constructors_validate_arguments(time, event, id)

  data_tbl <- dplyr::as_tibble(data)
  .task_constructors_check_columns(data_tbl, c(time, event, id))

  if (is.null(id)) {
    id <- ".task_id"
    data_tbl[[id]] <- seq_len(nrow(data_tbl))
  }

  features <- .task_constructors_resolve_features(data_tbl, features, c(time, event, id))
  .task_constructors_check_binary(data_tbl[[event]], event)

  time_vector <- data_tbl[[time]]
  event_vector <- data_tbl[[event]]
  if (!is.numeric(time_vector)) {
    rlang::abort("`time` column must be numeric.")
  }

  event_times <- time_vector[event_vector == 1]
  if (length(event_times) == 0) {
    rlang::abort("No observed events in data; unable to create survival task.")
  }
  max_event_time <- max(event_times, na.rm = TRUE)
  if (!is.finite(max_event_time)) {
    rlang::abort("Unable to determine maximum event time; check for missing or non-finite values.")
  }

  if (is.null(time_units)) {
    time_units <- NA_character_
  }
  time_range <- c(0, max_event_time)

  task <- list(
    data = data_tbl,
    id_col = id,
    time_col = time,
    event_col = event,
    cause_col = NULL,
    features = features,
    time_units = time_units,
    metadata = list(
      outcome = "survival",
      created = Sys.time()
    ),
    time_range = time_range
  )

  attr(task, "task_type") <- "survival"
  class(task) <- c("t2e_task_surv", "t2e_task")
  task
}

#' Create a competing risks analysis task
#'
#' @inheritParams ml4t2e_task_surv
#' @param status Name of the event indicator column (0 = censored).
#' @param cause Name of the event type column (factor/character/integer).
#'
#' @return An object of class `t2e_task_cr` (inherits from `t2e_task`).
#' @export
ml4t2e_task_cr <- function(data,
                           time,
                           status,
                           cause,
                           features = NULL,
                           id = NULL,
                           time_units = NULL) {
  .task_constructors_validate_arguments(time, status, id)
  if (!is.character(cause) || length(cause) != 1) {
    rlang::abort("`cause` must be a single column name.")
  }

  data_tbl <- dplyr::as_tibble(data)
  .task_constructors_check_columns(data_tbl, c(time, status, cause, id))

  if (is.null(id)) {
    id <- ".task_id"
    data_tbl[[id]] <- seq_len(nrow(data_tbl))
  }

  features <- .task_constructors_resolve_features(data_tbl, features, c(time, status, cause, id))
  .task_constructors_check_binary(data_tbl[[status]], status)

  time_vector <- data_tbl[[time]]
  if (!is.numeric(time_vector)) {
    rlang::abort("`time` column must be numeric.")
  }

  status_vector <- data_tbl[[status]]
  cause_vector <- data_tbl[[cause]]
  active_causes <- unique(cause_vector[status_vector == 1])
  if (length(active_causes) < 1) {
    rlang::abort("No observed events in data; unable to create competing risk task.")
  }

  cause_levels <- sort(unique(active_causes))
  numeric_levels <- suppressWarnings(as.numeric(cause_levels))
  if (!any(is.na(numeric_levels))) {
    cause_codes <- numeric_levels
  } else {
    cause_codes <- seq_along(cause_levels)
  }
  cause_map <- data.frame(
    cause = cause_levels,
    code = cause_codes,
    stringsAsFactors = FALSE
  )

  event_code_col <- ".event_code"
  event_codes <- rep(0, nrow(data_tbl))
  if (length(cause_levels) > 0) {
    match_idx <- match(as.character(cause_vector[status_vector == 1]), cause_map$cause)
    if (any(is.na(match_idx))) {
      rlang::abort("Unable to map all observed causes to numeric codes.")
    }
    event_codes[status_vector == 1] <- cause_map$code[match_idx]
  }
  data_tbl[[event_code_col]] <- event_codes

  max_event_time <- max(time_vector[status_vector == 1], na.rm = TRUE)
  if (!is.finite(max_event_time)) {
    rlang::abort("Unable to determine maximum event time; check for missing or non-finite values.")
  }

  if (is.null(time_units)) {
    time_units <- NA_character_
  }
  time_range <- c(0, max_event_time)

  task <- list(
    data = data_tbl,
    id_col = id,
    time_col = time,
    event_col = event_code_col,
    cause_col = cause,
    features = features,
    time_units = time_units,
    metadata = list(
      outcome = "competing_risk",
      created = Sys.time(),
      cause_map = cause_map,
      status_col = status
    ),
    time_range = time_range
  )

  attr(task, "task_type") <- "competing_risk"
  class(task) <- c("t2e_task_cr", "t2e_task")
  task
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------
.task_constructors_validate_arguments <- function(time, event, id) {
  if (!is.character(time) || length(time) != 1) {
    rlang::abort("`time` must be a single column name.")
  }
  if (!is.character(event) || length(event) != 1) {
    rlang::abort("Outcome columns must be single column names.")
  }
  if (!is.null(id) && (!is.character(id) || length(id) != 1)) {
    rlang::abort("`id` must be a single column name when provided.")
  }
}

.task_constructors_check_columns <- function(data, cols) {
  cols <- cols[!is.na(cols)]
  cols <- cols[cols != ""]
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    rlang::abort(paste0("Missing columns: ", paste(missing_cols, collapse = ", ")))
  }
}

.task_constructors_resolve_features <- function(data, features, reserved) {
  if (is.null(features)) {
    setdiff(colnames(data), reserved)
  } else {
    if (!is.character(features)) {
      rlang::abort("`features` must be a character vector of column names.")
    }
    missing_features <- setdiff(features, colnames(data))
    if (length(missing_features) > 0) {
      rlang::abort(paste0("`features` contains unknown columns: ",
                          paste(missing_features, collapse = ", ")))
    }
    features
  }
}

.task_constructors_check_binary <- function(x, name) {
  if (!is.numeric(x) && !is.integer(x) && !is.logical(x)) {
    rlang::abort(glue::glue("`{name}` must be numeric/integer/logical."))
  }
  unique_vals <- sort(unique(stats::na.omit(x)))
  if (!all(unique_vals %in% c(0, 1))) {
    rlang::abort(glue::glue("`{name}` must encode events as 0/1."))
  }
}
