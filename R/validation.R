#' Input Validation Utilities
#'
#' Comprehensive validation functions with informative error messages.
#' Uses checkmate for type validation when available, falls back to base R.
#'
#' @keywords internal
#' @noRd

#' Validate task object
#'
#' Ensures task has required structure and fields.
#'
#' @param task Task object to validate
#' @param required_fields Character vector of required field names
#' @param context Character; where this is being called (for error messages)
#' @return Invisible TRUE or aborts with informative error
#' @keywords internal
.validate_task <- function(task, required_fields = c("data", "time_col", "event_col"), context = "function") {
    if (is.null(task)) {
        rlang::abort(
            glue::glue("Task is NULL in {context}."),
            class = "invalid_task_error"
        )
    }

    if (!inherits(task, "t2e_task")) {
        rlang::abort(
            glue::glue(
                "Invalid task object in {context}.\n",
                "Expected: object of class 't2e_task'\n",
                "Actual: {paste(class(task), collapse=', ')}\n",
                "Hint: Create task with ml4t2e_task_surv() or ml4t2e_task_cr()"
            ),
            class = "invalid_task_error"
        )
    }

    missing_fields <- setdiff(required_fields, names(task))
    if (length(missing_fields) > 0) {
        rlang::abort(
            glue::glue(
                "Task is missing required fields: {paste(missing_fields, collapse=', ')}\n",
                "Context: {context}\n",
                "Available fields: {paste(names(task), collapse=', ')}"
            ),
            class = "invalid_task_error"
        )
    }

    # Validate data field specifically
    if (!is.null(task$data) && !is.data.frame(task$data)) {
        rlang::abort(
            glue::glue(
                "Task$data must be a data.frame, got {class(task$data)[1]}"
            ),
            class = "invalid_task_error"
        )
    }

    invisible(TRUE)
}

#' Validate model names against registry
#'
#' @param models Character vector of model names
#' @param outcome_type Character; "survival" or "competing_risk"
#' @param context Character; calling context
#' @return Invisible TRUE or aborts
#' @keywords internal
.validate_models <- function(models, outcome_type, context = "ml4t2e_fit") {
    if (length(models) == 0) {
        rlang::abort(
            glue::glue(
                "No models specified in {context}.\n",
                "At least one model name is required.\n",
                "See ml4t2e_list_models('{outcome_type}') for available options."
            ),
            class = "invalid_argument_error"
        )
    }

    if (!is.character(models)) {
        rlang::abort(
            glue::glue(
                "Models must be a character vector, got {class(models)[1]}\n",
                "Example: models = c('cox', 'random_forest')"
            ),
            class = "invalid_argument_error"
        )
    }

    # Check against registry
    available <- ml4t2e_list_models(outcome_type)$engine
    invalid <- setdiff(models, available)

    if (length(invalid) > 0) {
        rlang::abort(
            glue::glue(
                "Unknown model(s) for {outcome_type}: {paste(invalid, collapse=', ')}\n",
                "Available models: {paste(head(available, 10), collapse=', ')}",
                if (length(available) > 10) paste0(" (and ", length(available) - 10, " more)") else "",
                "\nSee: ml4t2e_list_models('{outcome_type}')"
            ),
            class = "invalid_model_error"
        )
    }

    invisible(TRUE)
}

#' Validate time grid
#'
#' @param times Numeric vector of times or NULL
#' @param allow_null Logical; whether NULL is acceptable
#' @param min_val Numeric; minimum allowed value (default 0)
#' @param max_val Numeric; maximum allowed value
#' @param context Character
#' @return Invisible TRUE or aborts
#' @keywords internal
.validate_times <- function(times, allow_null = TRUE, min_val = 0, max_val = Inf, context = "function") {
    if (is.null(times)) {
        if (!allow_null) {
            rlang::abort(
                glue::glue("Times cannot be NULL in {context}."),
                class = "invalid_argument_error"
            )
        }
        return(invisible(TRUE))
    }

    if (!is.numeric(times)) {
        rlang::abort(
            glue::glue(
                "Times must be numeric, got {class(times)[1]} in {context}."
            ),
            class = "invalid_argument_error"
        )
    }

    if (any(is.na(times))) {
        rlang::abort(
            glue::glue(
                "Times contains NA values in {context}.\n",
                "Position(s): {paste(which(is.na(times)), collapse=', ')}"
            ),
            class = "invalid_argument_error"
        )
    }

    if (any(times < min_val)) {
        bad_idx <- which(times < min_val)
        rlang::abort(
            glue::glue(
                "Times contains values below minimum ({min_val}) in {context}.\n",
                "Violating values at positions: {paste(head(bad_idx, 5), collapse=', ')}\n",
                "Values: {paste(head(times[bad_idx], 5), collapse=', ')}"
            ),
            class = "invalid_argument_error"
        )
    }

    if (any(times > max_val)) {
        bad_idx <- which(times > max_val)
        rlang::abort(
            glue::glue(
                "Times contains values above maximum ({max_val}) in {context}.\n",
                "Violating values at positions: {paste(head(bad_idx, 5), collapse=', ')}"
            ),
            class = "invalid_argument_error"
        )
    }

    invisible(TRUE)
}

#' Validate data quality
#'
#' Checks for common data issues in time-to-event analysis.
#'
#' @param data Data frame
#' @param required_cols Character vector of required column names
#' @param id_col Character; name of ID column (optional)
#' @param time_col Character; name of time column (optional, for range checks)
#' @param event_col Character; name of event column (optional, for checks)
#' @param context Character
#' @return Invisible TRUE, aborts on critical errors, warns on minor issues
#' @keywords internal
.validate_data_quality <- function(data, required_cols, id_col = NULL, time_col = NULL,
                                   event_col = NULL, context = "function") {
    if (!is.data.frame(data)) {
        rlang::abort(
            glue::glue("Data must be a data.frame, got {class(data)[1]} in {context}."),
            class = "invalid_argument_error"
        )
    }

    if (nrow(data) == 0) {
        rlang::abort(
            glue::glue("Data is empty (0 rows) in {context}."),
            class = "invalid_data_error"
        )
    }

    # Check required columns
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        rlang::abort(
            glue::glue(
                "Required column(s) missing in {context}: {paste(missing_cols, collapse=', ')}\n",
                "Available columns: {paste(names(data), collapse=', ')}"
            ),
            class = "missing_column_error"
        )
    }

    # Check ID uniqueness if specified
    if (!is.null(id_col) && id_col %in% names(data)) {
        ids <- data[[id_col]]
        if (any(duplicated(ids))) {
            dup_ids <- ids[duplicated(ids)]
            rlang::abort(
                glue::glue(
                    "Duplicate IDs found in {context}.\n",
                    "ID column: '{id_col}'\n",
                    "First few duplicates: {paste(head(unique(dup_ids), 5), collapse=', ')}\n",
                    "Count: {length(unique(dup_ids))} unique IDs are duplicated"
                ),
                class = "duplicate_id_error"
            )
        }
    }

    # Check time column
    if (!is.null(time_col) && time_col %in% names(data)) {
        times <- data[[time_col]]

        if (any(is.na(times))) {
            n_na <- sum(is.na(times))
            rlang::abort(
                glue::glue(
                    "Time column contains {n_na} NA value(s) in {context}.\n",
                    "Time column: '{time_col}'\n",
                    "Row positions: {paste(head(which(is.na(times)), 10), collapse=', ')}"
                ),
                class = "missing_data_error"
            )
        }

        if (any(times < 0, na.rm = TRUE)) {
            n_neg <- sum(times < 0, na.rm = TRUE)
            rlang::abort(
                glue::glue(
                    "Time column contains {n_neg} negative value(s) in {context}.\n",
                    "Time column: '{time_col}'\n",
                    "First few negative times: {paste(head(times[times < 0], 5), collapse=', ')}"
                ),
                class = "invalid_data_error"
            )
        }

        # Warn if all times are the same
        if (length(unique(times)) == 1) {
            rlang::warn(
                glue::glue(
                    "All times are identical ({unique(times)[1]}) in {context}.\n",
                    "This may cause model fitting issues."
                )
            )
        }
    }

    # Check event column
    if (!is.null(event_col) && event_col %in% names(data)) {
        events <- data[[event_col]]

        if (any(is.na(events))) {
            n_na <- sum(is.na(events))
            rlang::warn(
                glue::glue(
                    "Event column contains {n_na} NA value(s) in {context}.\n",
                    "These will be treated as censored observations.\n",
                    "Event column: '{event_col}'"
                )
            )
        }

        # Check if all censored
        if (all(events == 0 | is.na(events))) {
            rlang::abort(
                glue::glue(
                    "All observations are censored (event=0 or NA) in {context}.\n",
                    "Cannot fit time-to-event model without any events.\n",
                    "Event column: '{event_col}'"
                ),
                class = "no_events_error"
            )
        }

        # Warn if very few events
        n_events <- sum(events > 0, na.rm = TRUE)
        if (n_events < 5) {
            rlang::warn(
                glue::glue(
                    "Very few events ({n_events}) in {context}.\n",
                    "Model may be unstable. Consider collecting more data."
                )
            )
        }
    }

    invisible(TRUE)
}

#' Validate conformal alpha parameter
#'
#' @param alpha Numeric; coverage level (e.g., 0.1 for 90% confidence)
#' @param allow_null Logical
#' @param context Character
#' @return Invisible TRUE or aborts
#' @keywords internal
.validate_conformal_alpha <- function(alpha, allow_null = TRUE, context = "predict") {
    if (is.null(alpha)) {
        if (!allow_null) {
            rlang::abort(
                glue::glue("Conformal alpha cannot be NULL in {context}."),
                class = "invalid_argument_error"
            )
        }
        return(invisible(TRUE))
    }

    if (!is.numeric(alpha) || length(alpha) != 1) {
        rlang::abort(
            glue::glue(
                "Conformal alpha must be a single numeric value in {context}.\n",
                "Got: {class(alpha)[1]} of length {length(alpha)}\n",
                "Example: conformal_alpha = 0.1  # for 90% confidence bands"
            ),
            class = "invalid_argument_error"
        )
    }

    if (alpha <= 0 || alpha >= 1) {
        rlang::abort(
            glue::glue(
                "Conformal alpha must be between 0 and 1 (exclusive) in {context}.\n",
                "Got: {alpha}\n",
                "Common values: 0.05 (95% confidence), 0.1 (90% confidence), 0.2 (80% confidence)"
            ),
            class = "invalid_argument_error"
        )
    }

    invisible(TRUE)
}

#' Validate ensemble strategy
#'
#' @param strategy Character; ensemble strategy name
#' @param available_strategies Character vector of valid strategies
#' @param context Character
#' @return Invisible TRUE or aborts
#' @keywords internal
.validate_ensemble_strategy <- function(strategy,
                                        available_strategies = c("simple", "stack", "none", "auto"),
                                        context = "ml4t2e_fit") {
    if (!is.character(strategy) || length(strategy) != 1) {
        rlang::abort(
            glue::glue(
                "Ensemble strategy must be a single character value in {context}.\n",
                "Got: {class(strategy)[1]} of length {length(strategy)}"
            ),
            class = "invalid_argument_error"
        )
    }

    if (!tolower(strategy) %in% tolower(available_strategies)) {
        rlang::abort(
            glue::glue(
                "Unknown ensemble strategy '{strategy}' in {context}.\n",
                "Available: {paste(available_strategies, collapse=', ')}"
            ),
            class = "invalid_argument_error"
        )
    }

    invisible(TRUE)
}

#' Validate prediction output structure
#'
#' Ensures prediction objects have required fields and valid data.
#'
#' @param preds Prediction object to validate
#' @param expected_type Character; "survival" or "cif"
#' @param context Character
#' @return Invisible TRUE or aborts
#' @keywords internal
.validate_predictions <- function(preds, expected_type = c("survival", "cif"), context = "function") {
    expected_type <- match.arg(expected_type)

    if (!is.data.frame(preds)) {
        rlang::abort(
            glue::glue(
                "Predictions must be a data.frame, got {class(preds)[1]} in {context}."
            ),
            class = "invalid_prediction_error"
        )
    }

    required_cols <- switch(expected_type,
        survival = c("id", "time", "surv", "model"),
        cif = c("id", "time", "cause", "cif", "model")
    )

    missing <- setdiff(required_cols, names(preds))
    if (length(missing) > 0) {
        rlang::abort(
            glue::glue(
                "Prediction object missing required columns: {paste(missing, collapse=', ')}\n",
                "Expected type: {expected_type}\n",
                "Context: {context}"
            ),
            class = "invalid_prediction_error"
        )
    }

    # Check for invalid probabilities
    prob_col <- if (expected_type == "survival") "surv" else "cif"
    probs <- preds[[prob_col]]

    if (any(probs < 0 | probs > 1, na.rm = TRUE)) {
        bad_idx <- which(probs < 0 | probs > 1)
        rlang::warn(
            glue::glue(
                "Predictions contain invalid probabilities in {context}.\n",
                "Column: '{prob_col}'\n",
                "Invalid rows: {paste(head(bad_idx, 5), collapse=', ')}\n",
                "Values: {paste(head(probs[bad_idx], 5), collapse=', ')}"
            )
        )
    }

    invisible(TRUE)
}
