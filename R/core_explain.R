#' Explain model predictions via permutation importance
#'
#' @name ml4t2e_explain
#' @param fit A `t2e_fit` object.
#' @param newdata Optional data frame. If `NULL`, training data stored in the
#'   task is used.
#' @param method Explanation method. Only `"permutation"` is currently
#'   implemented; `"auto"` maps to `"permutation"`.
#' @param include Models to include (same semantics as `predict()`).
#' @param times Optional time grid used for predictions (defaults to fit's grid).
#' @param top_n Optionally limit the number of features returned (by absolute
#'   importance).
#'
#' @return A tibble with columns `model`, `feature`, `importance`, and
#'   `metric`.
#' @export

utils::globalVariables(c("feature", "importance", "model"))

ml4t2e_explain <- function(fit,
                           newdata = NULL,
                           method = c("auto", "permutation"),
                           include = "all",
                           times = NULL,
                           top_n = NULL) {
  if (!inherits(fit, "t2e_fit")) {
    rlang::abort("`fit` must be a `t2e_fit` object.")
  }
  outcome_type <- fit[["outcome_type"]]
  if (!identical(outcome_type, "survival")) {
    rlang::abort("Explanation is currently implemented for survival outcomes only.")
  }

  method <- match.arg(method)
  if (identical(method, "auto")) {
    method <- "permutation"
  }
  if (!identical(method, "permutation")) {
    rlang::abort("Unsupported explanation method.")
  }

  task <- fit[["task"]]
  if (is.null(newdata)) {
    data <- task$data
    eval_task <- task
  } else {
    data <- dplyr::as_tibble(newdata)
    .verify_new_data_columns(data, task)
    eval_task <- ml4t2e_task_surv(
      data = data,
      time = task$time_col,
      event = task$event_col,
      features = task$features,
      time_units = task$time_units
    )
    data <- eval_task$data
  }

  include <- .normalize_include(include, fit[["model_names"]])
  times <- .resolve_prediction_times(times, fit[["time_grid"]])

  baseline_preds <- predict(fit, newdata = data, times = times, include = include, type = "survival")
  baseline_metrics <- ml4t2e_evaluate(baseline_preds, task = eval_task, metrics = "c_index")

  features <- task$features
  results <- list()

  for (feature in features) {
    permuted <- data
    permuted[[feature]] <- sample(permuted[[feature]])
    perm_preds <- predict(fit, newdata = permuted, times = times, include = include, type = "survival")
    perm_metrics <- ml4t2e_evaluate(perm_preds, task = eval_task, metrics = "c_index")

    merged <- dplyr::left_join(
      baseline_metrics,
      perm_metrics,
      by = c("model", "metric"),
      suffix = c("_baseline", "_permuted")
    )
    merged$feature <- feature
    merged$importance <- merged$value_baseline - merged$value_permuted
    merged$method <- "permutation"
    merged$value_baseline <- NULL
    merged$value_permuted <- NULL
    results[[length(results) + 1]] <- merged
  }

  output <- dplyr::bind_rows(results)
  if (!is.null(top_n)) {
    pieces <- split(output, output$model)
    trimmed <- lapply(pieces, function(df) {
      dplyr::slice_max(df, order_by = abs(df$importance), n = top_n, with_ties = FALSE)
    })
    output <- dplyr::bind_rows(trimmed)
  }
  class(output) <- unique(c("t2e_explain", class(output)))
  output
}

#' Autoplot method for t2e_explain objects
#'
#' @param object A `t2e_explain` object.
#' @param what What to plot (currently only "importance" is supported).
#' @param ... Additional arguments passed to ggplot2.
#'
#' @return A ggplot2 object.
#' @export
autoplot.t2e_explain <- function(object,
                                 what = c("importance"),
                                 ...) {
  what <- match.arg(what)
  ggplot2::ggplot(object, ggplot2::aes(
    x = stats::reorder(feature, importance),
    y = importance,
    fill = model
  )) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Feature", y = "Importance (\u0394 c-index)", fill = "Model") +
    ggplot2::theme_minimal()
}
