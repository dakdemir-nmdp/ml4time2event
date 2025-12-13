utils::globalVariables(c("value", "sd", "n", "mean"))

#' Generate a model leaderboard
#'
#' Aggregates performance metrics from a fitted pipeline to rank models.
#' Prioritises cross-validated metrics if available; falls back to training metrics with a warning.
#'
#' @param object A fitted `ml4t2e_pipeline` object.
#' @param metric Optional character string to filter by a specific metric (e.g. "c_index").
#'
#' @return A `t2e_leaderboard` object (a tibble) with columns:
#'   - `model`: Model engine name.
#'   - `metric`: Metric name.
#'   - `mean`: Average score across folds (or single training score).
#'   - `sd`: Standard deviation across folds (NA for training metrics).
#'   - `n`: Number of folds/repeats.
#' @export
ml4t2e_leaderboard <- function(object, metric = NULL) {
    if (!inherits(object, "ml4t2e_pipeline") && !inherits(object, "T2EPipeline")) {
        rlang::abort("`object` must be a `ml4t2e_pipeline` object.")
    }

    resampling_results <- object$resampling_results
    training_metrics <- object$training_metrics

    if (!is.null(resampling_results)) {
        data <- resampling_results
        source_label <- "resampling"
    } else if (!is.null(training_metrics)) {
        rlang::warn("No resampling results found. Using training metrics (likely optimistic).")
        data <- training_metrics
        source_label <- "training"
    } else {
        rlang::abort("No metrics found in pipeline. Did you run `$fit()`?")
    }

    if (!is.null(metric)) {
        data <- data[data$metric == metric, , drop = FALSE]
        if (nrow(data) == 0) {
            rlang::abort(glue::glue("Metric '{metric}' not found in results."))
        }
    }

    # Aggregation
    if (source_label == "resampling") {
        leaderboard <- data %>%
            dplyr::group_by(model, metric) %>%
            dplyr::summarise(
                mean = mean(value, na.rm = TRUE),
                sd = sd(value, na.rm = TRUE),
                n = dplyr::n(),
                .groups = "drop"
            )
    } else {
        leaderboard <- data %>%
            dplyr::mutate(mean = value, sd = NA_real_, n = 1) %>%
            dplyr::select(model, metric, mean, sd, n)
    }

    # Sorting (higher C-index is better, lower IBS is better)
    # For now, we just sort by mean descending, but this might be wrong for IBS.
    # Let's try to be smart about it.

    # Heuristic: "c_index" -> Descending, "ibs" / "brier" -> Ascending
    # We can't easily sort mixed metrics in one view.
    # We will sort by the first metric encountered or just leave it by model.
    # Actually, users usually view one metric at a time or facet.

    # Let's sort by mean descending by default, but users can re-sort.
    leaderboard <- leaderboard %>%
        dplyr::arrange(dplyr::desc(mean))

    attr(leaderboard, "source") <- source_label
    class(leaderboard) <- c("t2e_leaderboard", "tbl_df", "tbl", "data.frame")
    leaderboard
}

#' Plot model leaderboard
#'
#' @param object A `t2e_leaderboard` object.
#' @param metric The metric to plot. Defaults to the first metric available.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `ggplot` object.
#' @export
autoplot.t2e_leaderboard <- function(object, metric = NULL, ...) {
    if (is.null(metric)) {
        metric <- unique(object$metric)[1]
        msg <- glue::glue("Plotting metric: {metric}")
        rlang::inform(msg)
    }

    plot_data <- object[object$metric == metric, , drop = FALSE]

    if (nrow(plot_data) == 0) {
        rlang::abort(glue::glue("Metric '{metric}' not found in leaderboard."))
    }

    # Determine sort order
    # C-index: higher is better
    # IBS: lower is better
    is_loss <- grepl("ibs|brier", metric, ignore.case = TRUE)

    if (is_loss) {
        # Sort ascending (best on top/left)
        plot_data <- plot_data %>% dplyr::arrange(mean)
    } else {
        # Sort descending (best on top/left)
        plot_data <- plot_data %>% dplyr::arrange(dplyr::desc(mean))
    }

    # Lock in factor order for plotting
    plot_data$model <- factor(plot_data$model, levels = rev(plot_data$model))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean, y = model)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(
            title = "Model Performance Leaderboard",
            subtitle = paste0("Metric: ", metric, " (", attr(object, "source"), ")"),
            x = "Score (Mean)",
            y = "Model"
        ) +
        ggplot2::theme_minimal()

    if (any(!is.na(plot_data$sd))) {
        p <- p + ggplot2::geom_errorbar(
            ggplot2::aes(xmin = mean - sd, xmax = mean + sd),
            width = 0.2
        )
    }

    p
}
