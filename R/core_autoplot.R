#' @importFrom utils globalVariables
utils::globalVariables(c("time", "surv", "model", "cause", "cif"))

#' @importFrom ggplot2 autoplot
#' @export
autoplot.t2e_pred <- function(object,
                              what = c("curves"),
                              include = "all",
                              alpha = 0.4,
                              ...) {
  what <- match.arg(what)
  pred_type <- attr(object, "pred_type", exact = TRUE)
  available_models <- unique(object[["model"]])
  include <- .normalize_include(include, available_models)
  data <- object[object[["model"]] %in% include, , drop = FALSE]

  if (nrow(data) == 0) {
    rlang::abort("No predictions available to plot for the requested models.")
  }

  if (identical(pred_type, "survival") || is.null(pred_type)) {
    plot_data <- data %>%
      dplyr::group_by(model, time) %>%
      dplyr::summarise(surv = mean(surv, na.rm = TRUE), .groups = "drop")

    ggplot2::ggplot(plot_data, ggplot2::aes(
      x = time,
      y = surv,
      colour = model
    )) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Time", y = "Survival probability", colour = "Model") +
      ggplot2::theme_minimal()
  } else if (identical(pred_type, "cif")) {
    plot_data <- data %>%
      dplyr::group_by(model, cause, time) %>%
      dplyr::summarise(cif = mean(cif, na.rm = TRUE), .groups = "drop")

    ggplot2::ggplot(plot_data, ggplot2::aes(
      x = time,
      y = cif,
      colour = model
    )) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~cause, scales = "free_y") +
      ggplot2::labs(x = "Time", y = "Cumulative incidence", colour = "Model", title = "CIF predictions") +
      ggplot2::theme_minimal()
  } else {
    rlang::abort(glue::glue("Unsupported prediction type '{pred_type}' for autoplot."))
  }
}

#' @export
autoplot.t2e_fit <- function(object,
                             what = c("curves"),
                             include = "all",
                             ...) {
  outcome <- object$outcome_type
  pred_type <- if (identical(outcome, "competing_risk")) "cif" else "survival"
  preds <- predict(object, include = include, type = pred_type, ...)
  autoplot(preds, what = what, include = include, ...)
}
