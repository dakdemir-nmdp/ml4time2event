#' @title SHAP Variable Importance Plot
#'
#' @description Creates a variable importance plot showing which features
#' contribute most to predictions based on SHAP values.
#'
#' @param shap_result An object of class "ml4t2e_shap" from [ml4t2e_calculate_shap()].
#' @param max_features Maximum number of features to display. If NULL, shows all features.
#'   Features are ranked by mean absolute SHAP value.
#' @param plot_type Type of plot: "beeswarm" for scatter plot showing distribution,
#'   or "bar" for simple bar chart. Default is "beeswarm".
#'
#' @return A ggplot2 object showing variable importance.
#'
#' @details
#' The importance plot shows features ranked by their average impact on predictions.
#' In a beeswarm plot, each point represents one observation, colored by the
#' feature value (red = high, blue = low). The horizontal position shows the
#' SHAP value (contribution to prediction).
#'
#' In a bar plot, features are shown with their mean absolute SHAP value.
#'
#' @examples
#' \dontrun{
#' # Calculate SHAP values
#' shap_result <- ml4t2e_calculate_shap(pipeline, data, time_horizon = 365)
#'
#' # Create importance plot
#' ml4t2e_shap_importance(shap_result)
#'
#' # Show top 10 features only
#' ml4t2e_shap_importance(shap_result, max_features = 10)
#' }
#'
#' @export
ml4t2e_shap_importance <- function(shap_result,
                                    max_features = NULL,
                                    plot_type = c("beeswarm", "bar")) {

  # Validate inputs
  if (!inherits(shap_result, "ml4t2e_shap")) {
    stop("'shap_result' must be an ml4t2e_shap object from ml4t2e_calculate_shap().",
         call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2')",
         call. = FALSE)
  }

  plot_type <- match.arg(plot_type)

  # Extract SHAP values and feature values
  shap_vals <- shap_result$shap_values
  feature_vals <- shap_result$feature_values

  # Calculate mean absolute SHAP value for each feature
  importance <- colMeans(abs(shap_vals))
  importance_order <- order(importance, decreasing = TRUE)

  # Limit to top features if specified
  if (!is.null(max_features)) {
    if (!is.numeric(max_features) || max_features < 1) {
      stop("'max_features' must be a positive integer.", call. = FALSE)
    }
    max_features <- min(max_features, length(importance_order))
    importance_order <- importance_order[1:max_features]
  }

  selected_features <- names(importance)[importance_order]

  if (plot_type == "bar") {
    # Simple bar plot of mean absolute SHAP values
    plot_data <- data.frame(
      feature = factor(selected_features, levels = rev(selected_features)),
      importance = importance[importance_order]
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = importance, y = feature)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(
        title = "SHAP Variable Importance",
        subtitle = paste("Time horizon:", shap_result$time_horizon),
        x = "Mean |SHAP value|",
        y = "Feature"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 11)
      )

  } else {
    # Beeswarm plot showing distribution of SHAP values
    # Prepare long-format data
    plot_data <- data.frame()

    for (feat in selected_features) {
      feat_shap <- shap_vals[, feat]
      feat_val <- feature_vals[[feat]]

      # Normalize feature values for coloring (0-1 scale)
      if (is.numeric(feat_val)) {
        feat_val_norm <- (feat_val - min(feat_val, na.rm = TRUE)) /
          (max(feat_val, na.rm = TRUE) - min(feat_val, na.rm = TRUE) + 1e-10)
      } else {
        # For categorical, use factor levels
        feat_val_norm <- as.numeric(as.factor(feat_val)) /
          (length(unique(feat_val)) + 1e-10)
      }

      # Keep feature values as characters for categorical features to avoid factor
      # level clashes when combining across predictors.
      feat_val_display <- if (is.numeric(feat_val)) feat_val else as.character(feat_val)

      plot_data <- rbind(
        plot_data,
        data.frame(
          feature = feat,
          shap_value = feat_shap,
          feature_value = feat_val_display,
          feature_value_norm = feat_val_norm,
          stringsAsFactors = FALSE
        )
      )
    }

    # Order features by importance
    plot_data$feature <- factor(plot_data$feature,
                                 levels = rev(selected_features))

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = shap_value, y = feature, color = feature_value_norm)) +
      ggplot2::geom_jitter(alpha = 0.6, height = 0.2, size = 2) +
      ggplot2::scale_color_gradient(low = "blue", high = "red",
                                     name = "Feature\nvalue",
                                     breaks = c(0, 1),
                                     labels = c("Low", "High")) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = "SHAP Variable Importance",
        subtitle = paste("Time horizon:", shap_result$time_horizon),
        x = "SHAP value (impact on prediction)",
        y = "Feature"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 11),
        legend.position = "right"
      )
  }

  return(p)
}


#' @title SHAP Dependence Plot
#'
#' @description Creates a dependence plot showing how a specific feature
#' affects predictions across different feature values.
#'
#' @param shap_result An object of class "ml4t2e_shap" from [ml4t2e_calculate_shap()].
#' @param feature Name of the feature to plot. Must be one of the features
#'   in the SHAP result.
#' @param color_by Optional name of another feature to use for coloring points.
#'   If NULL (default), automatically selects the feature with highest interaction
#'   with the main feature.
#'
#' @return A ggplot2 object showing the dependence plot.
#'
#' @details
#' A dependence plot shows the relationship between a feature's value and its
#' impact on predictions (SHAP value). Points are colored by another feature
#' to reveal interactions.
#'
#' The x-axis shows the feature value, the y-axis shows the SHAP value
#' (contribution to prediction). The color reveals interactions with other features.
#'
#' @examples
#' \dontrun{
#' shap_result <- ml4t2e_calculate_shap(pipeline, data, time_horizon = 365)
#'
#' # Dependence plot for age
#' ml4t2e_shap_dependence(shap_result, feature = "age")
#'
#' # Specify interaction feature
#' ml4t2e_shap_dependence(shap_result, feature = "age", color_by = "sex")
#' }
#'
#' @export
ml4t2e_shap_dependence <- function(shap_result,
                                    feature,
                                    color_by = NULL) {

  # Validate inputs
  if (!inherits(shap_result, "ml4t2e_shap")) {
    stop("'shap_result' must be an ml4t2e_shap object from ml4t2e_calculate_shap().",
         call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.",
         call. = FALSE)
  }

  if (missing(feature) || !is.character(feature) || length(feature) != 1) {
    stop("'feature' must be a single character string specifying the feature name.",
         call. = FALSE)
  }

  if (!feature %in% shap_result$feature_names) {
    stop("Feature '", feature, "' not found in SHAP result. ",
         "Available features: ", paste(shap_result$feature_names, collapse = ", "),
         call. = FALSE)
  }

  # Extract SHAP values and feature values for the main feature
  shap_vals <- shap_result$shap_values[, feature]
  feat_vals <- shap_result$feature_values[[feature]]

  # Determine color_by feature if not specified
  if (is.null(color_by)) {
    # Find feature with highest interaction (correlation with SHAP values of main feature)
    other_features <- setdiff(shap_result$feature_names, feature)

    if (length(other_features) > 0) {
      interactions <- sapply(other_features, function(f) {
        feat_val <- shap_result$feature_values[[f]]
        if (is.numeric(feat_val)) {
          abs(cor(feat_val, shap_vals, use = "complete.obs"))
        } else {
          # For categorical, use variance in SHAP across levels
          shap_by_level <- tapply(shap_vals, feat_val, var, na.rm = TRUE)
          mean(shap_by_level, na.rm = TRUE)
        }
      })

      color_by <- names(which.max(interactions))
    } else {
      color_by <- NULL
    }
  } else {
    if (!color_by %in% shap_result$feature_names) {
      stop("color_by feature '", color_by, "' not found in SHAP result.",
           call. = FALSE)
    }
  }

  # Create plot data
  if (!is.null(color_by)) {
    color_vals <- shap_result$feature_values[[color_by]]

    # Normalize color values if numeric
    if (is.numeric(color_vals)) {
      color_vals_norm <- (color_vals - min(color_vals, na.rm = TRUE)) /
        (max(color_vals, na.rm = TRUE) - min(color_vals, na.rm = TRUE) + 1e-10)
    } else {
      color_vals_norm <- as.numeric(as.factor(color_vals))
    }

    plot_data <- data.frame(
      feature_value = feat_vals,
      shap_value = shap_vals,
      color_value = color_vals,
      color_value_norm = color_vals_norm
    )

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = feature_value, y = shap_value, color = color_value_norm)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::scale_color_gradient(low = "blue", high = "red",
                                     name = paste0(color_by, "\nvalue")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

  } else {
    plot_data <- data.frame(
      feature_value = feat_vals,
      shap_value = shap_vals
    )

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = feature_value, y = shap_value)) +
      ggplot2::geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  }

  # Add smoothing line if feature is numeric
  if (is.numeric(feat_vals)) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, color = "black",
                                   linetype = "solid", linewidth = 0.8)
  }

  p <- p +
    ggplot2::labs(
      title = paste("SHAP Dependence Plot:", feature),
      subtitle = paste("Time horizon:", shap_result$time_horizon,
                       if (!is.null(color_by)) paste("| Colored by:", color_by) else ""),
      x = feature,
      y = "SHAP value (impact on prediction)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11)
    )

  return(p)
}


#' @title SHAP Waterfall Plot for Individual Prediction
#'
#' @description Creates a waterfall plot explaining an individual prediction
#' by showing how each feature contributes.
#'
#' @param shap_result An object of class "ml4t2e_shap" from [ml4t2e_calculate_shap()].
#' @param obs_id Integer index of the observation to explain (row number in the
#'   original data).
#' @param max_features Maximum number of features to display. If NULL, shows all.
#'   Features are ranked by absolute SHAP value.
#'
#' @return A ggplot2 object showing the waterfall plot.
#'
#' @details
#' A waterfall plot shows how the prediction for a single observation builds up
#' from the baseline (average prediction) by adding the contribution of each feature.
#'
#' The plot shows:
#' \itemize{
#'   \item Baseline: The average prediction across all observations
#'   \item Feature contributions: How each feature pushes the prediction higher or lower
#'   \item Final prediction: The actual prediction for this observation
#' }
#'
#' Formula: Prediction = Baseline + sum(SHAP values)
#'
#' @examples
#' \dontrun{
#' shap_result <- ml4t2e_calculate_shap(pipeline, data, time_horizon = 365)
#'
#' # Explain first observation
#' ml4t2e_shap_waterfall(shap_result, obs_id = 1)
#'
#' # Show top 10 features only
#' ml4t2e_shap_waterfall(shap_result, obs_id = 5, max_features = 10)
#' }
#'
#' @export
ml4t2e_shap_waterfall <- function(shap_result,
                                   obs_id,
                                   max_features = NULL) {

  # Validate inputs
  if (!inherits(shap_result, "ml4t2e_shap")) {
    stop("'shap_result' must be an ml4t2e_shap object from ml4t2e_calculate_shap().",
         call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.",
         call. = FALSE)
  }

  if (missing(obs_id) || !is.numeric(obs_id) || length(obs_id) != 1) {
    stop("'obs_id' must be a single integer.", call. = FALSE)
  }

  if (obs_id < 1 || obs_id > nrow(shap_result$shap_values)) {
    stop("'obs_id' out of range. Must be between 1 and ",
         nrow(shap_result$shap_values), ".", call. = FALSE)
  }

  if (obs_id != round(obs_id) || obs_id <= 0) {
    stop("'obs_id' must be a positive integer.", call. = FALSE)
  }

  # Extract SHAP values for this observation
  obs_shap <- shap_result$shap_values[obs_id, ]

  # Rank features by absolute SHAP value
  feature_order <- order(abs(obs_shap), decreasing = TRUE)

  # Limit to top features if specified
  if (!is.null(max_features)) {
    if (!is.numeric(max_features) || max_features < 1) {
      stop("'max_features' must be a positive integer.", call. = FALSE)
    }
    max_features <- min(max_features, length(feature_order))
    feature_order <- feature_order[1:max_features]

    # Sum remaining features as "Other"
    other_shap <- sum(obs_shap[-feature_order])
  } else {
    other_shap <- 0
  }

  # Get selected features and their SHAP values
  selected_features <- names(obs_shap)[feature_order]
  selected_shap <- obs_shap[feature_order]

  # Build waterfall data
  baseline <- shap_result$baseline
  final_pred <- shap_result$predictions[obs_id]

  # Create components: baseline, features, other (if any), final
  components <- c("Baseline", selected_features)
  values <- c(baseline, selected_shap)

  if (!is.null(max_features) && abs(other_shap) > 1e-10) {
    components <- c(components, "Other features")
    values <- c(values, other_shap)
  }

  components <- c(components, "Prediction")
  values <- c(values, final_pred)

  # Calculate cumulative values for waterfall
  n <- length(components)
  start_vals <- numeric(n)
  end_vals <- numeric(n)

  start_vals[1] <- 0
  end_vals[1] <- baseline

  for (i in 2:(n - 1)) {
    start_vals[i] <- end_vals[i - 1]
    end_vals[i] <- start_vals[i] + values[i]
  }

  # Final prediction
  start_vals[n] <- 0
  end_vals[n] <- final_pred

  # Create plot data
  plot_data <- data.frame(
    component = factor(components, levels = components),
    start = start_vals,
    end = end_vals,
    value = values,
    contribution = c(NA, values[2:(n - 1)], NA),
    is_positive = c(NA, values[2:(n - 1)] > 0, NA),
    type = c("baseline", rep("feature", n - 2), "prediction")
  )

  # Create color mapping
  colors <- ifelse(plot_data$type == "baseline", "gray",
                   ifelse(plot_data$type == "prediction", "purple",
                          ifelse(plot_data$is_positive, "red", "blue")))

  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_rect(ggplot2::aes(xmin = as.numeric(component) - 0.4,
                                    xmax = as.numeric(component) + 0.4,
                                    ymin = start,
                                    ymax = end,
                                    fill = component)) +
    ggplot2::scale_fill_manual(values = stats::setNames(colors, components),
                               guide = "none") +
    ggplot2::geom_segment(
      data = plot_data[2:n, ],
      ggplot2::aes(x = as.numeric(component) - 0.4,
                   xend = as.numeric(component) - 0.6,
                   y = start,
                   yend = start),
      linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = as.numeric(component),
                   y = (start + end) / 2,
                   label = sprintf("%.3f", value)),
      size = 3
    ) +
    ggplot2::labs(
      title = paste("SHAP Waterfall Plot - Observation", obs_id),
      subtitle = paste("Time horizon:", shap_result$time_horizon,
                       "| Baseline:", sprintf("%.3f", baseline),
                       "| Prediction:", sprintf("%.3f", final_pred)),
      x = "Component",
      y = "Expected Time Lost"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = ggplot2::element_blank()
    )

  return(p)
}
