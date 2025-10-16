#' @title plot_survival_curves
#'
#' @description Plot survival curves from one or more survival model predictions
#'
#' @param predictions A list of prediction objects from survival models, or a single prediction object.
#'   Each prediction object should have components 'Probs' (survival probabilities matrix) 
#'   and 'Times' (time points vector).
#' @param model_names Character vector of model names. If NULL, uses names from predictions list.
#' @param patients_to_plot Integer vector specifying which patients (row indices) to plot.
#'   If NULL, plots first 3 patients.
#' @param colors Named character vector of colors for each model. If NULL, uses default colors.
#' @param highlight_ensemble Logical, whether to highlight ensemble predictions in black (default TRUE).
#' @param title Character string for plot title. If NULL, uses default title.
#' @param subtitle Character string for plot subtitle. If NULL, uses default subtitle.
#' @param ncol_facets Integer, number of columns for facetting by patient (default 3).
#' @param add_median_line Logical, whether to add horizontal line at 0.5 survival (default FALSE).
#' @param legend_position Character, position of legend ("bottom", "right", "top", "left", "none").
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Single model
#' cox_pred <- Predict_SurvModel_Cox(cox_model, test_data)
#' plot_survival_curves(cox_pred, model_names = "Cox")
#' 
#' # Multiple models
#' predictions <- list(
#'   Cox = Predict_SurvModel_Cox(cox_model, test_data),
#'   RF = Predict_SurvModel_RF(rf_model, test_data)
#' )
#' plot_survival_curves(predictions, patients_to_plot = 1:5)
#' }
#'
#' @importFrom rlang .data
#' @importFrom grDevices rainbow
#' @importFrom stats median
#' @export
plot_survival_curves <- function(predictions,
                                 model_names = NULL,
                                 patients_to_plot = NULL,
                                 colors = NULL,
                                 highlight_ensemble = TRUE,
                                 title = NULL,
                                 subtitle = NULL,
                                 ncol_facets = 3,
                                 add_median_line = FALSE,
                                 legend_position = "bottom") {
  
  # Input validation
  if (is.null(predictions)) {
    stop("'predictions' cannot be NULL")
  }
  
  # Handle single prediction object
  if (!is.list(predictions) || (!is.null(predictions$Probs) && !is.null(predictions$Times))) {
    predictions <- list(Model1 = predictions)
  }
  
  # Validate prediction objects
  for (i in seq_along(predictions)) {
    pred <- predictions[[i]]
    if (is.null(pred$Probs) || is.null(pred$Times)) {
      stop("Each prediction object must have 'Probs' and 'Times' components")
    }
    if (!is.matrix(pred$Probs)) {
      stop("'Probs' component must be a matrix")
    }
    if (!is.numeric(pred$Times)) {
      stop("'Times' component must be numeric")
    }
  }
  
  # Set model names
  if (is.null(model_names)) {
    if (!is.null(names(predictions))) {
      model_names <- names(predictions)
    } else {
      model_names <- paste0("Model", seq_along(predictions))
      names(predictions) <- model_names
    }
  } else {
    if (length(model_names) != length(predictions)) {
      stop("Length of 'model_names' must match length of 'predictions'")
    }
    names(predictions) <- model_names
  }
  
  # Determine patients to plot
  n_patients <- ncol(predictions[[1]]$Probs)
  if (is.null(patients_to_plot)) {
    patients_to_plot <- seq_len(min(3, n_patients))
  } else {
    if (max(patients_to_plot) > n_patients) {
      stop("'patients_to_plot' contains indices larger than number of patients")
    }
    if (min(patients_to_plot) < 1) {
      stop("'patients_to_plot' must contain positive integers")
    }
  }
  
  # Create plot data
  plot_data <- data.frame()
  
  for (patient in patients_to_plot) {
    for (model_name in model_names) {
      pred <- predictions[[model_name]]
      model_data <- data.frame(
        Time = pred$Times,
        Survival = pred$Probs[, patient],
        Model = model_name,
        Patient = paste("Patient", patient),
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, model_data)
    }
  }
  
  # Set colors
  if (is.null(colors)) {
    unique_models <- unique(plot_data$Model)
    colors <- rainbow(length(unique_models))
    names(colors) <- unique_models
    
    # Highlight ensemble if requested
    if (highlight_ensemble && "Ensemble" %in% unique_models) {
      colors["Ensemble"] <- "black"
    }
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Survival, color = .data$Model)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors)
  
  # Add facetting if multiple patients
  if (length(patients_to_plot) > 1) {
    p <- p + ggplot2::facet_wrap(~Patient, ncol = ncol_facets)
  }
  
  # Add median line if requested
  if (add_median_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.7)
  }
  
  # Set title and subtitle
  if (is.null(title)) {
    title <- "Survival Curves"
    if (length(model_names) > 1) {
      title <- paste(title, "- Model Comparison")
    }
  }
  
  if (is.null(subtitle)) {
    subtitle <- paste("Comparing", length(model_names), "model(s)")
    if (length(patients_to_plot) > 1) {
      subtitle <- paste(subtitle, "across", length(patients_to_plot), "patients")
    }
  }
  
  p <- p + ggplot2::labs(title = title, subtitle = subtitle, 
                          x = "Time", y = "Survival Probability") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = legend_position,
                   plot.title = ggplot2::element_text(face = "bold"))
  
  # Adjust legend for multiple models
  if (length(model_names) > 4 && legend_position == "bottom") {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(ncol = min(5, length(model_names))))
  }
  
  return(p)
}


#' @title plot_cif_curves
#'
#' @description Plot cumulative incidence function (CIF) curves from one or more competing risks model predictions
#'
#' @param predictions A list of prediction objects from competing risks models, or a single prediction object.
#'   Each prediction object should have components 'CIFs' (cumulative incidence matrix) 
#'   and 'Times' (time points vector).
#' @param model_names Character vector of model names. If NULL, uses names from predictions list.
#' @param patients_to_plot Integer vector specifying which patients (row indices) to plot.
#'   If NULL, plots first 3 patients.
#' @param colors Named character vector of colors for each model. If NULL, uses default colors.
#' @param highlight_ensemble Logical, whether to highlight ensemble predictions in black (default TRUE).
#' @param title Character string for plot title. If NULL, uses default title.
#' @param subtitle Character string for plot subtitle. If NULL, uses default subtitle.
#' @param ncol_facets Integer, number of columns for facetting by patient (default 3).
#' @param legend_position Character, position of legend ("bottom", "right", "top", "left", "none").
#' @param event_label Character string describing the event of interest (default "Event").
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Single model
#' cox_pred <- Predict_CRModel_Cox(cox_model, test_data)
#' plot_cif_curves(cox_pred, model_names = "Cox", event_label = "Relapse")
#' 
#' # Multiple models
#' predictions <- list(
#'   Cox = Predict_CRModel_Cox(cox_model, test_data),
#'   FineGray = Predict_CRModel_FineGray(fg_model, test_data)
#' )
#' plot_cif_curves(predictions, patients_to_plot = 1:5, event_label = "Disease Progression")
#' }
#'
#' @importFrom rlang .data
#' @importFrom grDevices rainbow
#' @export
plot_cif_curves <- function(predictions,
                           model_names = NULL,
                           patients_to_plot = NULL,
                           colors = NULL,
                           highlight_ensemble = TRUE,
                           title = NULL,
                           subtitle = NULL,
                           ncol_facets = 3,
                           legend_position = "bottom",
                           event_label = "Event") {
  
  # Input validation
  if (is.null(predictions)) {
    stop("'predictions' cannot be NULL")
  }
  
  # Handle single prediction object
  if (!is.list(predictions) || (!is.null(predictions$CIFs) && !is.null(predictions$Times))) {
    predictions <- list(Model1 = predictions)
  }
  
  # Validate prediction objects
  for (i in seq_along(predictions)) {
    pred <- predictions[[i]]
    if (is.null(pred$CIFs) || is.null(pred$Times)) {
      stop("Each prediction object must have 'CIFs' and 'Times' components")
    }
    if (!is.matrix(pred$CIFs)) {
      stop("'CIFs' component must be a matrix")
    }
    if (!is.numeric(pred$Times)) {
      stop("'Times' component must be numeric")
    }
  }
  
  # Set model names
  if (is.null(model_names)) {
    if (!is.null(names(predictions))) {
      model_names <- names(predictions)
    } else {
      model_names <- paste0("Model", seq_along(predictions))
      names(predictions) <- model_names
    }
  } else {
    if (length(model_names) != length(predictions)) {
      stop("Length of 'model_names' must match length of 'predictions'")
    }
    names(predictions) <- model_names
  }
  
  # Determine patients to plot (CIFs matrix is times x observations)
  n_patients <- ncol(predictions[[1]]$CIFs)
  if (is.null(patients_to_plot)) {
    patients_to_plot <- seq_len(min(3, n_patients))
  } else {
    if (max(patients_to_plot) > n_patients) {
      stop("'patients_to_plot' contains indices larger than number of patients")
    }
    if (min(patients_to_plot) < 1) {
      stop("'patients_to_plot' must contain positive integers")
    }
  }
  
  # Create plot data
  plot_data <- data.frame()
  
  for (patient in patients_to_plot) {
    for (model_name in model_names) {
      pred <- predictions[[model_name]]
      model_data <- data.frame(
        Time = pred$Times,
        CIF = pred$CIFs[, patient],  # CIFs matrix is [times, observations]
        Model = model_name,
        Patient = paste("Patient", patient),
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, model_data)
    }
  }
  
  # Set colors
  if (is.null(colors)) {
    unique_models <- unique(plot_data$Model)
    colors <- rainbow(length(unique_models))
    names(colors) <- unique_models
    
    # Highlight ensemble if requested
    if (highlight_ensemble && "Ensemble" %in% unique_models) {
      colors["Ensemble"] <- "black"
    }
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$CIF, color = .data$Model)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors)
  
  # Add facetting if multiple patients
  if (length(patients_to_plot) > 1) {
    p <- p + ggplot2::facet_wrap(~Patient, ncol = ncol_facets)
  }
  
  # Set title and subtitle
  if (is.null(title)) {
    title <- paste("Cumulative Incidence Functions (CIF)")
    if (length(model_names) > 1) {
      title <- paste(title, "- Model Comparison")
    }
  }
  
  if (is.null(subtitle)) {
    subtitle <- paste("Event of Interest:", event_label, "|", 
                     "Comparing", length(model_names), "model(s)")
    if (length(patients_to_plot) > 1) {
      subtitle <- paste(subtitle, "across", length(patients_to_plot), "patients")
    }
  }
  
  p <- p + ggplot2::labs(title = title, subtitle = subtitle, 
                          x = "Time", y = paste("Cumulative Incidence of", event_label)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = legend_position,
                   plot.title = ggplot2::element_text(face = "bold"))
  
  # Adjust legend for multiple models
  if (length(model_names) > 4 && legend_position == "bottom") {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(ncol = min(5, length(model_names))))
  }
  
  return(p)
}


#' @title plot_model_performance
#'
#' @description Create performance comparison plots for survival or competing risks models
#'
#' @param performance_df Data frame with columns 'Model', 'Concordance', and optionally 'Brier_Score', 'Prediction_Error'
#' @param metric Character string specifying which metric to plot ("concordance", "brier", "prediction_error")
#' @param highlight_ensemble Logical, whether to highlight ensemble in black (default TRUE)
#' @param title Character string for plot title. If NULL, uses default title.
#' @param flip_coords Logical, whether to flip coordinates (horizontal bars, default TRUE)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' performance_df <- data.frame(
#'   Model = c("Cox", "RF", "Ensemble"),
#'   Concordance = c(0.72, 0.75, 0.78),
#'   Brier_Score = c(0.15, 0.12, 0.10)
#' )
#' 
#' # Concordance plot
#' plot_model_performance(performance_df, metric = "concordance")
#' 
#' # Brier score plot
#' plot_model_performance(performance_df, metric = "brier")
#' }
#'
#' @importFrom rlang .data
#' @export
plot_model_performance <- function(performance_df,
                                  metric = "concordance",
                                  highlight_ensemble = TRUE,
                                  title = NULL,
                                  flip_coords = TRUE) {
  
  # Input validation
  if (!is.data.frame(performance_df)) {
    stop("'performance_df' must be a data frame")
  }
  
  if (!"Model" %in% colnames(performance_df)) {
    stop("'performance_df' must have a 'Model' column")
  }
  

  metric <- match.arg(metric, c("concordance", "brier", "prediction_error", "expected_time_lost"))

  # Check if metric column exists
  metric_col <- switch(metric,
                      "concordance" = "Concordance",
                      "brier" = "Brier_Score", 
                      "prediction_error" = "Prediction_Error",
                      "expected_time_lost" = "Expected_Time_Lost")

  if (!metric_col %in% colnames(performance_df)) {
    stop(paste("'performance_df' must have a", metric_col, "column"))
  }
  
  # Remove rows with missing values for the metric
  performance_df <- performance_df[!is.na(performance_df[[metric_col]]), ]
  
  if (nrow(performance_df) == 0) {
    stop(paste("No valid values found for", metric_col))
  }
  
  # Set up plotting variables
  y_var <- metric_col
  fill_var <- if (highlight_ensemble) performance_df$Model == "Ensemble" else NULL
  
  # Create reordering for y-axis
  if (metric == "brier" || metric == "prediction_error" || metric == "expected_time_lost") {
    # For metrics where lower is better, reverse ordering
    x_order <- reorder(performance_df$Model, -performance_df[[metric_col]])
  } else {
    # For metrics where higher is better
    x_order <- reorder(performance_df$Model, performance_df[[metric_col]])
  }
  
  # Create plot
  p <- ggplot2::ggplot(performance_df, ggplot2::aes(x = x_order, y = .data[[y_var]]))
  
  if (highlight_ensemble) {
    p <- p + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$Model == "Ensemble")) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "steelblue"), guide = "none")
  } else {
    p <- p + ggplot2::geom_bar(stat = "identity", fill = "steelblue")
  }
  
  # Flip coordinates if requested
  if (flip_coords) {
    p <- p + ggplot2::coord_flip()
  }
  
  # Set title and labels

  if (is.null(title)) {
    title <- switch(metric,
                   "concordance" = "Model Performance: Concordance Index",
                   "brier" = "Model Performance: Brier Score", 
                   "prediction_error" = "Model Performance: Prediction Error",
                   "expected_time_lost" = "Model Performance: Expected Time Lost")
  }

  subtitle <- switch(metric,
                    "concordance" = "Higher is better",
                    "brier" = "Lower is better",
                    "prediction_error" = "Lower is better",
                    "expected_time_lost" = "Lower is better")

  if (highlight_ensemble) {
    subtitle <- paste(subtitle, "| Ensemble in BLACK")
  }

  y_label <- switch(metric,
                   "concordance" = "Concordance Index",
                   "brier" = "Brier Score",
                   "prediction_error" = "Prediction Error",
                   "expected_time_lost" = "Expected Time Lost")
  
  p <- p + ggplot2::labs(title = title, subtitle = subtitle,
                         x = "Model", y = y_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  return(p)
}