#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme_minimal theme element_text ggsave
#' @importFrom grDevices rainbow
.organize_predictions_for_plotting <- function(predictions_output,
                                              model_type = "survival",
                                              include_ensemble = TRUE) {

  # Validate input
  if (!is.list(predictions_output)) {
    stop("predictions_output must be a list from PredictSurvModels or PredictCRModels")
  }

  if (is.null(predictions_output$ModelPredictions) || is.null(predictions_output$NewTimes)) {
    stop("predictions_output must contain 'ModelPredictions' and 'NewTimes' components")
  }

  individual_predictions <- predictions_output$ModelPredictions
  common_times <- predictions_output$NewTimes
  ensemble_probs <- predictions_output$NewProbs

  # Determine the probability/CIF component name
  prob_name <- if (tolower(model_type) == "competing_risks" || tolower(model_type) == "cr") {
    "CIFs"
  } else {
    "Probs"
  }

  # Organize individual model predictions
  model_prediction_objects <- lapply(individual_predictions, function(mat) {
    pred_obj <- list(Times = common_times)
    pred_obj[[prob_name]] <- mat
    pred_obj
  })

  # Format model names
  raw_names <- names(individual_predictions)
  formatted_names <- format_model_name(raw_names, model_type = model_type)
  names(model_prediction_objects) <- formatted_names

  # Add ensemble if requested and available
  if (include_ensemble && !is.null(ensemble_probs)) {
    ensemble_obj <- list(Times = common_times)
    ensemble_obj[[prob_name]] <- ensemble_probs
    model_prediction_objects[["Ensemble"]] <- ensemble_obj
  }

  return(model_prediction_objects)
}


#'
#'
#'
#'
#'
#'
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

  # Check if this is output from PredictSurvModels()
  if (!is.null(predictions$ModelPredictions) && !is.null(predictions$NewTimes)) {
    predictions <- .organize_predictions_for_plotting(
      predictions,
      model_type = "survival",
      include_ensemble = highlight_ensemble
    )
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
    if (nrow(pred$Probs) != length(pred$Times)) {
      stop("Each prediction object's 'Times' vector must align with the number of rows in 'Probs'")
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


#'
#'
#'
#'
#'
#'
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

  # Check if this is output from PredictCRModels()
  if (!is.null(predictions$ModelPredictions) && !is.null(predictions$NewTimes)) {
    predictions <- .organize_predictions_for_plotting(
      predictions,
      model_type = "competing_risks",
      include_ensemble = highlight_ensemble
    )
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


#'
#'
#'
#'
#'
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
