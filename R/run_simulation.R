#' @title run_survival_simulation
#'
#' @description Run a complete simulation study for survival models with multiple replicates.
#' Generates simulated data, trains models, makes predictions, and evaluates performance.
#'
#' @param n_replicates Number of simulation replicates
#' @param n_train Sample size for training data
#' @param n_test Sample size for test data
#' @param p Number of predictors
#' @param data_params List of parameters for data generation (passed to generate_survival_data)
#' @param models_to_run Character vector of model names to run (NULL for all)
#' @param model_params List of model-specific parameters
#' @param eval_times Time points for evaluation
#' @param seeds Optional vector of random seeds (one per replicate)
#' @param generate_plots Logical, whether to generate and save plots (default TRUE)
#' @param generate_report Logical, whether to generate and save a summary report (default TRUE)
#'
#' @return A list containing aggregated results across all replicates:
#'   \item{calibration_results}{Calibration metrics aggregated across replicates}
#'   \item{accuracy_results}{Accuracy metrics aggregated across replicates}
#'   \item{replicate_results}{Individual results from each replicate (if save_intermediate=TRUE)}
#'   \item{simulation_params}{Parameters used for the simulation}
#'   \item{summary}{Overall summary statistics}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run simulation with 10 replicates
#' sim_results <- run_survival_simulation(
#'   n_replicates = 10,
#'   n_train = 500,
#'   n_test = 200,
#'   p = 5,
#'   data_params = list(
#'     baseline_hazard_type = "weibull",
#'     baseline_params = list(shape = 1.5, scale = 2),
#'     linear_effects = c(0.5, -0.3, 0.2, -0.1, 0.4),
#'     censoring_rate = 0.3
#'   ),
#'   models_to_run = c("Cox", "gbm", "xgboost"),
#'   eval_times = seq(0.1, 3, by = 0.2)
#' )
#' }
#'
#' @export
run_survival_simulation <- function(n_replicates = 10,
                                    n_train = 1000,
                                    n_test = 500,
                                    p = 10,
                                    data_params = list(),
                                    models_to_run = NULL,
                                    model_params = list(),
                                    eval_times = seq(0.1, 3, by = 0.2),
                                    seeds = NULL,
                                    save_intermediate = FALSE,
                                    results_dir = "simulation_results",
                                    generate_plots = TRUE,
                                    generate_report = TRUE,
                                    verbose = TRUE) {

  # Set default data parameters
  default_data_params <- list(
    baseline_hazard_type = "weibull",
    baseline_params = list(shape = 1.5, scale = 2),
    linear_effects = rep(0.1, p),
    nonlinear_effects = NULL,
    censoring_rate = 0.3
  )

  # Merge with user-provided parameters
  data_params <- modifyList(default_data_params, data_params)

  # Set seeds if not provided
  if (is.null(seeds)) {
    seeds <- 1:n_replicates
  } else if (length(seeds) != n_replicates) {
    stop("Length of seeds must equal n_replicates")
  }

  # Initialize results storage
  all_calibration <- list()
  all_accuracy <- list()
  replicate_results <- list()

  # Create results directory if generating plots or reports
  if ((generate_plots || generate_report) && !dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  # Run simulation replicates
  for (rep in 1:n_replicates) {
    if (verbose) cat("\n=== Replicate", rep, "of", n_replicates, "===\n")

    tryCatch({
      # Set seed for this replicate
      set.seed(seeds[rep])

      # Generate training data
      if (verbose) cat("Generating training data...\n")
      train_data <- do.call(generate_survival_data, c(
        list(n = n_train, p = p),
        data_params
      ))

      # Generate test data (same data generating process)
      if (verbose) cat("Generating test data...\n")
      test_data <- do.call(generate_survival_data, c(
        list(n = n_test, p = p),
        data_params
      ))

      # Train models
      if (verbose) cat("Training models...\n")
      model_results <- run_survival_models(
        data = train_data$data,
        expvars = paste0("x", 1:p),
        timevar = "time",
        eventvar = "event",
        models_to_run = models_to_run,
        model_params = model_params,
        verbose = verbose
      )

      # Make predictions
      if (verbose) cat("Making predictions...\n")
      predictions <- predict_survival_models(
        model_results = model_results,
        newdata = test_data$data,
        newtimes = eval_times
      )

      # Evaluate calibration
      if (verbose) cat("Evaluating calibration...\n")
      calib_results <- evaluate_calibration(
        predictions = predictions,
        sim_data = test_data,  # Use test data's true functions
        test_data = test_data$data,
        eval_times = eval_times
      )

      # Evaluate accuracy
      if (verbose) cat("Evaluating accuracy...\n")
      acc_results <- evaluate_accuracy(
        predictions = predictions,
        test_data = test_data$data,
        timevar = "time",
        eventvar = "event",
        eval_times = eval_times
      )

      # Store results
      all_calibration[[rep]] <- calib_results
      all_accuracy[[rep]] <- acc_results

      if (save_intermediate) {
        replicate_results[[rep]] <- list(
          train_data = train_data,
          test_data = test_data,
          model_results = model_results,
          predictions = predictions,
          calibration = calib_results,
          accuracy = acc_results,
          seed = seeds[rep]
        )

        # Save to file
        saveRDS(replicate_results[[rep]],
                file = file.path(results_dir, paste0("replicate_", rep, ".rds")))
      }

      if (verbose) cat("Replicate", rep, "completed successfully\n")

    }, error = function(e) {
      warning("Replicate ", rep, " failed: ", e$message)
      all_calibration[[rep]] <- list(error = e$message)
      all_accuracy[[rep]] <- list(error = e$message)
    })
  }

  # Aggregate results across replicates
  if (verbose) cat("\nAggregating results across replicates...\n")

  aggregated_results <- aggregate_simulation_results(
    all_calibration = all_calibration,
    all_accuracy = all_accuracy,
    n_replicates = n_replicates
  )

  # Create summary
  summary_stats <- create_simulation_summary(aggregated_results)

  # Final results object
  results <- list(
    calibration_results = aggregated_results$calibration,
    accuracy_results = aggregated_results$accuracy,
    replicate_results = if(save_intermediate) replicate_results else NULL,
    simulation_params = list(
      n_replicates = n_replicates,
      n_train = n_train,
      n_test = n_test,
      p = p,
      data_params = data_params,
      models_to_run = models_to_run,
      eval_times = eval_times,
      seeds = seeds
    ),
    summary = summary_stats
  )

  if (verbose) {
    cat("\nSimulation completed!\n")
    cat("Models evaluated:", paste(summary_stats$models_evaluated, collapse = ", "), "\n")
    cat("Best calibration model (by MAE):", summary_stats$best_calibration_model$by_mae, "\n")
    cat("Best accuracy model (by concordance):", summary_stats$best_accuracy_model$by_integrated_concordance, "\n")
    cat("Overall best model:", summary_stats$best_overall_model, "\n")
  }

  # Generate plots and reports if requested
  if (generate_plots || generate_report) {
    if (verbose) cat("\nGenerating plots and reports...\n")

    # Generate plots
    if (generate_plots) {
      tryCatch({
        generate_simulation_plots(results, results_dir, verbose = verbose)
      }, error = function(e) {
        warning("Failed to generate plots: ", e$message)
      })
    }

    # Generate report
    if (generate_report) {
      tryCatch({
        generate_simulation_report(results, results_dir, verbose = verbose)
      }, error = function(e) {
        warning("Failed to generate report: ", e$message)
      })
    }
  }

  return(results)
}


#' @title aggregate_simulation_results
#'
#' @description Aggregate calibration and accuracy results across simulation replicates
#'
#' @param all_calibration List of calibration results from each replicate
#' @param all_accuracy List of accuracy results from each replicate
#' @param n_replicates Total number of replicates
#'
#' @return List with aggregated calibration and accuracy results
#'
#' @noRd
aggregate_simulation_results <- function(all_calibration, all_accuracy, n_replicates) {

  # Get successful replicates - more robust check
  successful_reps <- which(sapply(all_calibration, function(x) {
    is.list(x) && !is.null(x$model_calibration) && is.list(x$model_calibration)
  }))
  n_successful <- length(successful_reps)

  if (n_successful == 0) {
    stop("No successful replicates to aggregate")
  }

  # Get model names from first successful replicate
  first_success <- successful_reps[1]
  model_names <- names(all_calibration[[first_success]]$model_calibration)

  # Aggregate calibration results
  calib_aggregated <- list()
  for (model in model_names) {
    # Extract metrics across replicates, handling failed models
    mae_values <- sapply(successful_reps, function(rep) {
      result <- all_calibration[[rep]]$model_calibration[[model]]
      if (is.list(result) && !is.null(result$mae)) result$mae else NA
    })

    rmse_values <- sapply(successful_reps, function(rep) {
      result <- all_calibration[[rep]]$model_calibration[[model]]
      if (is.list(result) && !is.null(result$rmse)) result$rmse else NA
    })

    ice_values <- sapply(successful_reps, function(rep) {
      result <- all_calibration[[rep]]$model_calibration[[model]]
      if (is.list(result) && !is.null(result$integrated_calibration_error)) result$integrated_calibration_error else NA
    })

    calib_aggregated[[model]] <- list(
      mae = list(
        mean = mean(mae_values, na.rm = TRUE),
        sd = sd(mae_values, na.rm = TRUE),
        median = median(mae_values, na.rm = TRUE),
        values = mae_values
      ),
      rmse = list(
        mean = mean(rmse_values, na.rm = TRUE),
        sd = sd(rmse_values, na.rm = TRUE),
        median = median(rmse_values, na.rm = TRUE),
        values = rmse_values
      ),
      integrated_calibration_error = list(
        mean = mean(ice_values, na.rm = TRUE),
        sd = sd(ice_values, na.rm = TRUE),
        median = median(ice_values, na.rm = TRUE),
        values = ice_values
      )
    )
  }

  # Aggregate accuracy results
  acc_aggregated <- list()
  for (model in model_names) {
    # Extract integrated concordance, handling failed models
    ic_values <- sapply(successful_reps, function(rep) {
      result <- all_accuracy[[rep]]$model_accuracy[[model]]
      if (is.list(result) && !is.null(result$integrated_concordance) && !is.na(result$integrated_concordance)) {
        result$integrated_concordance
      } else {
        NA
      }
    })

    # Extract integrated Brier score, handling failed models
    ib_values <- sapply(successful_reps, function(rep) {
      result <- all_accuracy[[rep]]$model_accuracy[[model]]
      if (is.list(result) && !is.null(result$integrated_brier_score) && !is.na(result$integrated_brier_score)) {
        result$integrated_brier_score
      } else {
        NA
      }
    })

    acc_aggregated[[model]] <- list(
      integrated_concordance = list(
        mean = mean(ic_values, na.rm = TRUE),
        sd = if (all(is.na(ic_values))) NA else sd(ic_values, na.rm = TRUE),
        median = median(ic_values, na.rm = TRUE),
        values = ic_values
      ),
      integrated_brier_score = list(
        mean = if (all(is.na(ib_values))) NA else mean(ib_values, na.rm = TRUE),
        sd = if (all(is.na(ib_values)) || sum(!is.na(ib_values)) < 2) NA else sd(ib_values, na.rm = TRUE),
        median = if (all(is.na(ib_values))) NA else median(ib_values, na.rm = TRUE),
        values = ib_values
      )
    )
  }

  return(list(
    calibration = calib_aggregated,
    accuracy = acc_aggregated,
    n_successful_replicates = n_successful,
    total_replicates = n_replicates
  ))
}


#' @title create_simulation_summary
#'
#' @description Create a summary of simulation results
#'
#' @param aggregated_results Output from aggregate_simulation_results
#'
#' @return List with summary statistics
#'
#' @noRd
create_simulation_summary <- function(aggregated_results) {

  model_names <- names(aggregated_results$calibration)

  # Best models by different metrics
  mae_means <- sapply(aggregated_results$calibration, function(x) x$mae$mean)
  rmse_means <- sapply(aggregated_results$calibration, function(x) x$rmse$mean)
  ice_means <- sapply(aggregated_results$calibration, function(x) x$integrated_calibration_error$mean)

  ic_means <- sapply(aggregated_results$accuracy, function(x) x$integrated_concordance$mean)
  ib_means <- sapply(aggregated_results$accuracy, function(x) x$integrated_brier_score$mean)

  best_calibration_mae <- names(which.min(mae_means))
  best_calibration_rmse <- names(which.min(rmse_means))
  best_calibration_ice <- names(which.min(ice_means))

  best_accuracy_ic <- names(which.max(ic_means))
  best_accuracy_ib <- names(which.min(ib_means))

  # Overall best (balancing calibration and accuracy)
  # Simple ranking: average of standardized scores
  standardize <- function(x) (x - min(x)) / (max(x) - min(x))

  calib_score <- (standardize(mae_means) + standardize(rmse_means) + standardize(ice_means)) / 3
  acc_score <- (1 - standardize(ic_means)) + standardize(ib_means)  # Lower is better for Brier, higher for C
  overall_score <- (calib_score + acc_score) / 2

  best_overall <- names(which.min(overall_score))

  return(list(
    models_evaluated = model_names,
    best_calibration_model = list(
      by_mae = best_calibration_mae,
      by_rmse = best_calibration_rmse,
      by_integrated_calibration_error = best_calibration_ice
    ),
    best_accuracy_model = list(
      by_integrated_concordance = best_accuracy_ic,
      by_integrated_brier = best_accuracy_ib
    ),
    best_overall_model = best_overall,
    n_successful_replicates = aggregated_results$n_successful_replicates,
    total_replicates = aggregated_results$total_replicates,
    success_rate = aggregated_results$n_successful_replicates / aggregated_results$total_replicates
  ))
}


#' @title generate_simulation_plots
#'
#' @description Generate and save plots from simulation results
#'
#' @param results Output from run_survival_simulation
#' @param results_dir Directory to save plots
#' @param verbose Logical, whether to print progress
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme_minimal theme element_text ggsave
#' @importFrom stats reorder
#'
#' @noRd
generate_simulation_plots <- function(results, results_dir, verbose = TRUE) {

  # Create plots subdirectory
  plots_dir <- file.path(results_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  # Extract aggregated results
  calib_results <- results$calibration_results
  acc_results <- results$accuracy_results

  # 1. Calibration performance plots
  if (verbose) cat("Generating calibration plots...\n")

  # MAE comparison - create data frame with expected column names
  mae_df <- data.frame(
    Model = names(calib_results),
    Prediction_Error = sapply(calib_results, function(x) x$mae$mean),
    stringsAsFactors = FALSE
  )
  mae_df <- mae_df[order(mae_df$Prediction_Error), ]

  mae_plot <- plot_model_performance(mae_df, metric = "prediction_error",
                                    title = "Calibration Performance: Mean Absolute Error")
  ggsave(file.path(plots_dir, "calibration_mae.png"), mae_plot, width = 8, height = 6, dpi = 300)

  # RMSE comparison
  rmse_df <- data.frame(
    Model = names(calib_results),
    Prediction_Error = sapply(calib_results, function(x) x$rmse$mean),
    stringsAsFactors = FALSE
  )
  rmse_df <- rmse_df[order(rmse_df$Prediction_Error), ]

  rmse_plot <- plot_model_performance(rmse_df, metric = "prediction_error",
                                     title = "Calibration Performance: Root Mean Square Error")
  ggsave(file.path(plots_dir, "calibration_rmse.png"), rmse_plot, width = 8, height = 6, dpi = 300)

  # 2. Accuracy performance plots
  if (verbose) cat("Generating accuracy plots...\n")

  # Concordance comparison
  conc_df <- data.frame(
    Model = names(acc_results),
    Concordance = sapply(acc_results, function(x) x$integrated_concordance$mean),
    stringsAsFactors = FALSE
  )
  conc_df <- conc_df[order(conc_df$Concordance, decreasing = TRUE), ]

  conc_plot <- plot_model_performance(conc_df, metric = "concordance",
                                     title = "Accuracy Performance: Integrated Concordance")
  ggsave(file.path(plots_dir, "accuracy_concordance.png"), conc_plot, width = 8, height = 6, dpi = 300)

  # Brier score comparison (only if available)
  brier_values <- sapply(acc_results, function(x) x$integrated_brier_score$mean)
  if (!all(is.na(brier_values))) {
    brier_df <- data.frame(
      Model = names(acc_results),
      Brier_Score = brier_values,
      stringsAsFactors = FALSE
    )
    brier_df <- brier_df[!is.na(brier_df$Brier_Score), ]
    brier_df <- brier_df[order(brier_df$Brier_Score), ]

    brier_plot <- plot_model_performance(brier_df, metric = "brier",
                                        title = "Accuracy Performance: Integrated Brier Score")
    ggsave(file.path(plots_dir, "accuracy_brier.png"), brier_plot, width = 8, height = 6, dpi = 300)
  }

  # 3. Combined performance plot
  if (verbose) cat("Generating combined performance plot...\n")

  combined_df <- data.frame(
    Model = names(calib_results),
    MAE = sapply(calib_results, function(x) x$mae$mean),
    Concordance = sapply(acc_results, function(x) x$integrated_concordance$mean)
  )

  # Create scatter plot
  combined_plot <- ggplot(combined_df, aes(x = MAE, y = Concordance, label = Model)) +
    geom_point(size = 3, color = "blue") +
    geom_text(vjust = -1, size = 3) +
    labs(title = "Model Performance: Calibration vs Accuracy",
         x = "Mean Absolute Error (lower is better)",
         y = "Integrated Concordance (higher is better)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(plots_dir, "combined_performance.png"), combined_plot, width = 8, height = 6, dpi = 300)

  if (verbose) cat("Plots saved to:", plots_dir, "\n")
}


#' @title generate_simulation_report
#'
#' @description Generate a summary report from simulation results
#'
#' @param results Output from run_survival_simulation
#' @param results_dir Directory to save report
#' @param verbose Logical, whether to print progress
#'
#' @noRd
generate_simulation_report <- function(results, results_dir, verbose = TRUE) {

  if (verbose) cat("Generating simulation report...\n")

  # Extract results
  summary_stats <- results$summary
  calib_results <- results$calibration_results
  acc_results <- results$accuracy_results
  sim_params <- results$simulation_params

  # Create report content
  report_content <- c(
    "# Survival Model Simulation Study Report",
    "",
    "## Overview",
    sprintf("This report summarizes the results of a simulation study comparing %d survival models.", length(calib_results)),
    sprintf("The simulation was run with %d replicates, training on %d observations and testing on %d observations.",
            sim_params$n_replicates, sim_params$n_train, sim_params$n_test),
    "",
    "## Simulation Parameters",
    sprintf("- **Models tested**: %s", paste(summary_stats$models_evaluated, collapse = ", ")),
    sprintf("- **Number of replicates**: %d", sim_params$n_replicates),
    sprintf("- **Training sample size**: %d", sim_params$n_train),
    sprintf("- **Test sample size**: %d", sim_params$n_test),
    sprintf("- **Number of predictors**: %d", sim_params$p),
    sprintf("- **Baseline hazard**: %s", sim_params$data_params$baseline_hazard_type),
    sprintf("- **Censoring rate**: %.1f%%", sim_params$data_params$censoring_rate * 100),
    "",
    "## Best Performing Models",
    "",
    "### Calibration (Prediction Accuracy)",
    sprintf("- **Best MAE**: %s (%.4f +/- %.4f)",
            summary_stats$best_calibration_model$by_mae,
            calib_results[[summary_stats$best_calibration_model$by_mae]]$mae$mean,
            calib_results[[summary_stats$best_calibration_model$by_mae]]$mae$sd),
    sprintf("- **Best RMSE**: %s (%.4f +/- %.4f)",
            summary_stats$best_calibration_model$by_rmse,
            calib_results[[summary_stats$best_calibration_model$by_rmse]]$rmse$mean,
            calib_results[[summary_stats$best_calibration_model$by_rmse]]$rmse$sd),
    "",
    "### Accuracy (Discrimination)",
    if (!is.null(summary_stats$best_accuracy_model$by_integrated_concordance) &&
        summary_stats$best_accuracy_model$by_integrated_concordance != "") {
      sprintf("- **Best concordance**: %s (%.4f +/- %.4f)",
              summary_stats$best_accuracy_model$by_integrated_concordance,
              acc_results[[summary_stats$best_accuracy_model$by_integrated_concordance]]$integrated_concordance$mean,
              acc_results[[summary_stats$best_accuracy_model$by_integrated_concordance]]$integrated_concordance$sd)
    } else {
      "- **Best concordance**: Not available (calculation failed)"
    },
    "",
    "## Detailed Results",
    "",
    "### Calibration Metrics",
    "| Model | MAE | RMSE | Integrated CE |",
    "|-------|-----|------|---------------|"
  )

  # Add calibration table rows
  for (model in names(calib_results)) {
    mae <- sprintf("%.4f +/- %.4f", calib_results[[model]]$mae$mean, calib_results[[model]]$mae$sd)
    rmse <- sprintf("%.4f +/- %.4f", calib_results[[model]]$rmse$mean, calib_results[[model]]$rmse$sd)
    ice <- sprintf("%.4f +/- %.4f", calib_results[[model]]$integrated_calibration_error$mean,
                   calib_results[[model]]$integrated_calibration_error$sd)
    report_content <- c(report_content, sprintf("| %s | %s | %s | %s |", model, mae, rmse, ice))
  }

  report_content <- c(report_content, "",
                     "### Accuracy Metrics",
                     "| Model | Integrated Concordance | Integrated Brier |",
                     "|-------|----------------------|------------------|")

  # Add accuracy table rows
  for (model in names(acc_results)) {
    conc <- if (is.na(acc_results[[model]]$integrated_concordance$mean)) {
      "N/A"
    } else {
      sprintf("%.4f +/- %.4f", acc_results[[model]]$integrated_concordance$mean,
              acc_results[[model]]$integrated_concordance$sd)
    }
    brier <- if (is.na(acc_results[[model]]$integrated_brier_score$mean)) {
      "N/A"
    } else {
      sprintf("%.4f +/- %.4f", acc_results[[model]]$integrated_brier_score$mean,
              acc_results[[model]]$integrated_brier_score$sd)
    }
    report_content <- c(report_content, sprintf("| %s | %s | %s |", model, conc, brier))
  }

  report_content <- c(report_content, "",
                     "## Conclusions",
                     sprintf("The simulation study evaluated %d survival models across %d replicates.",
                            length(calib_results), sim_params$n_replicates),
                     "Results show that different models excel in different aspects:",
                     "- Calibration performance measures how well predicted survival curves match true curves",
                     "- Accuracy performance measures how well models discriminate between subjects",
                     "",
                     sprintf("**Overall best model**: %s", summary_stats$best_overall_model),
                     "",
                     "---",
                     sprintf("*Report generated on %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  # Write report to file
  report_file <- file.path(results_dir, "simulation_report.md")
  writeLines(report_content, report_file)

  if (verbose) cat("Report saved to:", report_file, "\n")
}