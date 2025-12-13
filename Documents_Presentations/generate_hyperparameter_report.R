#!/usr/bin/env Rscript
# Generate visualizations and report from hyperparameter test results

library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

# Read results
results <- read.csv("bart_hyperparameter_results.csv", stringsAsFactors = FALSE)

# Create output directory for plots
dir.create("hyperparameter_plots", showWarnings = FALSE)

# ============================================================================
# Plot 1: Execution Time by Sample Size and Configuration
# ============================================================================

p1 <- ggplot(results, aes(x = n, y = total_time, color = config, shape = config)) +
    geom_point(size = 3) +
    geom_line() +
    facet_grid(model ~ p, labeller = label_both) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
        title = "BART Execution Time by Sample Size and Configuration",
        subtitle = "Faceted by model type and number of predictors",
        x = "Sample Size (log scale)",
        y = "Total Time (seconds, log scale)",
        color = "Configuration",
        shape = "Configuration"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("hyperparameter_plots/execution_time_by_sample_size.png", p1, width = 12, height = 8, dpi = 300)

# ============================================================================
# Plot 2: Fit vs Prediction Time
# ============================================================================

results_long <- results %>%
    pivot_longer(cols = c(fit_time, pred_time), names_to = "time_type", values_to = "time")

p2 <- ggplot(results_long, aes(x = config, y = time, fill = time_type)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(model ~ n, labeller = label_both, scales = "free_y") +
    labs(
        title = "BART Fit vs Prediction Time",
        subtitle = "Faceted by model type and sample size",
        x = "Configuration",
        y = "Time (seconds)",
        fill = "Time Type"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave("hyperparameter_plots/fit_vs_pred_time.png", p2, width = 12, height = 8, dpi = 300)

# ============================================================================
# Plot 3: Time per Observation
# ============================================================================

results$time_per_obs <- results$total_time / results$n_test

p3 <- ggplot(results, aes(x = p, y = time_per_obs, color = config, linetype = model)) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(~n, labeller = label_both) +
    labs(
        title = "BART Time per Test Observation",
        subtitle = "By number of predictors and sample size",
        x = "Number of Predictors",
        y = "Time per Observation (seconds)",
        color = "Configuration",
        linetype = "Model Type"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("hyperparameter_plots/time_per_observation.png", p3, width = 12, height = 6, dpi = 300)

# ============================================================================
# Plot 4: Configuration Comparison
# ============================================================================

config_summary <- results %>%
    group_by(model, config) %>%
    summarise(
        mean_time = mean(total_time),
        sd_time = sd(total_time),
        .groups = "drop"
    )

p4 <- ggplot(config_summary, aes(x = config, y = mean_time, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time),
        position = position_dodge(0.9),
        width = 0.2
    ) +
    labs(
        title = "Mean Execution Time by Configuration",
        subtitle = "Error bars show standard deviation across all test scenarios",
        x = "Configuration",
        y = "Mean Total Time (seconds)",
        fill = "Model Type"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("hyperparameter_plots/configuration_comparison.png", p4, width = 10, height = 6, dpi = 300)

# ============================================================================
# Plot 5: Scalability Analysis
# ============================================================================

# Calculate time ratio relative to minimal configuration
results_scaled <- results %>%
    group_by(model, n, p) %>%
    mutate(
        time_ratio = total_time / total_time[config == "Minimal"]
    ) %>%
    ungroup()

p5 <- ggplot(results_scaled, aes(x = n, y = time_ratio, color = config, shape = config)) +
    geom_point(size = 3) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    facet_grid(model ~ p, labeller = label_both) +
    scale_x_log10() +
    labs(
        title = "Execution Time Relative to Minimal Configuration",
        subtitle = "Faceted by model type and number of predictors",
        x = "Sample Size (log scale)",
        y = "Time Ratio (relative to Minimal)",
        color = "Configuration",
        shape = "Configuration"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("hyperparameter_plots/scalability_analysis.png", p5, width = 12, height = 8, dpi = 300)

# ============================================================================
# Generate Summary Tables
# ============================================================================

# Table 1: Overall summary by configuration
summary_by_config <- results %>%
    group_by(model, config) %>%
    summarise(
        n_tests = n(),
        mean_total_time = mean(total_time),
        sd_total_time = sd(total_time),
        mean_fit_time = mean(fit_time),
        mean_pred_time = mean(pred_time),
        .groups = "drop"
    ) %>%
    arrange(model, config)

write.csv(summary_by_config, "hyperparameter_plots/summary_by_configuration.csv", row.names = FALSE)

# Table 2: Summary by sample size
summary_by_n <- results %>%
    group_by(model, n, config) %>%
    summarise(
        mean_total_time = mean(total_time),
        mean_fit_time = mean(fit_time),
        mean_pred_time = mean(pred_time),
        .groups = "drop"
    ) %>%
    arrange(model, n, config)

write.csv(summary_by_n, "hyperparameter_plots/summary_by_sample_size.csv", row.names = FALSE)

# Table 3: Summary by number of predictors
summary_by_p <- results %>%
    group_by(model, p, config) %>%
    summarise(
        mean_total_time = mean(total_time),
        mean_fit_time = mean(fit_time),
        mean_pred_time = mean(pred_time),
        .groups = "drop"
    ) %>%
    arrange(model, p, config)

write.csv(summary_by_p, "hyperparameter_plots/summary_by_predictors.csv", row.names = FALSE)

# ============================================================================
# Print Summary
# ============================================================================

cat("\n", strrep("=", 80), "\n")
cat("VISUALIZATION AND SUMMARY GENERATION COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Generated Files:\n")
cat("  - hyperparameter_plots/execution_time_by_sample_size.png\n")
cat("  - hyperparameter_plots/fit_vs_pred_time.png\n")
cat("  - hyperparameter_plots/time_per_observation.png\n")
cat("  - hyperparameter_plots/configuration_comparison.png\n")
cat("  - hyperparameter_plots/scalability_analysis.png\n")
cat("  - hyperparameter_plots/summary_by_configuration.csv\n")
cat("  - hyperparameter_plots/summary_by_sample_size.csv\n")
cat("  - hyperparameter_plots/summary_by_predictors.csv\n\n")

# Print key findings
cat("KEY FINDINGS:\n")
cat(strrep("-", 80), "\n\n")

# Find fastest and slowest configurations
fastest <- summary_by_config %>%
    group_by(model) %>%
    slice_min(mean_total_time, n = 1)

slowest <- summary_by_config %>%
    group_by(model) %>%
    slice_max(mean_total_time, n = 1)

for (m in unique(results$model)) {
    cat(sprintf("%s BART:\n", m))
    fast <- fastest[fastest$model == m, ]
    slow <- slowest[slowest$model == m, ]

    cat(sprintf(
        "  Fastest: %s (%.2f ± %.2f sec)\n",
        fast$config, fast$mean_total_time, fast$sd_total_time
    ))
    cat(sprintf(
        "  Slowest: %s (%.2f ± %.2f sec)\n",
        slow$config, slow$mean_total_time, slow$sd_total_time
    ))
    cat(sprintf(
        "  Speed-up factor: %.2fx\n\n",
        slow$mean_total_time / fast$mean_total_time
    ))
}

# Scalability analysis
cat("SCALABILITY:\n")
cat(strrep("-", 80), "\n\n")

for (m in unique(results$model)) {
    for (cfg in unique(results$config)) {
        subset_data <- results[results$model == m & results$config == cfg, ]
        if (nrow(subset_data) > 0) {
            # Simple linear model: log(time) ~ log(n) + log(p)
            if (length(unique(subset_data$n)) > 1 && length(unique(subset_data$p)) > 1) {
                lm_model <- lm(log(total_time) ~ log(n) + log(p), data = subset_data)
                coefs <- coef(lm_model)

                cat(sprintf("%s - %s:\n", m, cfg))
                cat(sprintf("  Sample size scaling: O(n^%.2f)\n", coefs["log(n)"]))
                cat(sprintf("  Predictor scaling: O(p^%.2f)\n\n", coefs["log(p)"]))
            }
        }
    }
}

cat("\nAll visualizations and summaries have been generated successfully!\n")
