#!/usr/bin/env Rscript
# Comprehensive Hyperparameter Performance Testing for BART Models
# Tests survival and competing risks BART across different configurations

library(BART)
library(survival)

# Source the implementation files
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/general_utils.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/data_summary.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_interpolation.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/surv_bart.R")
source("/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/R/cr_bart.R")

# ============================================================================
# Helper Functions
# ============================================================================

#' Generate realistic survival data with predictor effects
generate_survival_data <- function(n, p, censoring_rate = 0.3, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Generate predictors
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("x", 1:p)

    # True coefficients (first half have effects, second half are noise)
    n_signal <- ceiling(p / 2)
    beta <- c(rnorm(n_signal, mean = 0, sd = 0.5), rep(0, p - n_signal))

    # Linear predictor
    eta <- as.vector(X %*% beta)

    # Generate survival times with Weibull distribution
    # Scale parameter depends on predictors
    scale <- exp(eta)
    shape <- 1.5 # Weibull shape parameter

    true_time <- rweibull(n, shape = shape, scale = scale)

    # Generate censoring times
    censor_time <- rexp(n, rate = -log(1 - censoring_rate) / median(true_time))

    # Observed time and event indicator
    time <- pmin(true_time, censor_time)
    event <- as.integer(true_time <= censor_time)

    # Create data frame
    data <- data.frame(time = time, event = event, X)

    # Return data and true parameters
    list(
        data = data,
        beta = beta,
        censoring_rate = mean(event == 0),
        median_time = median(time)
    )
}

#' Generate realistic competing risks data
generate_competing_risks_data <- function(n, p, censoring_rate = 0.2, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Generate predictors
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("x", 1:p)

    # True coefficients for each event type
    n_signal <- ceiling(p / 2)
    beta1 <- c(rnorm(n_signal, mean = 0.3, sd = 0.4), rep(0, p - n_signal))
    beta2 <- c(rnorm(n_signal, mean = -0.3, sd = 0.4), rep(0, p - n_signal))

    # Linear predictors
    eta1 <- as.vector(X %*% beta1)
    eta2 <- as.vector(X %*% beta2)

    # Generate event times for each cause
    time1 <- rweibull(n, shape = 1.5, scale = exp(eta1))
    time2 <- rweibull(n, shape = 1.3, scale = exp(eta2))

    # Generate censoring times
    median_time <- median(c(time1, time2))
    censor_time <- rexp(n, rate = -log(1 - censoring_rate) / median_time)

    # Determine observed event
    time <- pmin(time1, time2, censor_time)
    event <- ifelse(time == censor_time, 0,
        ifelse(time == time1, 1, 2)
    )

    # Create data frame
    data <- data.frame(time = time, event = event, X)

    list(
        data = data,
        beta1 = beta1,
        beta2 = beta2,
        censoring_rate = mean(event == 0),
        event1_rate = mean(event == 1),
        event2_rate = mean(event == 2)
    )
}

#' Calculate prediction correlation between two models
calc_prediction_correlation <- function(pred1, pred2) {
    # Handle both survival (Probs) and competing risks (CIFs)
    if (!is.null(pred1$Probs) && !is.null(pred2$Probs)) {
        # Survival predictions
        cor(as.vector(pred1$Probs), as.vector(pred2$Probs), use = "complete.obs")
    } else if (!is.null(pred1$CIFs) && !is.null(pred2$CIFs)) {
        # Competing risks predictions
        cor(as.vector(pred1$CIFs), as.vector(pred2$CIFs), use = "complete.obs")
    } else {
        NA
    }
}

# ============================================================================
# Test Configurations
# ============================================================================

# Sample sizes to test
sample_sizes <- c(100, 500, 1000)

# Number of predictors to test
n_predictors <- c(5, 10, 20)

# Hyperparameter configurations for survival BART
surv_configs <- list(
    minimal = list(K = 5, ntree = 20, ndpost = 100, nskip = 25, keepevery = 2, label = "Minimal"),
    default = list(K = 8, ntree = 50, ndpost = 200, nskip = 50, keepevery = 2, label = "Default"),
    enhanced = list(K = 10, ntree = 100, ndpost = 400, nskip = 100, keepevery = 2, label = "Enhanced")
)

# Hyperparameter configurations for competing risks BART
cr_configs <- list(
    minimal = list(K = 8, ntree = 50, ndpost = 500, nskip = 50, keepevery = 5, label = "Minimal"),
    default = list(K = 10, ntree = 100, ndpost = 1000, nskip = 100, keepevery = 10, label = "Default"),
    enhanced = list(K = 15, ntree = 150, ndpost = 2000, nskip = 200, keepevery = 10, label = "Enhanced")
)

# ============================================================================
# Run Survival BART Tests
# ============================================================================

cat(strrep("=", 80), "\n")
cat("SURVIVAL BART HYPERPARAMETER TESTING\n")
cat(strrep("=", 80), "\n\n")

surv_results <- list()
result_idx <- 1

for (n in sample_sizes) {
    for (p in n_predictors) {
        cat("\n", strrep("=", 80), "\n")
        cat(sprintf("Sample Size: %d, Predictors: %d\n", n, p))
        cat(strrep("=", 80), "\n")

        # Generate data
        sim_data <- generate_survival_data(n = n, p = p, censoring_rate = 0.3, seed = 123)
        train_data <- sim_data$data[1:floor(0.7 * n), ]
        test_data <- sim_data$data[(floor(0.7 * n) + 1):n, ]

        cat(sprintf("Censoring rate: %.1f%%\n", sim_data$censoring_rate * 100))
        cat(sprintf("Median time: %.2f\n\n", sim_data$median_time))

        # Store predictions for correlation analysis
        predictions <- list()

        for (config_name in names(surv_configs)) {
            config <- surv_configs[[config_name]]
            cat(sprintf("Testing %s configuration...\n", config$label))

            # Time the fitting
            start_time <- Sys.time()

            model <- tryCatch(
                {
                    SurvModel_BART(
                        data = train_data,
                        expvars = paste0("x", 1:p),
                        timevar = "time",
                        eventvar = "event",
                        K = config$K,
                        ntree = config$ntree
                    )
                },
                error = function(e) {
                    cat("  ERROR:", e$message, "\n")
                    NULL
                }
            )

            if (is.null(model)) {
                cat("  Skipping due to fitting error\n\n")
                next
            }

            fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

            # Time the prediction
            start_time <- Sys.time()

            pred <- tryCatch(
                {
                    Predict_SurvModel_BART(
                        modelout = model,
                        newdata = test_data,
                        new_times = quantile(train_data$time, probs = seq(0.1, 0.9, by = 0.1))
                    )
                },
                error = function(e) {
                    cat("  PREDICTION ERROR:", e$message, "\n")
                    NULL
                }
            )

            if (is.null(pred)) {
                cat("  Skipping due to prediction error\n\n")
                next
            }

            pred_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            total_time <- fit_time + pred_time

            # Store predictions
            predictions[[config_name]] <- pred

            # Store results
            surv_results[[result_idx]] <- list(
                n = n,
                p = p,
                config = config$label,
                K = config$K,
                ntree = config$ntree,
                ndpost = config$ndpost,
                fit_time = fit_time,
                pred_time = pred_time,
                total_time = total_time,
                n_test = nrow(test_data),
                n_times = length(pred$Times)
            )
            result_idx <- result_idx + 1

            cat(sprintf("  Fit time: %.2f sec\n", fit_time))
            cat(sprintf("  Prediction time: %.2f sec\n", pred_time))
            cat(sprintf("  Total time: %.2f sec\n\n", total_time))
        }

        # Calculate correlations between configurations
        if (length(predictions) >= 2) {
            cat("Prediction Correlations:\n")
            config_names <- names(predictions)
            for (i in 1:(length(config_names) - 1)) {
                for (j in (i + 1):length(config_names)) {
                    corr <- calc_prediction_correlation(
                        predictions[[config_names[i]]],
                        predictions[[config_names[j]]]
                    )
                    cat(sprintf(
                        "  %s vs %s: %.4f\n",
                        surv_configs[[config_names[i]]]$label,
                        surv_configs[[config_names[j]]]$label,
                        corr
                    ))
                }
            }
            cat("\n")
        }
    }
}

# ============================================================================
# Run Competing Risks BART Tests
# ============================================================================

cat("\n\n", strrep("=", 80), "\n")
cat("COMPETING RISKS BART HYPERPARAMETER TESTING\n")
cat(strrep("=", 80), "\n\n")

cr_results <- list()
result_idx <- 1

for (n in sample_sizes) {
    for (p in n_predictors) {
        cat("\n", strrep("=", 80), "\n")
        cat(sprintf("Sample Size: %d, Predictors: %d\n", n, p))
        cat(strrep("=", 80), "\n")

        # Generate data
        sim_data <- generate_competing_risks_data(n = n, p = p, censoring_rate = 0.2, seed = 456)
        train_data <- sim_data$data[1:floor(0.7 * n), ]
        test_data <- sim_data$data[(floor(0.7 * n) + 1):n, ]

        cat(sprintf("Censoring rate: %.1f%%\n", sim_data$censoring_rate * 100))
        cat(sprintf("Event 1 rate: %.1f%%\n", sim_data$event1_rate * 100))
        cat(sprintf("Event 2 rate: %.1f%%\n\n", sim_data$event2_rate * 100))

        # Store predictions for correlation analysis
        predictions <- list()

        for (config_name in names(cr_configs)) {
            config <- cr_configs[[config_name]]
            cat(sprintf("Testing %s configuration...\n", config$label))

            # Time the fitting
            start_time <- Sys.time()

            model <- tryCatch(
                {
                    CRModel_BART(
                        data = train_data,
                        expvars = paste0("x", 1:p),
                        timevar = "time",
                        eventvar = "event",
                        event_codes = 1,
                        K = config$K,
                        ntree = config$ntree,
                        ndpost = config$ndpost,
                        nskip = config$nskip,
                        keepevery = config$keepevery,
                        verbose = FALSE
                    )
                },
                error = function(e) {
                    cat("  ERROR:", e$message, "\n")
                    NULL
                }
            )

            if (is.null(model)) {
                cat("  Skipping due to fitting error\n\n")
                next
            }

            fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

            # Time the prediction
            start_time <- Sys.time()

            pred <- tryCatch(
                {
                    Predict_CRModel_BART(
                        modelout = model,
                        newdata = test_data,
                        new_times = quantile(train_data$time, probs = seq(0.1, 0.9, by = 0.1))
                    )
                },
                error = function(e) {
                    cat("  PREDICTION ERROR:", e$message, "\n")
                    NULL
                }
            )

            if (is.null(pred)) {
                cat("  Skipping due to prediction error\n\n")
                next
            }

            pred_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            total_time <- fit_time + pred_time

            # Store predictions
            predictions[[config_name]] <- pred

            # Store results
            cr_results[[result_idx]] <- list(
                n = n,
                p = p,
                config = config$label,
                K = config$K,
                ntree = config$ntree,
                ndpost = config$ndpost,
                fit_time = fit_time,
                pred_time = pred_time,
                total_time = total_time,
                n_test = nrow(test_data),
                n_times = length(pred$Times)
            )
            result_idx <- result_idx + 1

            cat(sprintf("  Fit time: %.2f sec\n", fit_time))
            cat(sprintf("  Prediction time: %.2f sec\n", pred_time))
            cat(sprintf("  Total time: %.2f sec\n\n", total_time))
        }

        # Calculate correlations between configurations
        if (length(predictions) >= 2) {
            cat("Prediction Correlations:\n")
            config_names <- names(predictions)
            for (i in 1:(length(config_names) - 1)) {
                for (j in (i + 1):length(config_names)) {
                    corr <- calc_prediction_correlation(
                        predictions[[config_names[i]]],
                        predictions[[config_names[j]]]
                    )
                    cat(sprintf(
                        "  %s vs %s: %.4f\n",
                        cr_configs[[config_names[i]]]$label,
                        cr_configs[[config_names[j]]]$label,
                        corr
                    ))
                }
            }
            cat("\n")
        }
    }
}

# ============================================================================
# Save Results
# ============================================================================

# Convert to data frames
surv_df <- do.call(rbind, lapply(surv_results, function(x) {
    data.frame(
        model = "Survival",
        n = x$n,
        p = x$p,
        config = x$config,
        K = x$K,
        ntree = x$ntree,
        ndpost = x$ndpost,
        fit_time = x$fit_time,
        pred_time = x$pred_time,
        total_time = x$total_time,
        n_test = x$n_test,
        n_times = x$n_times,
        stringsAsFactors = FALSE
    )
}))

cr_df <- do.call(rbind, lapply(cr_results, function(x) {
    data.frame(
        model = "Competing Risks",
        n = x$n,
        p = x$p,
        config = x$config,
        K = x$K,
        ntree = x$ntree,
        ndpost = x$ndpost,
        fit_time = x$fit_time,
        pred_time = x$pred_time,
        total_time = x$total_time,
        n_test = x$n_test,
        n_times = x$n_times,
        stringsAsFactors = FALSE
    )
}))

# Combine results
all_results <- rbind(surv_df, cr_df)

# Save to CSV
output_file <- "/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event/Documents_Presentations/bart_hyperparameter_results.csv"
write.csv(all_results, output_file, row.names = FALSE)

cat("\n\n", strrep("=", 80), "\n")
cat("RESULTS SAVED TO:", output_file, "\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n", strrep("=", 80), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 80), "\n\n")

cat("SURVIVAL BART:\n")
cat(strrep("-", 80), "\n")
for (config in unique(surv_df$config)) {
    subset_data <- surv_df[surv_df$config == config, ]
    cat(sprintf("\n%s Configuration:\n", config))
    cat(sprintf(
        "  Mean total time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$total_time), sd(subset_data$total_time)
    ))
    cat(sprintf(
        "  Mean fit time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$fit_time), sd(subset_data$fit_time)
    ))
    cat(sprintf(
        "  Mean pred time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$pred_time), sd(subset_data$pred_time)
    ))
}

cat("\n\nCOMPETING RISKS BART:\n")
cat(strrep("-", 80), "\n")
for (config in unique(cr_df$config)) {
    subset_data <- cr_df[cr_df$config == config, ]
    cat(sprintf("\n%s Configuration:\n", config))
    cat(sprintf(
        "  Mean total time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$total_time), sd(subset_data$total_time)
    ))
    cat(sprintf(
        "  Mean fit time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$fit_time), sd(subset_data$fit_time)
    ))
    cat(sprintf(
        "  Mean pred time: %.2f sec (SD: %.2f)\n",
        mean(subset_data$pred_time), sd(subset_data$pred_time)
    ))
}

cat("\n\nTesting complete!\n")
