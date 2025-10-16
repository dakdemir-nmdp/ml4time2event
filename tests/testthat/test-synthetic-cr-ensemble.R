## ============================================================================
## Test file for updated ensemble methods for competing risks models
## ============================================================================

library(testthat)
library(ml4time2event)

test_that("ensemble methods work with competing risks synthetic data", {
  # For CIF testing
  eval_times <- seq(0.1, 1, by=0.2)
  n_times <- length(eval_times)
  
  # Create test data
  set.seed(123)
  n_obs <- 20
  test_data <- data.frame(
    x1 = rnorm(n_obs),
    x2 = rnorm(n_obs),
    group = sample(c("low", "high"), n_obs, replace=TRUE)
  )
  
  # Create synthetic prediction matrices for two models
  # Create a list to store valid predictions
  individual_preds <- list()
  
  message("Creating synthetic predictions for testing")
  
  # Model 1: High values for 'low' group, low values for 'high' group
  pred_mat1 <- matrix(0, nrow=n_times, ncol=n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i/n_times
    pred_mat1[i, test_data$group == "low"] <- 0.8 * time_factor  
    pred_mat1[i, test_data$group == "high"] <- 0.2 * time_factor
  }
  individual_preds[["Model1"]] <- pred_mat1
  
  # Model 2: Similar pattern but more extreme
  pred_mat2 <- matrix(0, nrow=n_times, ncol=n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i/n_times
    pred_mat2[i, test_data$group == "low"] <- 0.9 * time_factor
    pred_mat2[i, test_data$group == "high"] <- 0.1 * time_factor
  }
  individual_preds[["Model2"]] <- pred_mat2
  
  # Test models have expected dimensions
  message(paste("Model1 dimensions:", paste(dim(pred_mat1), collapse="x")))
  message(paste("Model2 dimensions:", paste(dim(pred_mat2), collapse="x")))
  
  # Test simple averaging with our synthetic data
  ensemble_pred <- cifMatListAveraging(individual_preds, type="prob")
  
  # Verify predictions
  expect_true(is.matrix(ensemble_pred))
  expect_equal(nrow(ensemble_pred), n_times)
  expect_equal(ncol(ensemble_pred), n_obs)
  
  # Stratify predictions by group for last time point (final CIF value)
  low_indices <- which(test_data$group == "low")
  high_indices <- which(test_data$group == "high")
  
  low_probs <- ensemble_pred[n_times, low_indices]
  high_probs <- ensemble_pred[n_times, high_indices]
  
  # Compare group means - low should have higher risk
  mean_low <- mean(low_probs)
  mean_high <- mean(high_probs)
  
  message(sprintf("Ensemble prediction - Mean low: %.3f, Mean high: %.3f", mean_low, mean_high))
  
  # Low group should have higher risk (CIF values)
  expect_true(mean_low > mean_high, 
              label = sprintf("Ensemble prediction failed: mean low = %.3f, mean high = %.3f", 
                            mean_low, mean_high))
  
  # The mean for low group should be close to average of models (0.8 + 0.9)/2 = 0.85
  expect_true(abs(mean_low - 0.85) < 0.05, 
              label = sprintf("Low group mean (%.3f) not close to expected (0.85)", mean_low))
  
  # The mean for high group should be close to average of models (0.2 + 0.1)/2 = 0.15
  expect_true(abs(mean_high - 0.15) < 0.05,
              label = sprintf("High group mean (%.3f) not close to expected (0.15)", mean_high))
})