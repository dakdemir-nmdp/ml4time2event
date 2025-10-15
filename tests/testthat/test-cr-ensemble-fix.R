#' Test file for the ensemble predictions
#' 
#' This file contains synthetic tests to ensure the ensemble prediction
#' functionality works correctly with well-formed inputs.

library(testthat)
library(ml4time2event)

test_that("ensemble predictions work with synthetic data", {
  # Create synthetic test data
  set.seed(123)
  n_obs <- 20
  n_times <- 5
  eval_times <- seq(0.1, 1, by=0.2)
  
  # Create test data frame
  test_data <- data.frame(
    x1 = rnorm(n_obs),
    x2 = rnorm(n_obs),
    group = sample(c("low", "high"), n_obs, replace=TRUE)
  )
  
  # Create synthetic prediction matrices for two models
  # Model 1: High risk for 'low' group, low risk for 'high' group
  pred_mat1 <- matrix(0, nrow=n_times, ncol=n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i/n_times
    pred_mat1[i, test_data$group == "low"] <- 0.8 * time_factor  
    pred_mat1[i, test_data$group == "high"] <- 0.2 * time_factor 
  }
  
  # Model 2: Similar pattern but more extreme
  pred_mat2 <- matrix(0, nrow=n_times, ncol=n_obs)
  for (i in 1:n_times) {
    # Increasing CIF over time
    time_factor <- i/n_times
    pred_mat2[i, test_data$group == "low"] <- 0.9 * time_factor
    pred_mat2[i, test_data$group == "high"] <- 0.1 * time_factor
  }
  
  # Create prediction list
  model_predictions <- list(
    "Model1" = pred_mat1,
    "Model2" = pred_mat2
  )
  
  # Test simple averaging
  ensemble_pred <- cifMatListAveraging(model_predictions, type="prob")
  
  # Test ensemble structure
  expect_true(is.matrix(ensemble_pred))
  expect_equal(nrow(ensemble_pred), n_times)
  expect_equal(ncol(ensemble_pred), n_obs)
  
  # Verify ensemble predictions average the two models
  # For low group, should be around 0.85
  # For high group, should be around 0.15
  low_indices <- which(test_data$group == "low")
  high_indices <- which(test_data$group == "high")
  
  # Check final time point predictions
  expect_true(all(abs(ensemble_pred[n_times, low_indices] - 0.85) < 0.01))
  expect_true(all(abs(ensemble_pred[n_times, high_indices] - 0.15) < 0.01))
  
  # Test weighted averaging
  weights <- c(Model1 = 0.75, Model2 = 0.25)
  ensemble_weighted <- cifMatWeightedAveraging(model_predictions, weights, type="prob")
  
  # Test weighted ensemble structure
  expect_true(is.matrix(ensemble_weighted))
  expect_equal(dim(ensemble_weighted), dim(pred_mat1))
  
  # Verify weighted predictions match the weights
  # For low group: 0.75*0.8 + 0.25*0.9 = 0.6 + 0.225 = 0.825
  # For high group: 0.75*0.2 + 0.25*0.1 = 0.15 + 0.025 = 0.175
  expect_true(all(abs(ensemble_weighted[n_times, low_indices] - 0.825) < 0.01))
  expect_true(all(abs(ensemble_weighted[n_times, high_indices] - 0.175) < 0.01))
})