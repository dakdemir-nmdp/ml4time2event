test_that("Ensemble end-to-end test", {
  # Skip on CRAN
  skip_on_cran()
  
  # Create test data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  status <- sample(0:1, n, replace=TRUE, prob=c(0.3, 0.7))
  time <- rexp(n) * 10
  
  # Create a simple dataset
  test_data <- data.frame(
    x1 = x1,
    x2 = x2,
    status = status,
    time = time
  )
  
  # Build Cox model
  cox_model <- SurvModel_Cox(
    data = test_data,
    timevar = "time",
    eventvar = "status",
    expvars = c("x1", "x2")
  )
  
  # Build survival regression model
  survreg_model <- SurvModel_SurvReg(
    data = test_data,
    timevar = "time",
    eventvar = "status",
    expvars = c("x1", "x2")
  )
  
  # Make predictions
  times <- 1:10
  
  # Cox predictions
  cox_pred <- Predict_SurvModel_Cox(
    modelout = cox_model,
    newdata = test_data,
    new_times = times
  )
  message("Prediction successful for cox")
  
  # SurvReg predictions
  survreg_pred <- Predict_SurvModel_SurvReg(
    modelout = survreg_model,
    newdata = test_data,
    new_times = times
  )
  message("Prediction successful for survreg")
  
  # Convert to risk
  cox_risk <- cox_pred$Probs
  survreg_risk <- survreg_pred$Probs[2:11,]  # Skip time 0
  
  # Create prediction list
  individual_preds <- list(
    cox = cox_risk,
    survreg = survreg_risk
  )
  
  message(paste("Got valid predictions for", length(individual_preds), "models:", 
            paste(names(individual_preds), collapse=", ")))
  
  # Average the predictions
  message("Performing ensemble averaging...")
  ensemble_pred <- (cox_risk + survreg_risk) / 2
  
  # Create plots directory
  plots_dir <- file.path(tempdir(), "ml4time2event_plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Sample a few subjects to plot
  samples <- sample(seq_len(ncol(ensemble_pred)), min(5, ncol(ensemble_pred)))
  
  for (subject_idx in samples) {
    pdf_file <- file.path(plots_dir, paste0("subject_", subject_idx, "_predictions.pdf"))
    pdf(pdf_file, width = 8, height = 6)
    
    plot(NULL, xlim = c(0, max(times)), ylim = c(0, 1), 
         xlab = "Time", ylab = "Risk Probability", 
         main = paste("Predictions for Subject", subject_idx),
         cex.main = 1.5, cex.lab = 1.2)
    
    # Plot individual model predictions
    colors <- c("blue", "green", "red")
    lines(times, cox_risk[, subject_idx], col = colors[1], lwd = 2, lty = 2)
    lines(times, survreg_risk[, subject_idx], col = colors[2], lwd = 2, lty = 2)
    
    # Plot ensemble prediction
    lines(times, ensemble_pred[, subject_idx], col = colors[3], lwd = 3)
    
    # Add legend
    legend("topleft", legend = c("Cox", "SurvReg", "Ensemble"), 
           col = colors, 
           lwd = c(2, 2, 3),
           lty = c(2, 2, 1),
           cex = 1.2)
    
    dev.off()
    message(paste("Saved plot to", pdf_file))
  }
  
  message("Plots saved to directory: ", plots_dir)
  
  # Verify ensemble averaging works correctly
  expect_true(is.matrix(ensemble_pred), info = "Ensemble prediction should be a matrix")
  expect_equal(dim(ensemble_pred), c(length(times), ncol(cox_risk)), 
               info = "Ensemble prediction should have the correct dimensions")
})
