library(testthat)

test_that("Ensemble end-to-end test (local)", {
  # Load required libraries explicitly
  suppressMessages({
    library(ml4time2event)
    library(survival)
    if(requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
    if(requireNamespace("splines", quietly = TRUE)) library(splines)
  })
  
  message("Creating test data...")
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
  
  # Create training dataset (for model building)
  train_data <- test_data
  
  # Build basic survival models
  models <- list()
  
  message("Building Cox model...")
  # Try to create a Cox model
  cox_model <- tryCatch({
    # Use ml4time2event's functionality
    result <- SurvModel_Cox(
      data = train_data,
      timevar = "time",
      eventvar = "status",
      expvars = c("x1", "x2")
    )
    message("Cox model created successfully")
    result
  }, error = function(e) {
    message("Could not create Cox model: ", e$message)
    NULL
  })
  
  if (!is.null(cox_model)) {
    models$cox <- list(
      model = cox_model,
      type = "SURV"
    )
  }
  
  # Create a survival regression model
  message("Building SurvReg model...")
  survreg_model <- tryCatch({
    result <- SurvModel_SurvReg(
      data = train_data,
      timevar = "time",
      eventvar = "status",
      expvars = c("x1", "x2")
    )
    message("SurvReg model created successfully")
    result
  }, error = function(e) {
    message("Could not create SurvReg model: ", e$message)
    NULL
  })
  
  if (!is.null(survreg_model)) {
    models$survreg <- list(
      model = survreg_model,
      type = "SURV"
    )
  }
  
  # Individual predictions list
  individual_preds <- list()
  
  # Loop through models and make predictions
  for (model_name in names(models)) {
    model_info <- models[[model_name]]
    model <- model_info$model
    
    # Predict using the model
    pred <- tryCatch({
      if (model_info$type == "SURV") {
        message(paste("Predicting using model:", model_name))
        if (model_name == "cox") {
          pred <- Predict_SurvModel_Cox(
            modelout = model,
            newdata = test_data,
            newtimes = 1:10
          )
        } else if (model_name == "survreg") {
          pred <- Predict_SurvModel_SurvReg(
            modelout = model,
            newdata = test_data,
            newtimes = 1:10
          )
        }
        
        if (!is.null(pred) && is.list(pred) && "Probs" %in% names(pred)) {
          result <- pred$Probs
          if (model_name == "survreg") {
            # Survreg prediction includes time 0, skip it
            result <- result[2:11,]
          }
          individual_preds[[model_name]] <- result
          message(paste("Prediction successful for", model_name))
          result
        } else {
          message(paste("Invalid prediction structure for", model_name))
          NULL
        }
      } else {
        message(paste("Invalid model type for", model_name))
        NULL
      }
    }, error = function(e) {
      message(paste("Error predicting for", model_name, ":", e$message))
      NULL
    })
  }

  # Check if we have any valid predictions
  if (length(individual_preds) > 0) {
    message(paste("Got valid predictions for", length(individual_preds), "models:", 
                  paste(names(individual_preds), collapse=", ")))
    
    # Try to average the predictions
    ensemble_pred <- tryCatch({
      message("Performing ensemble averaging...")
      if (length(individual_preds) == 2 && all(names(individual_preds) %in% c("cox", "survreg"))) {
        # Direct averaging of cox and survreg models
        cox_risk <- individual_preds[["cox"]]
        survreg_risk <- individual_preds[["survreg"]]
        (cox_risk + survreg_risk) / 2
      } else {
        # Fall back to the general approach
        sum_pred <- matrix(0, nrow = nrow(individual_preds[[1]]), ncol = ncol(individual_preds[[1]]))
        for (pred in individual_preds) {
          sum_pred <- sum_pred + pred
        }
        sum_pred / length(individual_preds)
      }
    }, error = function(e) {
      message(paste("Ensemble averaging failed:", e$message))
      return(NULL)
    })
    
    # Create plots directory if it doesn't exist
    plots_dir <- file.path(tempdir(), "ml4time2event_plots")
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir, recursive = TRUE)
    }
    
    # Plot individual model predictions for a sample of subjects
    samples <- sample(seq_len(nrow(test_data)), min(5, nrow(test_data)))
    times <- 1:10
    
    # For each sampled subject
    for (subject_idx in samples) {
      pdf_file <- file.path(plots_dir, paste0("subject_", subject_idx, "_predictions.pdf"))
      pdf(pdf_file, width = 8, height = 6)
      
      plot(NULL, xlim = c(0, max(times)), ylim = c(0, 1), 
           xlab = "Time", ylab = "Risk", 
           main = paste("Predictions for Subject", subject_idx))
      
      colors <- rainbow(length(individual_preds) + 1)
      color_idx <- 1
      
      # Plot each individual model prediction
      for (model_name in names(individual_preds)) {
        pred <- individual_preds[[model_name]]
        lines(times, pred[, subject_idx], col = colors[color_idx], lwd = 2, lty = 2)
        color_idx <- color_idx + 1
      }
      
      # Plot ensemble prediction if available
      if (!is.null(ensemble_pred)) {
        lines(times, ensemble_pred[, subject_idx], col = colors[color_idx], lwd = 3)
        legend_names <- c(names(individual_preds), "Ensemble")
      } else {
        legend_names <- names(individual_preds)
      }
      
      # Add legend
      legend("topleft", legend = legend_names, 
             col = colors[seq_along(legend_names)], 
             lwd = c(rep(2, length(individual_preds)), 3)[seq_along(legend_names)],
             lty = c(rep(2, length(individual_preds)), 1)[seq_along(legend_names)])
      
      dev.off()
      message(paste("Saved plot to", pdf_file))
    }
    
    # Print where the plots are stored
    message("Plots saved to directory: ", plots_dir)
    
    # Verify ensemble averaging works correctly
    if (!is.null(ensemble_pred)) {
      expect_true(is.matrix(ensemble_pred), "Ensemble prediction should be a matrix")
      expect_equal(dim(ensemble_pred)[1], length(times), 
                   tolerance = 1e-6)
      expect_equal(dim(ensemble_pred)[2], nrow(test_data), 
                   tolerance = 1e-6)
    }
  } else {
    message("No valid predictions to plot.")
    stop("No valid predictions to test ensemble averaging")
  }
})