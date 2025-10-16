library(testthat)
library(ml4time2event)
library(survival)

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
cat("Building Cox model...\n")
cox_model <- SurvModel_Cox(
  data = test_data,
  timevar = "time",
  eventvar = "status",
  expvars = c("x1", "x2")
)
cat("Cox model created successfully\n")

# Build survival regression model
cat("\nBuilding SurvReg model...\n")
survreg_model <- SurvModel_SurvReg(
  data = test_data,
  timevar = "time",
  eventvar = "status",
  expvars = c("x1", "x2")
)
cat("SurvReg model created successfully\n")

# Make predictions
times <- 1:10
cat("\nMaking predictions with Cox model...\n")
cox_pred <- Predict_SurvModel_Cox(
  modelout = cox_model,
  newdata = test_data,
  newtimes = times
)
cat("Prediction successful for cox\n")

cat("\nMaking predictions with SurvReg model...\n")
survreg_pred <- Predict_SurvModel_SurvReg(
  modelout = survreg_model,
  newdata = test_data,
  newtimes = times
)
cat("Prediction successful for survreg\n")

# Convert to risk
cox_risk <- cox_pred$Probs
survreg_risk <- survreg_pred$Probs[2:11,]  # Skip time 0

# Create prediction list
individual_preds <- list(
  cox = cox_risk,
  survreg = survreg_risk
)

cat(paste("\nGot valid predictions for", length(individual_preds), "models:", 
            paste(names(individual_preds), collapse=", "), "\n"))

# Average the predictions
cat("\nPerforming ensemble averaging...\n")
ensemble_pred <- (cox_risk + survreg_risk) / 2

# Create plots directory
plots_dir <- file.path(tempdir(), "ml4time2event_plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Plot a few subjects
samples <- c(100, 70, 3, 20, 4) # Specific samples for reproducibility
for (subject_idx in samples) {
  pdf_file <- file.path(plots_dir, paste0("subject_", subject_idx, "_predictions.pdf"))
  pdf(pdf_file, width = 8, height = 6)
  
  plot(NULL, xlim = c(0, max(times)), ylim = c(0, 1), 
       xlab = "Time", ylab = "Risk Probability", 
       main = paste("Predictions for Subject", subject_idx))
  
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
         lty = c(2, 2, 1))
  
  dev.off()
  cat(paste("Saved plot to", pdf_file, "\n"))
}

cat("Plots saved to directory: ", plots_dir, "\n")

# Simple test
expect_true(is.matrix(ensemble_pred))