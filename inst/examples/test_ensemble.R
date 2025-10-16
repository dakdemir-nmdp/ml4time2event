# Test script for ensemble methods
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

# Average the predictions
cat("\nAveraging predictions...\n")
ensemble_pred <- (cox_risk + survreg_risk) / 2
cat("Ensemble prediction dimensions:", paste(dim(ensemble_pred), collapse=" x "), "\n")

# Create plots directory
plots_dir <- file.path(getwd(), "ensemble_plot.pdf")
pdf(plots_dir, width = 10, height = 6)

# Plot a specific subject
subject_idx <- 42

plot(NULL, xlim = c(0, max(times)), ylim = c(0, 1), 
     xlab = "Time", ylab = "Risk Probability", 
     main = paste("Predictions for Subject", subject_idx),
     cex.main = 1.5, cex.lab = 1.2)

# Plot individual model predictions
colors <- c("blue", "green", "red")

# Plot individual model predictions
lines(times, cox_risk[, subject_idx], col = colors[1], lwd = 2, lty = 2)
lines(times, survreg_risk[, subject_idx], col = colors[2], lwd = 2, lty = 2)

# Plot ensemble prediction
lines(times, ensemble_pred[, subject_idx], col = colors[3], lwd = 3)
legend_names <- c("Cox", "SurvReg", "Ensemble")

# Add legend
legend("topleft", legend = legend_names, 
       col = colors[1:length(legend_names)], 
       lwd = c(rep(2, 2), 3),
       lty = c(rep(2, 2), 1),
       cex = 1.2)

dev.off()

cat("\nPlot saved to:", plots_dir, "\n")