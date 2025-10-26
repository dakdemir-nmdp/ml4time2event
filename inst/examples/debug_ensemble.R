# Load required libraries
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
print(summary(cox_model))

# Build survival regression model
cat("\nBuilding SurvReg model...\n")
# Get the parameter names for SurvModel_SurvReg
print(args(SurvModel_SurvReg))
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
print(args(Predict_SurvModel_Cox))
cox_pred <- Predict_SurvModel_Cox(
  modelout = cox_model,
  newdata = test_data,
  new_times = times
)
print(str(cox_pred))

cat("\nMaking predictions with SurvReg model...\n")
print(args(Predict_SurvModel_SurvReg))
survreg_pred <- Predict_SurvModel_SurvReg(
  modelout = survreg_model,
  newdata = test_data,
  new_times = times
)
print(str(survreg_pred))

# Convert to risk (looks like Probs is already risk probabilities, not survival)
cox_risk <- cox_pred$Probs
survreg_risk <- survreg_pred$Probs

# Make sure dimensions are correct
cat("\nDimensions of predictions:\n")
print(paste("Cox Risk dimensions:", paste(dim(cox_risk), collapse=" x ")))
print(paste("SurvReg Risk dimensions:", paste(dim(survreg_risk), collapse=" x ")))

# Trim survreg_risk to match cox_risk in time points
survreg_risk <- survreg_risk[2:11, ]  # Skip time 0
print(paste("Adjusted SurvReg Risk dimensions:", paste(dim(survreg_risk), collapse=" x ")))

# Create prediction list
individual_preds <- list(
  cox = cox_risk,
  survreg = survreg_risk
)

# Average the predictions
cat("\nAveraging predictions...\n")
ensemble_pred <- EnsemblePredictions(individual_preds, type="prob")
print(dim(ensemble_pred))

# Create plots directory
plots_dir <- file.path(getwd(), "ensemble_plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Plot a sample subject
subject_idx <- 42
pdf_file <- file.path(plots_dir, paste0("subject_", subject_idx, "_predictions.pdf"))
pdf(pdf_file, width = 8, height = 6)

plot(NULL, xlim = c(0, max(times)), ylim = c(0, 1), 
     xlab = "Time", ylab = "Risk", 
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
cat("\nPlot saved to:", pdf_file, "\n")