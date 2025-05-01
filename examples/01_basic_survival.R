# Example 01: Basic Survival Analysis using ml4time2event

# --- 1. Load Libraries ---
# Ensure survival package is available for data and core functions
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
library(survival)
library(dplyr)
library(ml4time2event)

# --- Load the ml4time2event package ---
# IMPORTANT: This script assumes the 'ml4time2event' package functions are available.
# During development, use devtools::load_all() in the R console.
# If installed, uncomment the line below:
# library(ml4time2event)

# --- 2. Load and Prepare Data ---
# Using the Primary Biliary Cirrhosis (PBC) dataset
data(pbc, package = "survival")

# Use baseline data only, basic cleaning
pbc_baseline <- pbc[!duplicated(pbc$id), ] %>%
  filter(time > 0) %>%
  mutate(
    status_event = ifelse(status == 2, 1, 0), # 1 for death, 0 for censored
    # Convert relevant columns to factors or numeric
    sex = factor(sex),
    stage = factor(stage)
    # Add other conversions as needed by VariableProfile or models
  ) %>%
  select(-id, -status) # Remove original id and status

# Define outcome and predictor variables
time_var <- "time"
event_var <- "status_event"
predictor_vars <- setdiff(names(pbc_baseline), c(time_var, event_var))

# Handle missing values (simple median/mode imputation for demonstration)
# A more robust approach might use recipes or dedicated imputation packages
for (col in predictor_vars) {
  if (is.numeric(pbc_baseline[[col]])) {
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- median(pbc_baseline[[col]], na.rm = TRUE)
  } else if (is.factor(pbc_baseline[[col]])) {
    # Simple mode imputation for factors
    mode_val <- names(which.max(table(pbc_baseline[[col]])))
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
  }
  # Handle character columns if any remain
  if (is.character(pbc_baseline[[col]])) {
     mode_val <- names(which.max(table(pbc_baseline[[col]])))
     pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
     pbc_baseline[[col]] <- factor(pbc_baseline[[col]]) # Convert to factor after imputation
  }
}

# --- 3. Data Splitting (Optional but Recommended) ---
set.seed(123)
train_indices <- sample(1:nrow(pbc_baseline), size = floor(0.75 * nrow(pbc_baseline)))
train_data <- pbc_baseline[train_indices, ]
test_data  <- pbc_baseline[-train_indices, ]

# --- 4. Train Cox Model using ml4time2event ---
# Assuming SurvModel_Cox and VariableProfile are loaded from the package
cat("Training Cox model...\n")
# Note: VariableProfile function needs to be defined or loaded from the package.
# If it's not available, you might need to remove the varprof part or mock it.
# For now, assuming it exists.
cox_model_output <- SurvModel_Cox(
  data = train_data,
  expvars = predictor_vars,
  timevar = time_var,
  eventvar = event_var
)

# --- 5. Predict using the Trained Model ---
cat("Predicting on test data...\n")
# Assuming Predict_SurvModel_Cox is loaded
cox_predictions <- Predict_SurvModel_Cox(
  modelout = cox_model_output,
  newdata = test_data
)

# --- 6. Evaluate the Model ---
cat("Evaluating model performance (Time-dependent C-index)...\n")
# Assuming timedepConcordance is loaded
# Extract observed times and events from the test set
test_times <- test_data[[time_var]]
test_events <- test_data[[event_var]]

# Calculate time-dependent concordance
# Note: Ensure prediction times align reasonably with evaluation needs.
# Using prediction times directly for evaluation here.
c_index_results <- timedepConcordance(
  predsurv = cox_predictions$Probs,
  predsurvtimes = cox_predictions$Times,
  obstimes = test_times,
  obsevents = test_events
  # ctimes = Optional specific times for C-index calculation
)

# Display summary of C-index over time
print(summary(c_index_results$AppErr$matrix)) # Example: Print summary of apparent error C-index

# Plot C-index over time
plot(c_index_results)
title("Time-dependent Concordance Index (Cox Model)")

# --- 7. Basic Visualization of Predictions (Example for one patient) ---
patient_index <- 1
plot(cox_predictions$Times, cox_predictions$Probs[, patient_index], type = 's',
     xlab = "Time (days)", ylab = "Predicted Survival Probability",
     main = paste("Predicted Survival Curve for Patient", patient_index),
     ylim = c(0, 1))

# --- End of Example 01 ---
