# Example 02: Feature Selection using Random Forest Importance

# --- 1. Load Libraries ---
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
if (!requireNamespace("randomForestSRC", quietly = TRUE)) { # RF model dependency
  install.packages("randomForestSRC")
}
library(survival)
library(dplyr)
library(randomForestSRC) # Needed for SurvModel_RF likely

# --- Load the ml4time2event package ---
# IMPORTANT: Assumes package functions are available.
# Use devtools::load_all() or library(ml4time2event).

# --- 2. Load and Prepare Data ---
# Using the Primary Biliary Cirrhosis (PBC) dataset
data(pbc, package = "survival")
pbc_baseline <- pbc[!duplicated(pbc$id), ] %>%
  filter(time > 0) %>%
  mutate(
    status_event = ifelse(status == 2, 1, 0),
    sex = factor(sex), stage = factor(stage)
    # Add other conversions as needed
  ) %>%
  select(-id, -status)

time_var <- "time"
event_var <- "status_event"
predictor_vars <- setdiff(names(pbc_baseline), c(time_var, event_var))

# Handle missing values (simple median/mode imputation)
for (col in predictor_vars) {
  if (is.numeric(pbc_baseline[[col]])) {
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- median(pbc_baseline[[col]], na.rm = TRUE)
  } else if (is.factor(pbc_baseline[[col]])) {
    mode_val <- names(which.max(table(pbc_baseline[[col]])))
    pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
  }
  if (is.character(pbc_baseline[[col]])) {
     mode_val <- names(which.max(table(pbc_baseline[[col]])))
     pbc_baseline[[col]][is.na(pbc_baseline[[col]])] <- mode_val
     pbc_baseline[[col]] <- factor(pbc_baseline[[col]])
  }
}

# --- 3. Data Splitting ---
set.seed(123)
train_indices <- sample(1:nrow(pbc_baseline), size = floor(0.75 * nrow(pbc_baseline)))
train_data <- pbc_baseline[train_indices, ]
test_data  <- pbc_baseline[-train_indices, ]

# --- 4. Train Random Forest to Get Importance ---
# Assuming SurvModel_RF is loaded and uses randomForestSRC::rfsrc
cat("Training Random Forest model (all variables) for importance...\n")
# Need to know the actual function name and parameters (e.g., ntree)
# Assuming SurvModel_RF returns an object containing importance scores
# Check R/models/surv_random_forest.R for the actual implementation details.
# Placeholder call:
rf_model_full <- SurvModel_RF(
  data = train_data,
  expvars = predictor_vars,
  timevar = time_var,
  eventvar = event_var,
  ntree = 100 # Example parameter, adjust as needed
  # Add other necessary parameters like samplesize if required by the function
)

# --- 5. Extract and Select Top Features ---
n_top_vars <- 10 # Number of top variables to select

if (!is.null(rf_model_full) && !is.null(rf_model_full$hd.obj$importance)) {
  importance_scores <- rf_model_full$hd.obj$importance
  # Handle matrix or vector importance output
  if (is.matrix(importance_scores)) {
    # Assuming the first column is the relevant importance score
    sorted_importance <- sort(importance_scores[, 1], decreasing = TRUE)
  } else {
    sorted_importance <- sort(importance_scores, decreasing = TRUE)
  }
  top_predictor_vars <- names(sorted_importance)[1:min(length(sorted_importance), n_top_vars)]
  cat("\nTop", n_top_vars, "variables selected based on RF importance:\n")
  print(top_predictor_vars)

  # Visualize Importance (optional)
  dotchart(sorted_importance[1:min(length(sorted_importance), 20)],
           main = "Variable Importance (Random Forest)", pch=19)

} else {
  warning("Could not extract importance scores from Random Forest model. Using all variables.")
  top_predictor_vars <- predictor_vars
}

# --- 6. Train Model with Selected Features ---
# Example: Train a Cox model using only the top features
cat("\nTraining Cox model with top", length(top_predictor_vars), "variables...\n")
cox_model_top_vars <- SurvModel_Cox(
  data = train_data,
  expvars = top_predictor_vars,
  timevar = time_var,
  eventvar = event_var
)

# --- 7. Predict and Evaluate Model with Selected Features ---
cat("Predicting and evaluating Cox model (top variables)...\n")
cox_predictions_top_vars <- Predict_SurvModel_Cox(
  modelout = cox_model_top_vars,
  newdata = test_data
)

test_times <- test_data[[time_var]]
test_events <- test_data[[event_var]]

c_index_top_vars <- timedepConcordance(
  predsurv = cox_predictions_top_vars$Probs,
  predsurvtimes = cox_predictions_top_vars$Times,
  obstimes = test_times,
  obsevents = test_events
)

cat("\n--- Evaluation Summary (Cox Model with Top Variables) ---\n")
print(summary(c_index_top_vars$AppErr$matrix))

# Plot C-index for comparison (if desired, compare with C-index from all vars)
# plot(c_index_results, col = "blue") # From Example 01 (all vars)
# lines(c_index_top_vars$times, c_index_top_vars$AppErr$matrix, col = "red") # Top vars
# legend("bottomright", legend = c("All Vars", "Top Vars"), col = c("blue", "red"), lty = 1)


# --- End of Example 02 ---
