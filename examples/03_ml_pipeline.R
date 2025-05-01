# Example 03: Machine Learning Pipeline for Survival Analysis

# --- 1. Load Libraries ---
# Install necessary tidymodels packages if not present
if (!requireNamespace("tidymodels", quietly = TRUE)) {
  install.packages("tidymodels")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
# Specific model packages might be needed depending on the chosen engine
# e.g., install.packages("randomForestSRC") for random forest survival
if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
  install.packages("randomForestSRC")
}
if (!requireNamespace("ranger", quietly = TRUE)) { # Add ranger
  install.packages("ranger")
}

library(tidymodels)
library(survival)
library(dplyr)
# library(randomForestSRC) # Comment out or remove
library(ranger) # Add ranger

# Load the package being developed (assuming it's loaded or installed)
# library(ml4time2event) # Uncomment if needed and package is installed/loaded

# --- 2. Load and Prepare Data ---
data(pbc, package = "survival")
pbc_baseline <- pbc[!duplicated(pbc$id), ] %>%
  filter(time > 0) %>%
  mutate(
    status_event = ifelse(status == 2, 1, 0),
    # Convert potential factors explicitly for the recipe
    sex = factor(sex),
    ascites = factor(ascites),
    hepato = factor(hepato),
    spiders = factor(spiders),
    edema = factor(edema),
    stage = factor(stage)
  ) %>%
  select(-status, -id) # Remove original status and id

# --- 3. Data Splitting ---
# Stratify by the event status BEFORE creating the Surv object and removing original columns
set.seed(123)
split <- initial_split(pbc_baseline, prop = 0.75, strata = status_event)
train_data_orig <- training(split)
test_data_orig  <- testing(split)

# Create the Surv object and remove original time/event columns AFTER splitting
train_data <- train_data_orig %>%
  mutate(surv_obj = Surv(time, status_event)) %>%
  select(-time, -status_event)
test_data <- test_data_orig %>%
  mutate(surv_obj = Surv(time, status_event)) %>%
  select(-time, -status_event)


# --- 4. Define Preprocessing Recipe ---
# Use the Surv object as the outcome in the recipe
pbc_recipe <- recipe(surv_obj ~ ., data = train_data) %>%
  # Impute missing values
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  # Log transform skewed numeric predictors
  step_log(bili, base = 10, offset = 0.001) %>%
  step_log(protime, base = 10, offset = 0.001) %>%
  # Create dummy variables for categorical predictors
  step_dummy(all_nominal_predictors()) %>%
  # Remove zero-variance predictors
  step_zv(all_predictors()) %>%
  # Normalize numeric predictors
  step_normalize(all_numeric_predictors())

# --- 5. Define the Model (Using Cox PH only for this simplified example) ---
# # Using Random Survival Forest as an example
# # Other models like Cox (survival::coxph) or Boosted Trees (mboost::glmboost) can also be used
# # Note: Ensure the chosen engine supports Surv objects
# rf_spec <- rand_forest(trees = 100, mtry = tune(), min_n = tune()) %>%
#   set_engine("ranger", importance = "permutation") %>% # Use ranger package
#   set_mode("censored regression") # Specify the mode for survival analysis
#
# # --- 6. Create Workflow ---
# # Combine the recipe and model specification into a workflow
# rf_workflow <- workflow() %>%
#   add_recipe(pbc_recipe) %>%
#   add_model(rf_spec)
#
# # --- 7. Set up Cross-Validation ---
# set.seed(456)
# cv_folds <- vfold_cv(train_data, v = 5, strata = surv_obj) # 5-fold CV
#
# # --- 8. Hyperparameter Tuning (Optional but Recommended) ---
# # Define a grid for tuning mtry and min_n
# rf_grid <- grid_regular(
#   mtry(range = c(2, 10)), # Adjust range based on number of predictors after recipe
#   min_n(range = c(5, 15)),
#   levels = 3
# )
#
# # Define the metric for tuning (e.g., Concordance Index)
# # Note: tidymodels uses yardstick; concordance might need custom setup or use available metrics
# # For simplicity here, we'll skip explicit tuning metric setup and focus on fitting.
# # In practice, use tune_grid with appropriate metrics.
#
# # --- 9. Fit Model(s) using Cross-Validation ---
# # For simplicity, fit without tuning first. Replace with tune_grid for hyperparameter optimization.
# # We need to finalize the model spec if not tuning
# final_rf_spec <- rand_forest(trees = 100, mtry = 5, min_n = 10) %>% # Example fixed parameters
#   set_engine("ranger", importance = "permutation") %>% # Use ranger package
#   set_mode("censored regression")
#
# final_rf_workflow <- workflow() %>%
#   add_recipe(pbc_recipe) %>%
#   add_model(final_rf_spec)
#
# # Fit the model to the cross-validation folds
# # Note: Direct fitting with fit_resamples might not directly give survival-specific metrics easily.
# # Often, manual looping or specific survival packages integrated with tidymodels are used.
# # Let's fit the final workflow to the full training set for now.
# final_fit <- fit(final_rf_workflow, data = train_data)
#
# # --- 10. Evaluate on Test Set ---
# # Predict on the test set
# test_predictions <- predict(final_fit, new_data = test_data, type = "time") # Predict survival time
#
# # Combine predictions with actual test data outcomes
# test_results <- bind_cols(
#   test_data %>% select(surv_obj), # Actual survival object
#   test_predictions %>% rename(predicted_time = .pred_time) # Predicted time
# )
#
# # Calculate performance metrics (e.g., Concordance Index on test set)
# # This requires extracting time and status from the Surv object
# test_time <- test_results$surv_obj[, 1]
# test_status <- test_results$surv_obj[, 2]
#
# # Concordance calculation often needs predicted risk scores, not just time.
# # randomForestSRC can provide mortality scores (risk). Let's refit slightly differently
# # to get risk predictions if possible, or use a simpler model like Cox for C-index demo.

# --- Using Cox for simpler C-index demo ---
cox_spec <- proportional_hazards() %>%
  set_engine("survival") %>%
  set_mode("censored regression")

cox_workflow <- workflow() %>%
  add_recipe(pbc_recipe) %>%
  add_model(cox_spec)

cox_fit <- fit(cox_workflow, data = train_data)

# Predict risk (linear predictors) on the test set
cox_test_pred <- predict(cox_fit, new_data = test_data, type = "linear_pred")

# Calculate Concordance
# Need to extract actual time and status from the test set Surv object
# Use the test_data_orig which still has the original time/status columns
test_surv_obj_actual <- Surv(test_data_orig$time, test_data_orig$status_event)
c_index_cox <- concordance(
  test_surv_obj_actual ~ cox_test_pred$.pred_linear_pred,
  reverse = TRUE # Higher score means higher risk (lower survival time)
)$concordance

cat("\n--- Test Set Evaluation (Cox Model) ---\n")
cat("Concordance Index (C-index):", c_index_cox, "\n")


# --- 11. Final Model ---
# The 'final_fit' (or 'cox_fit') object is the trained pipeline ready for prediction on new data.
# You can save this object using saveRDS().

# --- End of Example 03 ---
