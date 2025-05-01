# Example 02: Feature Engineering for Survival Data

# --- 1. Load Libraries ---
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
if (!requireNamespace("recipes", quietly = TRUE)) {
  install.packages("recipes")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(survival)
library(recipes)
library(dplyr)

# Load the package being developed (assuming it's loaded or installed)
# library(ml4time2event) # Uncomment if needed and package is installed/loaded

# --- 2. Load Data ---
# Using the Primary Biliary Cirrhosis (PBC) dataset from the 'survival' package
data(pbc, package = "survival")

# Use baseline data only
pbc_baseline <- pbc[!duplicated(pbc$id), ]

# Define outcome (ensure time > 0, status is 0/1)
pbc_baseline <- pbc_baseline %>%
  filter(time > 0) %>%
  mutate(status_event = ifelse(status == 2, 1, 0)) %>%
  select(-status) # Remove original status column

# --- 3. Initial Data Split (Optional but good practice) ---
# For demonstration, we'll use the whole dataset for preprocessing steps.
# In a real scenario, split into training and testing sets *before* preprocessing.
# set.seed(123)
# split <- initial_split(pbc_baseline, prop = 0.75, strata = status_event)
# train_data <- training(split)
# test_data  <- testing(split)
# For now, use pbc_baseline as the data to preprocess
data_to_prep <- pbc_baseline

# --- 4. Define Predictors and Outcome ---
outcome_vars <- c("time", "status_event")
# Identify potential predictors (excluding id and outcome)
predictor_vars <- setdiff(names(data_to_prep), c("id", outcome_vars))

# --- 5. Create a Recipe for Preprocessing ---
# A recipe defines the sequence of steps for feature engineering.
# Specify predictors in vars, outcomes separately, and roles for predictors
pbc_recipe <- recipe(x = data_to_prep, vars = predictor_vars, roles = rep("predictor", length(predictor_vars)), outcomes = outcome_vars) %>%
  # Step 1: Handle Data Types (convert factors, etc.)
  # 'sex', 'ascites', 'hepato', 'spiders', 'edema', 'stage' are potential factors
  step_mutate(across(c(sex, ascites, hepato, spiders, edema, stage), as.factor)) %>%
  # Step 2: Impute Missing Values
  # Using median for numeric and mode for nominal (categorical)
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  # Step 3: Create Derived Features (Example: log transform skewed variables)
  step_log(bili, base = 10, offset = 0.001) %>% # Offset for values near 0
  step_log(protime, base = 10, offset = 0.001) %>%
  # Step 4: Dummy Code Categorical Variables
  # Creates binary indicator variables for each factor level
  step_dummy(all_nominal_predictors()) %>%
  # Step 5: Remove Zero-Variance Predictors (created if a factor level had no variance)
  step_zv(all_predictors()) %>%
  # Step 6: Normalize Numeric Predictors (center and scale)
  step_normalize(all_numeric_predictors())

# --- 6. Prepare (Train) the Recipe ---
# The prep() function estimates parameters needed for the steps (e.g., medians, means, std devs)
# Normally, you prep() on the training data only.
prepped_recipe <- prep(pbc_recipe, training = data_to_prep)

# --- 7. Apply (Bake) the Recipe ---
# The bake() function applies the prepared recipe to new data.
processed_data <- bake(prepped_recipe, new_data = data_to_prep)

# --- 8. Inspect Processed Data ---
cat("\n--- Processed Data Summary ---\n")
summary(processed_data)
cat("\n--- Processed Data Structure ---\n")
glimpse(processed_data) # dplyr's version of str()

# Note: The processed data now contains the outcome variables ('time', 'status_event')
# and the engineered predictor features, ready for model training.

# --- End of Example 02 ---
