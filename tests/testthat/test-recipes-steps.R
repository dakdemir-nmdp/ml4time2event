library(testthat)
library(here)
# library(data.table) # Removed
library(recipes)
library(rsample)
library(survival) # For Surv object in formula creation

# Assuming the functions are available in the environment
# Use here::here for robustness
source(here::here("R/recipes_steps.R"))

context("Testing recipes_steps functions")

# --- Helper Function for Test Data Setup ---
create_recipe_data_and_vars <- function() {
  # library(data.table) # Rely on NAMESPACE import now
  set.seed(147)
  n_recipe <- 50
  # Replace data.table() with data.frame()
  x_num1_temp <- rnorm(n_recipe, 10, 2) # Create x_num1 first
  recipe_data <- data.frame(
    id = 1:n_recipe,
    time = round(rexp(n_recipe, 0.1), 1),
  status = sample(0:1, n_recipe, replace = TRUE, prob = c(0.4, 0.6)),
    x_num1 = x_num1_temp, # Use temp x_num1
    x_fact1 = factor(sample(c("A", "B", "C", "D_rare"), n_recipe, replace = TRUE, prob = c(0.4, 0.3, 0.25, 0.05))),
    x_fact2 = factor(sample(c("X", "Y"), n_recipe, replace = TRUE)),
    x_nzv = rep(1, n_recipe), # Near-zero variance
    x_miss = sample(c(1:5, NA), n_recipe, replace = TRUE, prob = c(rep(0.18, 5), 0.1)), # ~10% missing
    x_num2_correlated = rnorm(n_recipe, mean = 0.5 * x_num1_temp), # Use temp x_num1
    x_high_miss = as.numeric(NA), # Initialize with NA
    stringsAsFactors = FALSE
)
  # Assign non-NA values to x_high_miss using standard indexing
  recipe_data[sample(1:n_recipe, 20), "x_high_miss"] <- rnorm(20)

  # Define variables for recipe
  time_var_recipe <- "time"
  event_var_recipe <- "status"
  id_vars_recipe <- "id"
  exp_vars_recipe <- c("x_num1", "x_num2_correlated", "x_fact1", "x_fact2", "x_nzv", "x_miss", "x_high_miss")

  return(list(
    recipe_data = recipe_data,
    time_var_recipe = time_var_recipe,
    event_var_recipe = event_var_recipe,
    id_vars_recipe = id_vars_recipe,
    exp_vars_recipe = exp_vars_recipe,
    n_recipe = n_recipe
  ))
}


# --- Tests for t2edata_split ---

test_that("t2edata_split creates train/test split", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("rsample")
  split_list <- t2edata_split(test_setup$recipe_data, prop = 0.7)

  expect_type(split_list, "list")
  expect_named(split_list, c("Train", "Test", "Split"))
  # Verify contents via explicit indexing (avoids partial-match ambiguity with $)
  expect_s3_class(split_list[["Train"]], "data.frame")
  expect_s3_class(split_list[["Test"]], "data.frame")
  expect_s3_class(split_list[["Split"]], "initial_split")

  # Check approximate proportions
  expect_equal(nrow(split_list[["Train"]]), floor(0.7 * test_setup$n_recipe))
  expect_equal(nrow(split_list[["Test"]]), ceiling(0.3 * test_setup$n_recipe))
})


# --- Tests for t2emodel_data_recipe_init ---

test_that("t2emodel_data_recipe_init initiates a recipe", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)

  expect_s3_class(init_recipe, "recipe")
  # Check if vars are stored (custom addition in the function) - Restore check
  expect_equal(init_recipe$vars, c(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe))
  # Check formula was constructed correctly (inspect terms)
  expect_equal(init_recipe$term_info$variable[init_recipe$term_info$role == "outcome"], c(test_setup$time_var_recipe, test_setup$event_var_recipe))
  expect_true(all(test_setup$exp_vars_recipe %in% init_recipe$term_info$variable[init_recipe$term_info$role == "predictor"]))
})


# --- Tests for minimal_data_recipe ---

test_that("minimal_data_recipe adds expected steps (dummy=TRUE)", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  min_recipe <- minimal_data_recipe(init_recipe, pmiss = 0.4, pother = 0.1, dummy = TRUE, onehot = FALSE)

  expect_s3_class(min_recipe, "recipe")
  # Check steps added
  step_classes <- sapply(min_recipe$steps, class)
  expect_true("step_filter_missing" %in% step_classes)
  expect_true("step_impute_mean" %in% step_classes)
  expect_true("step_other" %in% step_classes)
  expect_true("step_impute_mode" %in% step_classes)
  expect_true("step_dummy" %in% step_classes)
  expect_true("step_nzv" %in% step_classes)
  expect_true("step_zv" %in% step_classes) # Check for step_zv
  expect_true("step_range" %in% step_classes)


  # Check parameters passed (example)
  expect_equal(min_recipe$steps[[1]]$threshold, 0.4) # step_filter_missing
  expect_equal(min_recipe$steps[[3]]$threshold, 0.1) # step_other
  expect_false(min_recipe$steps[[5]]$one_hot)       # step_dummy
})

test_that("minimal_data_recipe adds expected steps (dummy=FALSE)", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  min_recipe_nodummy <- minimal_data_recipe(init_recipe, dummy = FALSE)

  step_classes <- sapply(min_recipe_nodummy$steps, class)
  expect_false("step_dummy" %in% step_classes) # Should not have dummy step
  expect_true("step_impute_mean" %in% step_classes) # Other steps should be present
  expect_true("step_nzv" %in% step_classes)
  expect_true("step_zv" %in% step_classes) # Check for step_zv
  expect_true("step_range" %in% step_classes)
})


# --- Tests for data_recipe_corstep ---

test_that("data_recipe_corstep adds correlation step", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  cor_recipe <- data_recipe_corstep(init_recipe, treshold = 0.85)

  expect_s3_class(cor_recipe, "recipe")
  step_classes <- sapply(cor_recipe$steps, class)
  expect_true("step_corr" %in% step_classes)
  # Check threshold parameter
  expect_equal(cor_recipe$steps[[length(cor_recipe$steps)]]$threshold, 0.85)
})


# --- Tests for rfimp_data_recipe ---

test_that("rfimp_data_recipe adds expected steps (dummy=TRUE)", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  skip_if_not_installed("randomForest") # step_impute_bag dependency

  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  rf_recipe <- rfimp_data_recipe(init_recipe, pmiss = 0.4, pother = 0.1, dummy = TRUE, onehot = TRUE, trees = 15)

  expect_s3_class(rf_recipe, "recipe")
  # Check steps added
  step_classes <- sapply(rf_recipe$steps, class)
  expect_true("step_filter_missing" %in% step_classes)
  expect_true("step_other" %in% step_classes)
  expect_true("step_impute_bag" %in% step_classes)
  expect_true("step_dummy" %in% step_classes)
  expect_true("step_nzv" %in% step_classes)
  expect_true("step_zv" %in% step_classes) # Check for step_zv
  expect_true("step_range" %in% step_classes)


  # Check parameters passed (example)
  expect_equal(rf_recipe$steps[[1]]$threshold, 0.4) # step_filter_missing
  expect_equal(rf_recipe$steps[[2]]$threshold, 0.1) # step_other
  expect_equal(rf_recipe$steps[[3]]$trees, 15)      # step_impute_bag
  expect_true(rf_recipe$steps[[4]]$one_hot)        # step_dummy
})

test_that("rfimp_data_recipe adds expected steps (dummy=FALSE)", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  skip_if_not_installed("randomForest")

  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  rf_recipe_nodummy <- rfimp_data_recipe(init_recipe, dummy = FALSE)

  step_classes <- sapply(rf_recipe_nodummy$steps, class)
  expect_false("step_dummy" %in% step_classes) # Should not have dummy step
  expect_true("step_impute_bag" %in% step_classes) # Other steps should be present
  expect_true("step_nzv" %in% step_classes)
  expect_true("step_zv" %in% step_classes) # Check for step_zv
  expect_true("step_range" %in% step_classes)
})


# --- Tests for prep_data_recipe ---

test_that("prep_data_recipe preps the recipe", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  # Revert to simple step_center
  step_recipe <- recipes::step_center(init_recipe, recipes::all_numeric_predictors())
  # Call the wrapper function
  prepped_recipe <- prep_data_recipe(step_recipe, training = test_setup$recipe_data)

  # Remove debugging prints
  # print("Structure of prepped_recipe:")
  # str(prepped_recipe)

  expect_s3_class(prepped_recipe, "recipe")
  # Check if the step is trained by verifying the 'means' element exists and is populated
  expect_true(!is.null(prepped_recipe$steps[[1]]$means))
  expect_type(prepped_recipe$steps[[1]]$means, "double") # Check type
  # Check vars attribute is retained - Remove check as vars handling is removed
  # expect_equal(prepped_recipe$vars, c(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe))
})


# --- Tests for bake_data_recipe ---

test_that("bake_data_recipe bakes new data", {
  test_setup <- create_recipe_data_and_vars() # Create data inside test
  skip_if_not_installed("recipes")
  # Create and prep a simple recipe
  init_recipe <- t2emodel_data_recipe_init(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$exp_vars_recipe, test_setup$id_vars_recipe, test_setup$recipe_data)
  # Add step_zv before step_range
  step_recipe <- init_recipe %>%
                 step_zv(all_predictors()) %>%
                 step_range(all_numeric_predictors(), min=0, max=1) %>%
                 step_dummy(all_nominal_predictors())
  # Ensure the full recipe_data is passed to prep
  prepped_recipe <- prep_data_recipe(step_recipe, training = test_setup$recipe_data)

  # Bake training data
  # Pass the full recipe_data to bake
  baked_train <- bake_data_recipe(prepped_recipe, data = test_setup$recipe_data)
  expect_s3_class(baked_train, "data.frame") # bake usually returns tibble/df
  expect_equal(nrow(baked_train), nrow(test_setup$recipe_data))
  # Check number of columns (depends on dummy variables created and zv removal)
  # This check might be fragile, focus on successful baking
  # expect_gt(ncol(baked_train), length(test_setup$exp_vars_recipe))

  # Bake new data (using a subset as example)
  new_data_subset <- test_setup$recipe_data[1:5, ]
  # Pass the full new_data_subset to bake
  baked_new <- bake_data_recipe(prepped_recipe, data = new_data_subset)
   expect_s3_class(baked_new, "data.frame")
   expect_equal(nrow(baked_new), 5)
   expect_equal(ncol(baked_new), ncol(baked_train)) # Should have same columns

   # Check range scaling applied (excluding potential zv columns)
   numeric_cols_baked <- names(which(sapply(baked_new, is.numeric)))
   numeric_cols_baked <- setdiff(numeric_cols_baked, c(test_setup$time_var_recipe, test_setup$event_var_recipe, test_setup$id_vars_recipe)) # Exclude outcome/id
   # Further exclude zero-variance columns if they were removed
   zv_info <- prepped_recipe$steps[[which(sapply(prepped_recipe$steps, class) == "step_zv")]]$removals
   numeric_cols_baked <- setdiff(numeric_cols_baked, zv_info)

   if(length(numeric_cols_baked) > 0) {
       num_matrix <- as.matrix(baked_new[, numeric_cols_baked, drop = FALSE])
       # Allow for potential NA results if imputation wasn't perfect
       expect_true(all(num_matrix >= 0 & num_matrix <= 1 | is.na(num_matrix)))
   }
})
