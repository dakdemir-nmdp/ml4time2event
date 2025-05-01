library(testthat)
library(here)
# library(data.table) # Removed

# Assuming the functions are available in the environment
source(here("R/data/data_transform.R"))

context("Testing data_transform functions")

# --- Test Data Setup ---
test_df_transform <- data.frame( # Renamed and replaced data.table()
  id = 1:5,
  num1 = c(10, 20, 30, 40, 50),
  num2 = c(-5, 0, 5, 10, 15),
  num_const = rep(5, 5),
  num_na = c(1, 2, NA, 4, 5),
  char_col = letters[1:5],
  factor_col = factor(letters[1:5]),
  stringsAsFactors = FALSE
)

# --- Tests for ZeroOneScalerData ---

test_that("ZeroOneScalerData scales numeric columns to [0, 1]", {
  scaled_result <- ZeroOneScalerData(test_df_transform) # Removed copy()
  scaled_data <- scaled_result$data

  # Check structure of result
  expect_type(scaled_result, "list")
  expect_named(scaled_result, c("data", "minxvec", "maxxvec"))
  expect_s3_class(scaled_data, "data.frame") # Check for data.frame

  # Check scaled values for num1 (range 10-50)
  expect_equal(scaled_data$num1, c(0, 0.25, 0.5, 0.75, 1))
  expect_equal(scaled_result$minxvec["num1"], 10)
  expect_equal(scaled_result$maxxvec["num1"], 50)

  # Check scaled values for num2 (range -5-15)
  expect_equal(scaled_data$num2, c(0, 0.25, 0.5, 0.75, 1))
  expect_equal(scaled_result$minxvec["num2"], -5)
  expect_equal(scaled_result$maxxvec["num2"], 15)

  # Check constant column (should be scaled to 0 with a warning)
  expect_warning(scaled_result_warn <- ZeroOneScalerData(test_df_transform), "zero range") # Removed copy()
  expect_equal(scaled_result_warn$data$num_const, rep(0, 5))
  expect_equal(scaled_result_warn$minxvec["num_const"], 5)
  expect_equal(scaled_result_warn$maxxvec["num_const"], 5)

  # Check column with NA
  expect_equal(scaled_data$num_na[c(1, 2, 4, 5)], c(0, 0.25, 0.75, 1))
  expect_true(is.na(scaled_data$num_na[3]))
  expect_equal(scaled_result$minxvec["num_na"], 1)
  expect_equal(scaled_result$maxxvec["num_na"], 5)

  # Check non-numeric columns are untouched
  expect_equal(scaled_data$id, test_df_transform$id)
  expect_equal(scaled_data$char_col, test_df_transform$char_col)
  expect_equal(scaled_data$factor_col, test_df_transform$factor_col)

  # Check min/max vectors have correct names and include NAs for non-numeric
  expect_equal(names(scaled_result$minxvec), colnames(test_df_transform))
  expect_true(is.na(scaled_result$minxvec["id"]))
  expect_true(is.na(scaled_result$minxvec["char_col"]))
  expect_true(is.na(scaled_result$minxvec["factor_col"]))
})

test_that("ZeroOneScalerData handles data.frames", {
  test_df <- as.data.frame(test_df_transform)
  scaled_result <- ZeroOneScalerData(test_df)
  expect_s3_class(scaled_result$data, "data.frame") # Already checking data.frame
  expect_equal(scaled_result$data$num1, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("ZeroOneScalerData handles empty input", {
  # empty_dt <- data.table() # Removed data.table call
  # scaled_result <- ZeroOneScalerData(empty_dt) # Removed data.table test
  # expect_equal(scaled_result$data, empty_dt) # Removed data.table test
  # expect_length(scaled_result$minxvec, 0) # Removed data.table test
  # expect_length(scaled_result$maxxvec, 0) # Removed data.table test

  empty_df <- data.frame()
  scaled_result_df <- ZeroOneScalerData(empty_df)
  expect_equal(scaled_result_df$data, empty_df)
  expect_length(scaled_result_df$minxvec, 0)
  expect_length(scaled_result_df$maxxvec, 0)
})


# --- Tests for ZeroOneScalerApplierData ---

# Get scaling parameters from original data
scaling_params <- ZeroOneScalerData(test_df_transform) # Removed copy()
mins <- scaling_params$minxvec
maxs <- scaling_params$maxxvec

# New data to apply scaling to
new_data <- data.frame( # Replaced data.table()
  id = 6:7,
  num1 = c(15, 35), # Within original range [10, 50]
  num2 = c(-10, 20), # Outside original range [-5, 15]
  num_const = c(5, 5), # Constant value
  num_na = c(NA, 3), # Has NA, value within range [1, 5]
  char_col = c("f", "g"),
  factor_col = factor(c("f", "g")),
  extra_col = c(100, 200), # Not in original scaling params
  stringsAsFactors = FALSE
)

test_that("ZeroOneScalerApplierData applies scaling correctly", {
  applied_data <- ZeroOneScalerApplierData(new_data, mins, maxs) # Removed copy()

  # Check num1: (15-10)/(50-10)=0.125, (35-10)/(50-10)=0.625
  expect_equal(applied_data$num1, c(0.125, 0.625))

  # Check num2: (-10 - -5)/(15 - -5) = -0.25, (20 - -5)/(15 - -5) = 1.25
  # Note: Function doesn't cap by default
  expect_equal(applied_data$num2, c(-0.25, 1.25))

  # Check num_const: range is 0, should become 0
  expect_equal(applied_data$num_const, c(0, 0))

  # Check num_na: NA remains NA, (3-1)/(5-1) = 0.5
  expect_true(is.na(applied_data$num_na[1]))
  expect_equal(applied_data$num_na[2], 0.5)

  # Check untouched columns
  expect_equal(applied_data$id, new_data$id)
  expect_equal(applied_data$char_col, new_data$char_col)
  expect_equal(applied_data$factor_col, new_data$factor_col)

  # Check column not in scaling params (should be skipped with warning)
  expect_warning(applied_data_warn <- ZeroOneScalerApplierData(new_data, mins, maxs), # Removed copy()
                 "not provided for numeric column 'extra_col'")
  expect_equal(applied_data_warn$extra_col, new_data$extra_col)
})

test_that("ZeroOneScalerApplierData handles missing/NA scaling parameters", {
  mins_missing <- mins[!names(mins) %in% "num1"]
  maxs_missing <- maxs[!names(maxs) %in% "num1"]
  expect_warning(applied_data <- ZeroOneScalerApplierData(new_data, mins_missing, maxs), # Removed copy()
                 "not provided for numeric column 'num1'")
  expect_equal(applied_data$num1, new_data$num1) # Should be unchanged

  mins_na <- mins
  mins_na["num1"] <- NA
  expect_warning(applied_data_na <- ZeroOneScalerApplierData(new_data, mins_na, maxs), # Removed copy()
                 "Missing min/max value for numeric column 'num1'")
  expect_equal(applied_data_na$num1, new_data$num1) # Should be unchanged
})

test_that("ZeroOneScalerApplierData handles empty input", {
   # empty_dt <- data.table() # Removed data.table call
   # expect_equal(ZeroOneScalerApplierData(empty_dt, mins, maxs), empty_dt) # Removed data.table test
   empty_df <- data.frame()
   expect_equal(ZeroOneScalerApplierData(empty_df, mins, maxs), empty_df)
})


# --- Tests for UndoZeroOneScalerApplierData ---

# Use the scaled data from the first test
scaled_data_to_undo <- scaling_params$data

test_that("UndoZeroOneScalerApplierData reverses scaling correctly", {
  undone_data <- UndoZeroOneScalerApplierData(scaled_data_to_undo, mins, maxs) # Removed copy()

  # Check numeric columns are restored (allowing for float precision issues)
  expect_equal(undone_data$num1, test_df_transform$num1)
  expect_equal(undone_data$num2, test_df_transform$num2)
  expect_equal(undone_data$num_const, test_df_transform$num_const) # Should revert to original constant
  expect_equal(undone_data$num_na, test_df_transform$num_na) # NAs should match

  # Check non-numeric columns are untouched
  expect_equal(undone_data$id, test_df_transform$id)
  expect_equal(undone_data$char_col, test_df_transform$char_col)
  expect_equal(undone_data$factor_col, test_df_transform$factor_col)
})

test_that("UndoZeroOneScalerApplierData handles missing/NA scaling parameters", {
  mins_missing <- mins[!names(mins) %in% "num1"]
  expect_warning(undone_data <- UndoZeroOneScalerApplierData(scaled_data_to_undo, mins_missing, maxs), # Removed copy()
                 "not provided for numeric column 'num1'")
  expect_equal(undone_data$num1, scaled_data_to_undo$num1) # Should be unchanged (still scaled)

  maxs_na <- maxs
  maxs_na["num2"] <- NA
  expect_warning(undone_data_na <- UndoZeroOneScalerApplierData(scaled_data_to_undo, mins, maxs_na), # Removed copy()
                 "Missing min/max value for numeric column 'num2'")
  expect_equal(undone_data_na$num2, scaled_data_to_undo$num2) # Should be unchanged
})

test_that("UndoZeroOneScalerApplierData handles empty input", {
   # empty_dt <- data.table() # Removed data.table call
   # expect_equal(UndoZeroOneScalerApplierData(empty_dt, mins, maxs), empty_dt) # Removed data.table test
   empty_df <- data.frame()
   expect_equal(UndoZeroOneScalerApplierData(empty_df, mins, maxs), empty_df)
})


# --- Tests for NumVarstCatsData ---

test_df_categorize <- data.frame( # Replaced data.table()
  id = 1:20,
  numeric_wide = seq(1, 100, length.out = 20), # Wide range
  numeric_narrow = rnorm(20, 5, 0.5),        # Narrow range
  numeric_skewed = c(rep(1, 15), 10, 20, 30, 40, 50), # Skewed
  numeric_few_unique = rep(c(1, 5, 10), length.out = 20), # Only 3 unique
  numeric_with_na = c(1:5, NA, 7:10, NA, 12:20),
  char_col = letters[1:20],
  stringsAsFactors = FALSE
)

test_that("NumVarstCatsData categorizes using numgroups (quantiles)", {
  # Categorize into 4 groups (quartiles)
  categorized_data <- NumVarstCatsData(test_df_categorize, numgroups = 4, min_unique_vals = 5) # Removed copy()

  # Check numeric_wide
  expect_s3_class(categorized_data$numeric_wide, "factor")
  expect_length(levels(categorized_data$numeric_wide), 4) # Expect 4 levels for quartiles
  # Check counts per level (should be roughly equal for uniform data)
  expect_equal(as.numeric(table(categorized_data$numeric_wide)), rep(5, 4))

  # Check numeric_narrow
  expect_s3_class(categorized_data$numeric_narrow, "factor")
  expect_length(levels(categorized_data$numeric_narrow), 4)

  # Check numeric_skewed (quantiles might be duplicated)
  expect_s3_class(categorized_data$numeric_skewed, "factor")
  # The number of levels might be less than 4 if quantiles are duplicated
  expect_true(length(levels(categorized_data$numeric_skewed)) <= 4)
  # Check that the high values are in the last category
  last_level <- levels(categorized_data$numeric_skewed)[length(levels(categorized_data$numeric_skewed))]
  expect_true(all(categorized_data$numeric_skewed[16:20] %in% last_level))


  # Check numeric_few_unique (should be skipped)
  expect_type(categorized_data$numeric_few_unique, "double") # Remains numeric

  # Check numeric_with_na
  expect_s3_class(categorized_data$numeric_with_na, "factor")
  expect_length(levels(categorized_data$numeric_with_na), 4)
  expect_equal(sum(is.na(categorized_data$numeric_with_na)), 2) # NAs preserved

  # Check non-numeric untouched
  expect_equal(categorized_data$char_col, test_df_categorize$char_col)
})

test_that("NumVarstCatsData categorizes using specific cuts", {
  cuts_defined <- c(20, 50, 80)
  categorized_data <- NumVarstCatsData(test_df_categorize, cuts = cuts_defined, min_unique_vals = 5) # Removed copy()

  # Check numeric_wide (1-100)
  expect_s3_class(categorized_data$numeric_wide, "factor")
  # Levels should be like [-Inf, 20), [20, 50), [50, 80), [80, Inf)
  expect_length(levels(categorized_data$numeric_wide), 4)
  # Check a value falls in the correct bin, e.g., value 30 should be in [20, 50)
  expect_true(grepl("\\[20, 50\\)", categorized_data$numeric_wide[test_df_categorize$numeric_wide > 20 & test_df_categorize$numeric_wide < 50][1]))

  # Check numeric_narrow (all values likely < 20)
  expect_s3_class(categorized_data$numeric_narrow, "factor")
  expect_length(levels(categorized_data$numeric_narrow), 4) # Still creates all levels based on cuts
  expect_true(all(categorized_data$numeric_narrow %in% levels(categorized_data$numeric_narrow)[1])) # All should be in first bin

  # Check numeric_few_unique (skipped)
  expect_type(categorized_data$numeric_few_unique, "double")
})

test_that("NumVarstCatsData respects min_unique_vals", {
  # Use high min_unique_vals - should skip most columns
  categorized_data <- NumVarstCatsData(test_df_categorize, numgroups = 4, min_unique_vals = 15) # Removed copy()

  expect_s3_class(categorized_data$numeric_wide, "factor") # Has 20 unique vals
  expect_s3_class(categorized_data$numeric_narrow, "factor") # Likely has >= 15 unique vals
  expect_type(categorized_data$numeric_skewed, "double") # Only 6 unique vals
  expect_type(categorized_data$numeric_few_unique, "double") # Only 3 unique vals
  expect_s3_class(categorized_data$numeric_with_na, "factor") # 18 unique non-NA vals
})

test_that("NumVarstCatsData handles duplicated quantiles", {
  # Data where quantiles will likely be duplicated
  df_dup_quant <- data.frame(x = c(rep(1, 10), rep(2, 5), rep(10, 5)), stringsAsFactors = FALSE) # Replaced data.table()
  # Quartiles (0, 0.25, 0.5, 0.75, 1) might be (1, 1, 1, 2, 10)
  expect_warning(categorized_data <- NumVarstCatsData(df_dup_quant, numgroups = 4, min_unique_vals = 3)) # Removed copy()
  expect_s3_class(categorized_data$x, "factor")
  # Expect fewer than 4 levels due to duplicate quantiles being treated as single breaks
  expect_equal(length(levels(categorized_data$x)), 3) # Expect levels like [1, 1], (1, 2], (2, 10] -> simplified to [1, 2), [2, 10) ? Check function logic.
  # Based on the internal quantcat logic, it uses unique quantiles as breaks: 1, 2, 10
  # Breaks become: 1, 2, 10 -> Levels: [1, 2), [2, 10]
  expect_equal(levels(categorized_data$x), c("[1, 2)", "[2, 10]"))
})


test_that("NumVarstCatsData requires numgroups or cuts", {
  expect_error(NumVarstCatsData(test_df_categorize), "Must provide either 'numgroups' or 'cuts'")
})

test_that("NumVarstCatsData warns if both numgroups and cuts provided", {
   expect_warning(NumVarstCatsData(test_df_categorize, numgroups = 4, cuts = c(25, 75)),
                  "'cuts' provided, 'numgroups' will be ignored")
   # Check that cuts were actually used
   categorized_data <- suppressWarnings(NumVarstCatsData(test_df_categorize, numgroups = 4, cuts = c(25, 75))) # Removed copy()
   expect_length(levels(categorized_data$numeric_wide), 3) # Cuts c(25, 75) -> 3 levels
})

test_that("NumVarstCatsData handles empty input", {
   # empty_dt <- data.table() # Removed data.table call
   # expect_equal(NumVarstCatsData(empty_dt, numgroups = 4), empty_dt) # Removed data.table test
   empty_df <- data.frame()
   expect_equal(NumVarstCatsData(empty_df, cuts = c(10, 20)), empty_df)
})
