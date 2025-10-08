# Comprehensive Test Script for Data Functions in ml4time2event Package
# This script tests all functions in data.R, data_io.R, data_cleaning.R,
# data_transform.R, and data_summary.R

cat("==========================================================\n")
cat("Testing Data Functions in ml4time2event Package\n")
cat("==========================================================\n\n")

# Load the package from source
suppressMessages(devtools::load_all())

# Set seed for reproducibility
set.seed(123)

# Create test datasets
create_test_data <- function() {
  data.frame(
    id = 1:100,
    age = c(rnorm(95, 50, 10), rep(NA, 5)),
    height = c(rnorm(90, 170, 15), rep(200, 5), rep(NA, 5)),
    weight = c(rnorm(85, 70, 12), rep(150, 10), rep(NA, 5)),
    gender = factor(c(rep("M", 45), rep("F", 45), rep(NA, 10))),
    race = factor(c(rep("A", 30), rep("B", 30), rep("C", 30), rep("D", 2), rep(NA, 8))),
    status = c(rep(0, 50), rep(1, 48), rep(NA, 2)),
    constant_var = rep(5, 100),
    all_na_var = rep(NA_real_, 100),
    binary_var = c(rep(0, 98), rep(1, 2)),
    char_col = c(rep("  spaces  ", 20), rep("normal", 70), rep(NA_character_, 10)),
    stringsAsFactors = FALSE
  )
}

test_data <- create_test_data()

cat("\n--- Test Data Created ---\n")
cat("Dimensions:", dim(test_data), "\n")
cat("Columns:", paste(colnames(test_data), collapse = ", "), "\n\n")

#==============================================================================
# TEST 1: data_cleaning.R functions
#==============================================================================
cat("\n### TEST 1: RemoveMonoVarsData ###\n")
tryCatch({
  result <- RemoveMonoVarsData(test_data)
  cat("PASS: RemoveMonoVarsData executed successfully\n")
  cat("Removed columns:", setdiff(colnames(test_data), colnames(result)), "\n")
  cat("Remaining columns:", ncol(result), "\n")
}, error = function(e) {
  cat("FAIL: RemoveMonoVarsData -", e$message, "\n")
})

cat("\n### TEST 2: RemoveAllNAVars ###\n")
tryCatch({
  result <- RemoveAllNAVars(test_data)
  cat("PASS: RemoveAllNAVars executed successfully\n")
  cat("Removed columns:", setdiff(colnames(test_data), colnames(result)), "\n")
  cat("Remaining columns:", ncol(result), "\n")
}, error = function(e) {
  cat("FAIL: RemoveAllNAVars -", e$message, "\n")
})

cat("\n### TEST 3: RemoveMissinVarsData ###\n")
tryCatch({
  result <- RemoveMissinVarsData(test_data, maxprop = 0.05)
  cat("PASS: RemoveMissinVarsData executed successfully\n")
  cat("Removed columns:", setdiff(colnames(test_data), colnames(result)), "\n")
  cat("Remaining columns:", ncol(result), "\n")
}, error = function(e) {
  cat("FAIL: RemoveMissinVarsData -", e$message, "\n")
})

cat("\n### TEST 4: RemoveMissinRecordsData ###\n")
tryCatch({
  result <- RemoveMissinRecordsData(test_data, maxprop = 0.2)
  cat("PASS: RemoveMissinRecordsData executed successfully\n")
  cat("Removed rows:", nrow(test_data) - nrow(result), "\n")
  cat("Remaining rows:", nrow(result), "\n")
}, error = function(e) {
  cat("FAIL: RemoveMissinRecordsData -", e$message, "\n")
})

cat("\n### TEST 5: getcharcolsData ###\n")
tryCatch({
  result <- getcharcolsData(test_data)
  cat("PASS: getcharcolsData executed successfully\n")
  cat("Character columns:", paste(result, collapse = ", "), "\n")
}, error = function(e) {
  cat("FAIL: getcharcolsData -", e$message, "\n")
})

cat("\n### TEST 6: RemoveRareCategoriesData ###\n")
tryCatch({
  result <- RemoveRareCategoriesData(test_data, minfreq = 0.05)
  cat("PASS: RemoveRareCategoriesData executed successfully\n")
}, error = function(e) {
  cat("FAIL: RemoveRareCategoriesData -", e$message, "\n")
})

cat("\n### TEST 7: RemoveRareBinaryVarsData ###\n")
tryCatch({
  result <- RemoveRareBinaryVarsData(test_data, minfreq = 0.05)
  cat("PASS: RemoveRareBinaryVarsData executed successfully\n")
}, error = function(e) {
  cat("FAIL: RemoveRareBinaryVarsData -", e$message, "\n")
})

cat("\n### TEST 8: CollapseRareCategoriestoOtherData ###\n")
tryCatch({
  result <- CollapseRareCategoriestoOtherData(test_data, minfreq = 0.05)
  cat("PASS: CollapseRareCategoriestoOtherData executed successfully\n")
}, error = function(e) {
  cat("FAIL: CollapseRareCategoriestoOtherData -", e$message, "\n")
})

cat("\n### TEST 9: droplevelsoffactorsData ###\n")
tryCatch({
  test_data_copy <- test_data
  test_data_copy$gender <- factor(test_data_copy$gender, levels = c("M", "F", "X"))
  result <- droplevelsoffactorsData(test_data_copy)
  cat("PASS: droplevelsoffactorsData executed successfully\n")
}, error = function(e) {
  cat("FAIL: droplevelsoffactorsData -", e$message, "\n")
})

cat("\n### TEST 10: findvarsnamesthatrepeatData ###\n")
tryCatch({
  test_data_names <- test_data
  colnames(test_data_names)[1:3] <- c("age", "age_new", "age_old")
  result <- findvarsnamesthatrepeatData(test_data_names)
  cat("PASS: findvarsnamesthatrepeatData executed successfully\n")
}, error = function(e) {
  cat("FAIL: findvarsnamesthatrepeatData -", e$message, "\n")
})

cat("\n### TEST 11: ReplaceOutlierNumValsData ###\n")
tryCatch({
  result <- ReplaceOutlierNumValsData(test_data, multIQR = 1.5, minnumgroup = 10)
  cat("PASS: ReplaceOutlierNumValsData executed successfully\n")
}, error = function(e) {
  cat("FAIL: ReplaceOutlierNumValsData -", e$message, "\n")
})

cat("\n### TEST 12: MakeTestDataConfWithTrainData ###\n")
tryCatch({
  train_data <- test_data[1:70, ]
  test_data_2 <- test_data[71:100, ]
  test_data_2$extra_col <- rnorm(30)
  result <- MakeTestDataConfWithTrainData(train_data, test_data_2)
  cat("PASS: MakeTestDataConfWithTrainData executed successfully\n")
  cat("Test data conformed to", ncol(result), "columns\n")
}, error = function(e) {
  cat("FAIL: MakeTestDataConfWithTrainData -", e$message, "\n")
})

cat("\n### TEST 13: RemoveEmptySpacesData ###\n")
tryCatch({
  result <- RemoveEmptySpacesData(test_data)
  cat("PASS: RemoveEmptySpacesData executed successfully\n")
}, error = function(e) {
  cat("FAIL: RemoveEmptySpacesData -", e$message, "\n")
})

#==============================================================================
# TEST 2: data_transform.R functions
#==============================================================================
cat("\n\n### TEST 14: ZeroOneScalerData ###\n")
tryCatch({
  result <- ZeroOneScalerData(test_data)
  cat("PASS: ZeroOneScalerData executed successfully\n")
  cat("Returned list components:", paste(names(result), collapse = ", "), "\n")
}, error = function(e) {
  cat("FAIL: ZeroOneScalerData -", e$message, "\n")
})

cat("\n### TEST 15: ZeroOneScalerApplierData ###\n")
tryCatch({
  scaler_result <- ZeroOneScalerData(test_data)
  new_data <- create_test_data()
  result <- ZeroOneScalerApplierData(new_data,
                                     scaler_result$minxvec,
                                     scaler_result$maxxvec)
  cat("PASS: ZeroOneScalerApplierData executed successfully\n")
}, error = function(e) {
  cat("FAIL: ZeroOneScalerApplierData -", e$message, "\n")
})

cat("\n### TEST 16: UndoZeroOneScalerApplierData ###\n")
tryCatch({
  scaler_result <- ZeroOneScalerData(test_data)
  scaled_data <- scaler_result$data
  result <- UndoZeroOneScalerApplierData(scaled_data,
                                         scaler_result$minxvec,
                                         scaler_result$maxxvec)
  cat("PASS: UndoZeroOneScalerApplierData executed successfully\n")
}, error = function(e) {
  cat("FAIL: UndoZeroOneScalerApplierData -", e$message, "\n")
})

cat("\n### TEST 17: NumVarstCatsData with numgroups ###\n")
tryCatch({
  result <- NumVarstCatsData(test_data, numgroups = 4)
  cat("PASS: NumVarstCatsData (numgroups) executed successfully\n")
}, error = function(e) {
  cat("FAIL: NumVarstCatsData (numgroups) -", e$message, "\n")
})

cat("\n### TEST 18: NumVarstCatsData with cuts ###\n")
tryCatch({
  result <- NumVarstCatsData(test_data, cuts = c(50, 100, 150))
  cat("PASS: NumVarstCatsData (cuts) executed successfully\n")
}, error = function(e) {
  cat("FAIL: NumVarstCatsData (cuts) -", e$message, "\n")
})

#==============================================================================
# TEST 3: data_summary.R functions
#==============================================================================
cat("\n\n### TEST 19: pairwiserelationshipsDataSummmary ###\n")
tryCatch({
  small_data <- test_data[1:50, c("age", "height", "weight", "status")]
  small_data <- na.omit(small_data)
  result <- pairwiserelationshipsDataSummmary(small_data)
  cat("PASS: pairwiserelationshipsDataSummmary executed successfully\n")
  cat("Correlation matrix dimensions:", dim(result), "\n")
}, error = function(e) {
  cat("FAIL: pairwiserelationshipsDataSummmary -", e$message, "\n")
})

cat("\n### TEST 20: gethighcorvarsDataSummmary ###\n")
tryCatch({
  small_data <- test_data[1:50, c("age", "height", "weight", "status")]
  small_data <- na.omit(small_data)
  cormat <- pairwiserelationshipsDataSummmary(small_data)
  result <- gethighcorvarsDataSummmary(cormat, corcutoff = 0.3)
  cat("PASS: gethighcorvarsDataSummmary executed successfully\n")
}, error = function(e) {
  cat("FAIL: gethighcorvarsDataSummmary -", e$message, "\n")
})

cat("\n### TEST 21: OneAgainstRestCorDataSummmary ###\n")
tryCatch({
  small_data <- test_data[1:50, c("age", "height", "weight", "status")]
  small_data <- na.omit(small_data)
  result <- OneAgainstRestCorDataSummmary(small_data)
  cat("PASS: OneAgainstRestCorDataSummmary executed successfully\n")
  cat("Correlation vector length:", length(result), "\n")
}, error = function(e) {
  cat("FAIL: OneAgainstRestCorDataSummmary -", e$message, "\n")
})

cat("\n### TEST 22: SummaryTableDataSummmary ###\n")
tryCatch({
  if (requireNamespace("gtsummary", quietly = TRUE)) {
    result <- SummaryTableDataSummmary(test_data,
                                       UseVars = c("age", "height", "gender"))
    cat("PASS: SummaryTableDataSummmary executed successfully\n")
  } else {
    cat("SKIP: SummaryTableDataSummmary - gtsummary not installed\n")
  }
}, error = function(e) {
  cat("FAIL: SummaryTableDataSummmary -", e$message, "\n")
})

#==============================================================================
# TEST 4: data_io.R functions
#==============================================================================
cat("\n\n### TEST 23: readData (CSV) ###\n")
tryCatch({
  # Create a temporary CSV file
  temp_csv <- tempfile(fileext = ".csv")
  write.csv(test_data, temp_csv, row.names = FALSE)
  result <- readData(file = temp_csv, filetype = "csv")
  unlink(temp_csv)
  cat("PASS: readData (CSV) executed successfully\n")
  cat("Read", nrow(result), "rows and", ncol(result), "columns\n")
}, error = function(e) {
  cat("FAIL: readData (CSV) -", e$message, "\n")
})

cat("\n### TEST 24: readData (XLSX) ###\n")
tryCatch({
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    temp_xlsx <- tempfile(fileext = ".xlsx")
    openxlsx::write.xlsx(test_data, temp_xlsx)
    result <- readData(file = temp_xlsx, filetype = "xlsx")
    unlink(temp_xlsx)
    cat("PASS: readData (XLSX) executed successfully\n")
    cat("Read", nrow(result), "rows and", ncol(result), "columns\n")
  } else {
    cat("SKIP: readData (XLSX) - openxlsx not installed\n")
  }
}, error = function(e) {
  cat("FAIL: readData (XLSX) -", e$message, "\n")
})

#==============================================================================
# TEST 5: Edge Cases
#==============================================================================
cat("\n\n### TEST 25: Edge Case - Empty Data Frame ###\n")
tryCatch({
  empty_df <- data.frame()
  result1 <- RemoveMonoVarsData(empty_df)
  result2 <- RemoveAllNAVars(empty_df)
  result3 <- getcharcolsData(empty_df)
  cat("PASS: Empty data frame handling\n")
}, error = function(e) {
  cat("FAIL: Empty data frame -", e$message, "\n")
})

cat("\n### TEST 26: Edge Case - All NA Column ###\n")
tryCatch({
  all_na_df <- data.frame(
    x = rep(NA_real_, 10),
    y = 1:10
  )
  result <- RemoveAllNAVars(all_na_df)
  cat("PASS: All NA column handling\n")
  cat("Remaining columns:", ncol(result), "\n")
}, error = function(e) {
  cat("FAIL: All NA column -", e$message, "\n")
})

cat("\n### TEST 27: Edge Case - Single Value Column ###\n")
tryCatch({
  single_val_df <- data.frame(
    x = rep(5, 10),
    y = 1:10
  )
  result <- RemoveMonoVarsData(single_val_df)
  cat("PASS: Single value column handling\n")
  cat("Remaining columns:", ncol(result), "\n")
}, error = function(e) {
  cat("FAIL: Single value column -", e$message, "\n")
})

#==============================================================================
# TEST 6: ImputeMissinRecordsData (requires missRanger)
#==============================================================================
cat("\n### TEST 28: ImputeMissinRecordsData ###\n")
tryCatch({
  if (requireNamespace("missRanger", quietly = TRUE)) {
    # Create data with missing values
    impute_data <- test_data[, c("age", "height", "weight", "gender")]
    result <- ImputeMissinRecordsData(impute_data, dontuse = "gender")
    cat("PASS: ImputeMissinRecordsData executed successfully\n")
    cat("NA count before:", sum(is.na(impute_data)), "\n")
    cat("NA count after:", sum(is.na(result)), "\n")
  } else {
    cat("SKIP: ImputeMissinRecordsData - missRanger not installed\n")
  }
}, error = function(e) {
  cat("FAIL: ImputeMissinRecordsData -", e$message, "\n")
})

#==============================================================================
# SUMMARY
#==============================================================================
cat("\n\n==========================================================\n")
cat("Testing Complete!\n")
cat("==========================================================\n")
cat("\nPlease review the output above for any FAIL messages.\n")
cat("All tests should show PASS or SKIP (if optional dependencies are missing).\n\n")
