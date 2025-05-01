library(testthat)
library(here) # Using here() for robust path construction

# Assuming the functions readData and readDict are available in the environment
# If they are part of a package, use library(yourpackagename)
source(here("R/data/data_io.R")) # Or source the file directly if not packaged

# Define paths to fixture files relative to the project root
fixture_data_path <- here("tests", "fixtures", "test_data.csv")
fixture_dict_path <- here("tests", "fixtures", "test_dict.csv")
non_existent_path <- here("tests", "fixtures", "non_existent_file.csv")

context("Testing data_io functions")

# --- Tests for readData ---

test_that("readData reads CSV correctly", {
  # Check if fixture file exists before running the test
  skip_if_not(file.exists(fixture_data_path), "Test data fixture not found")

  data <- readData(fixture_data_path)

  # Basic checks
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 5)
  expect_equal(ncol(data), 7)

  # Check column names
  expect_equal(colnames(data), c("ID", "Age", "Sex", "Value1", "Value2", "EventTime", "EventStatus"))

  # Check column types (adjust based on expected read behavior)
  expect_type(data$ID, "integer") # Assuming readr/data.table infers integer
  expect_type(data$Age, "integer")
  expect_type(data$Sex, "character")
  expect_type(data$Value1, "double")
  expect_type(data$Value2, "character")
  expect_type(data$EventTime, "integer")
  expect_type(data$EventStatus, "integer")

  # Check for NA value
  expect_true(is.na(data$Value1[3]))
})

test_that("readData handles non-existent file", {
  expect_error(readData(non_existent_path), regexp = "exist") # Check for error message containing "exist"
})

# --- Tests for readDict ---

test_that("readDict reads CSV correctly", {
  # Check if fixture file exists before running the test
  skip_if_not(file.exists(fixture_dict_path), "Test dictionary fixture not found")

  dict <- readDict(fixture_dict_path)

  # Basic checks
  expect_s3_class(dict, "data.frame")
  expect_equal(nrow(dict), 7)
  expect_equal(ncol(dict), 5)

  # Check column names
  expect_equal(colnames(dict), c("Variable", "Type", "Label", "MissingValue", "Format"))

  # Check column types (all should be character as read by default)
  expect_type(dict$Variable, "character")
  expect_type(dict$Type, "character")
  expect_type(dict$Label, "character")
  expect_type(dict$MissingValue, "character") # Even NAs are read as strings initially
  expect_type(dict$Format, "character")

  # Check specific values
  expect_equal(dict$Variable[1], "ID")
  expect_equal(dict$Type[2], "Num")
  expect_equal(dict$Label[3], "Biological Sex")
  expect_equal(dict$MissingValue[4], "NA") # readr reads NA as string "NA" by default
  expect_equal(dict$Format[4], "%.1f")
})

test_that("readDict handles non-existent file", {
  expect_error(readDict(non_existent_path), regexp = "exist") # Check for error message containing "exist"
})
