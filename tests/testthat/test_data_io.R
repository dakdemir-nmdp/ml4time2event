library(testthat)
library(here) # Using here() for robust path construction

# Source the data_io functions
source(here("R/data_io.R"))

# Define paths to fixture files relative to the project root
fixture_data_path <- here("tests", "fixtures", "test_data.csv")
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

  # Check column names (readData converts to lowercase)
  expect_equal(colnames(data), c("id", "age", "sex", "value1", "value2", "eventtime", "eventstatus"))

  # Check column types (adjust based on expected read behavior)
  expect_type(data$id, "integer")
  expect_type(data$age, "integer")
  expect_type(data$sex, "character")
  expect_type(data$value1, "double")
  expect_type(data$value2, "character")
  expect_type(data$eventtime, "integer")
  expect_type(data$eventstatus, "integer")

  # Check for NA value
  expect_true(is.na(data$value1[3]))
})

test_that("readData handles non-existent file", {
  expect_error(readData(non_existent_path), regexp = "exist") # Check for error message containing "exist"
})
