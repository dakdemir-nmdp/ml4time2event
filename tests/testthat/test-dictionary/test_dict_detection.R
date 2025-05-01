library(testthat)
library(here)
# library(data.table) # Removed

# Assuming the functions are available in the environment
source(here("R/dictionary/dict_detection.R"))

context("Testing dict_detection functions")

# --- Test Data Setup ---
test_dict_detection <- data.frame(
  Variable = c("PatientID", "Study_ID", "Record_Identifier", "Age", "Sex", "TreatmentDate", "Outcome_DTE", "Status", "Score", "Notes", "Pseudo_Code", "CaseNum", "DT_Start", "End_DT"),
  Label = c("Patient Identifier", "Study ID", "Record ID", "Age in Years", "Biological Sex", "Date of Treatment", "Outcome Date", "Patient Status", "Risk Score", "Clinical Notes", "Pseudonymized Code", "Case Number", "Start Date", "End Date"),
  Value = c(NA, NA, NA, NA, "M", NA, NA, "Unknown", "999", NA, NA, NA, NA, NA),
  Type = c("Num", "Char", "Char", "Num", "Cat", "Date", "Date", "Cat", "Num", "Char", "Char", "Num", "Date", "Date"),
  stringsAsFactors = FALSE
)

test_dict_na <- data.frame(
    Variable = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"),
    Label    = c("Regular Value", "Missing Data", "Unknown Status", "Not Answered", "Valid Label", "Valid Label", "Valid Label", "Valid Label", "Valid Label", "Valid Label", "Refused Response", "N/A Value", "Valid Label", "Declined"),
    Value    = c("1", "Missing", "Unknown", "Not Answ", ".", ".b", "99", "888", "<NA>", "", "Refused", "N/A", "Valid", "Declined"),
    Type     = rep("Char", 14),
    stringsAsFactors = FALSE
)


# --- Tests for DetectIdVarsDict ---

test_that("DetectIdVarsDict identifies ID variables by pattern", {
  detected_ids <- DetectIdVarsDict(test_dict_detection)
  expect_true(all(c("PatientID", "Study_ID", "Record_Identifier", "Pseudo_Code", "CaseNum") %in% detected_ids))
  expect_false("Age" %in% detected_ids)
  expect_false("TreatmentDate" %in% detected_ids)
})

test_that("DetectIdVarsDict includes user-provided idvars", {
  detected_ids <- DetectIdVarsDict(test_dict_detection, idvars = c("Age", "Status")) # Add non-ID vars for testing
  expect_true(all(c("PatientID", "Study_ID", "Record_Identifier", "Pseudo_Code", "CaseNum", "Age", "Status") %in% detected_ids))
})

test_that("DetectIdVarsDict handles case-insensitivity", {
  dict_case <- data.frame(Variable = c("patientid", "age", "NOTES"), stringsAsFactors = FALSE)
  detected_ids <- DetectIdVarsDict(dict_case)
  expect_equal(detected_ids, "patientid")
})

test_that("DetectIdVarsDict handles no ID variables found", {
  dict_no_ids <- data.frame(Variable = c("Age", "Score", "Notes"), stringsAsFactors = FALSE)
  detected_ids <- DetectIdVarsDict(dict_no_ids)
  expect_length(detected_ids, 0)
  expect_type(detected_ids, "character")
})

test_that("DetectIdVarsDict handles empty dictionary", {
  empty_dict <- data.frame(Variable = character(0), stringsAsFactors = FALSE)
  expect_length(DetectIdVarsDict(empty_dict), 0)
  empty_dict_no_col <- data.frame(Label = character(0), stringsAsFactors = FALSE)
  expect_error(DetectIdVarsDict(empty_dict_no_col), "must contain a 'Variable' column")
})


# --- Tests for DetectNAStringsDict ---

test_that("DetectNAStringsDict identifies NA strings and adds NACats column", {
  dict_na_detected <- DetectNAStringsDict(test_dict_na) # Use standard assignment

  expect_true("NACats" %in% colnames(dict_na_detected))
  expect_equal(nrow(dict_na_detected), nrow(test_dict_na))

  # Check specific rows based on default NA patterns
  # Labels: "Missing Data", "Unknown Status", "Not Answered", "Refused Response", "N/A Value", "Declined"
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "Missing Data"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "Unknown Status"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "Not Answered"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "Refused Response"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "N/A Value"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Label == "Declined"], "NAVAL")

  # Values: ".", ".b", "99", "888", "<NA>", ""
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == "."], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == ".b"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == "99"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == "888"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == "<NA>"], "NAVAL")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == ""], "NAVAL")

  # Check a regular value
  expect_equal(dict_na_detected$NACats[dict_na_detected$Variable == "Var1"], "VALUE")
  expect_equal(dict_na_detected$NACats[dict_na_detected$Value == "Valid"], "VALUE")

  # Count total NAVALs expected
  expected_naval_count <- length(unique(c(
      which(test_dict_na$Label %in% c("Missing Data", "Unknown Status", "Not Answered", "Refused Response", "N/A Value", "Declined")),
      which(test_dict_na$Value %in% c(".", ".b", "99", "888", "<NA>", "", "Missing", "Unknown", "Not Answ", "Refused", "N/A", "Declined")) # Values can also match label patterns
  )))
   expect_equal(sum(dict_na_detected$NACats == "NAVAL"), expected_naval_count)
   # Manual check: Rows 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14 -> 12 rows
   expect_equal(sum(dict_na_detected$NACats == "NAVAL"), 12)
})

test_that("DetectNAStringsDict uses knownNaLabels and knownNaValues", {
  dict_na_custom <- DetectNAStringsDict(test_dict_na, # Use standard assignment
                                        knownNaLabels = c("Valid Label"), # Mark "Valid Label" as NA
                                        knownNaValues = c("Valid"))       # Mark "Valid" as NA

  # Check rows previously marked as VALUE are now NAVAL
  expect_equal(dict_na_custom$NACats[dict_na_custom$Label == "Valid Label"], rep("NAVAL", 8))
  expect_equal(dict_na_custom$NACats[dict_na_custom$Value == "Valid"], "NAVAL")

  # Check previously NAVAL rows are still NAVAL
  expect_equal(dict_na_custom$NACats[dict_na_custom$Label == "Missing Data"], "NAVAL")
  expect_equal(dict_na_custom$NACats[dict_na_custom$Value == "."], "NAVAL")
})

test_that("DetectNAStringsDict handles no NA strings found", {
  dict_no_na <- data.frame(Variable = "VarA", Label = "Good", Value = "1", Type = "Num", stringsAsFactors = FALSE)
  dict_no_na_detected <- DetectNAStringsDict(dict_no_na) # Use standard assignment
  expect_true("NACats" %in% colnames(dict_no_na_detected))
  expect_equal(dict_no_na_detected$NACats, "VALUE")
})

test_that("DetectNAStringsDict handles empty dictionary", {
  empty_dict <- data.frame(Variable = character(0), Label = character(0), Value = character(0), stringsAsFactors = FALSE)
  detected_empty <- DetectNAStringsDict(empty_dict)
  expect_equal(colnames(detected_empty), c("Variable", "Label", "Value", "NACats"))
  expect_equal(nrow(detected_empty), 0)

  empty_dict_no_cols <- data.frame(Variable = character(0), stringsAsFactors = FALSE)
  expect_error(DetectNAStringsDict(empty_dict_no_cols), "must contain 'Label' and 'Value' columns")
})


# --- Tests for finddatevars ---

test_that("finddatevars identifies date variables by pattern", {
  var_names <- colnames(test_dict_detection)
  detected_dates <- finddatevars(var_names)
  expect_true(all(c("TreatmentDate", "Outcome_DTE", "DT_Start", "End_DT") %in% detected_dates))
  expect_false("Age" %in% detected_dates)
  expect_false("PatientID" %in% detected_dates)
})

test_that("finddatevars includes user-provided knowndatevars", {
  var_names <- colnames(test_dict_detection)
  detected_dates <- finddatevars(var_names, knowndatevars = c("Age", "Score")) # Add non-date vars
  expect_true(all(c("TreatmentDate", "Outcome_DTE", "DT_Start", "End_DT", "Age", "Score") %in% detected_dates))
})

test_that("finddatevars handles case-insensitivity", {
  var_names_case <- c("treatmentdate", "outcome_dte", "age", "DT_START")
  detected_dates <- finddatevars(var_names_case)
  expect_equal(detected_dates, c("treatmentdate", "outcome_dte", "DT_START"))
})

test_that("finddatevars handles no date variables found", {
  var_names_no_dates <- c("PatientID", "Age", "Score", "Status")
  detected_dates <- finddatevars(var_names_no_dates)
  expect_length(detected_dates, 0)
  expect_type(detected_dates, "character")
})

test_that("finddatevars handles empty input vector", {
  empty_names <- character(0)
  expect_length(finddatevars(empty_names), 0)
})
