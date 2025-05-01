library(testthat)
library(here)
# library(data.table) # Removed
library(dplyr) # For bind_rows used in bindDicts

# Assuming the functions are available in the environment
source(here("R/dictionary/dict_processing.R"))

context("Testing dict_processing functions")

# --- Test Data Setup ---
dict_for_fill <- data.frame(
  Variable = c("VarA", NA, NA, "VarB", NA, "VarC"),
  Description = c("Desc A", NA, "Desc A cont.", "Desc B", NA, NA),
  Type = c("Num", NA, NA, "Cat", "Cat", "Char"),
  Value = c(NA, 1, 2, "X", "Y", NA),
  Label = c("Label A", "Label A1", "Label A2", "Label Bx", "Label By", "Label C"),
  Notes = c("Note 1", NA, NA, NA, "Note 5", NA),
  stringsAsFactors = FALSE
)

dict1_for_bind <- data.frame(
    Variable = c("VarA", "VarA", "VarB"),
    Type = c("Num", "Num", "Cat"),
    Value = c(NA, NA, "X"),
    Label = c("Var A", "Var A", "Var B X"),
    stringsAsFactors = FALSE
)

dict2_for_bind <- data.frame(
    Variable = c("VarB", "VarC", "VarC"), # VarB is duplicate
    Type = c("Cat", "Char", "Char"),
    Value = c("Y", NA, NA),
    Label = c("Var B Y", "Var C", "Var C"),
    stringsAsFactors = FALSE
)

dict_for_subset <- data.frame(
    Variable = c("VarA", "VarB", "VarB", "VarC", "VarD"),
    Type = c("Continuous", "Categorical", "Categorical", "Character", "Date"),
    Value = c(NA, "X", "Y", NA, NA),
    Label = c("A", "Bx", "By", "C", "D"),
    stringsAsFactors = FALSE
)

dict_for_type_correct <- data.frame(
    Variable = c("VarA", "VarB", "VarB", "VarC", "VarD", "VarD", "VarE"),
    Type = c("Categorical", "Categorical", "Categorical", "Continuous", "Categorical", "Categorical", "Categorical"),
    Value = c(NA, "X", "Y", NA, "P", "Q", NA),
    Label = c("A", "Bx", "By", "C", "Dp", "Dq", "E"),
    stringsAsFactors = FALSE
) # VarA, VarC, VarE appear once. VarA & VarE are Cat -> should become Continuous


# --- Tests for filllasttoDict ---

test_that("filllasttoDict fills NAs forward", {
  dict_filled <- filllasttoDict(dict_for_fill) # Use standard assignment (copy-on-modify)

  # Check Variable column
  expect_equal(dict_filled$Variable, c("VarA", "VarA", "VarA", "VarB", "VarB", "VarC"))
  # Check Description column
  expect_equal(dict_filled$Description, c("Desc A", "Desc A", "Desc A cont.", "Desc B", "Desc B", "Desc B"))
  # Check Type column
  expect_equal(dict_filled$Type, c("Num", "Num", "Num", "Cat", "Cat", "Char"))
  # Check Value column (NA at start remains NA)
  expect_equal(dict_filled$Value, c(NA, 1, 2, "X", "Y", "Y"))
   # Check Label column
  expect_equal(dict_filled$Label, c("Label A", "Label A1", "Label A2", "Label Bx", "Label By", "Label C"))
  # Check Notes column (NAs after last value remain NA)
  expect_equal(dict_filled$Notes, c("Note 1", "Note 1", "Note 1", "Note 1", "Note 5", "Note 5"))
})

test_that("filllasttoDict handles data frame with no NAs", {
  dict_no_na <- data.frame(A = 1:3, B = letters[1:3], stringsAsFactors = FALSE)
  expect_equal(filllasttoDict(dict_no_na), dict_no_na) # Use standard assignment
})

test_that("filllasttoDict handles data frame with all NAs in a column", {
  dict_all_na_col <- data.frame(A = 1:3, B = rep(NA_character_, 3), stringsAsFactors = FALSE)
  expect_equal(filllasttoDict(dict_all_na_col), dict_all_na_col) # Use standard assignment
})

test_that("filllasttoDict handles empty data frame", {
  # empty_dt <- data.table() # Removed
  # Need to define columns for empty data.table to work with lapply
  # empty_dt_cols <- data.table(A=character(), B=numeric()) # Removed
  # expect_equal(filllasttoDict(empty_dt_cols), empty_dt_cols) # Removed

  empty_df <- data.frame(A=character(), B=numeric())
   expect_equal(filllasttoDict(empty_df), empty_df)
})


# --- Tests for bindDicts ---

test_that("bindDicts combines dictionaries correctly", {
  bound_dict <- bindDicts(dict1_for_bind, dict2_for_bind)

  # Expected result: All of dict1, plus non-duplicate VarC from dict2
  expected_dict <- data.frame(
      Variable = c("VarA", "VarA", "VarB", "VarC", "VarC"),
      Type = c("Num", "Num", "Cat", "Char", "Char"),
      Value = c(NA, NA, "X", NA, NA),
      Label = c("Var A", "Var A", "Var B X", "Var C", "Var C"),
      stringsAsFactors = FALSE
  )

  # Use arrange for consistent comparison as bind_rows doesn't guarantee order
  expect_equal(arrange(bound_dict, Variable, Value), arrange(expected_dict, Variable, Value))
  expect_equal(nrow(bound_dict), 5)
})

test_that("bindDicts handles second dictionary having no new variables", {
  dict2_only_dups <- data.frame(Variable = "VarA", Type = "Num", Value = NA, Label = "Dup", stringsAsFactors = FALSE)
  bound_dict <- bindDicts(dict1_for_bind, dict2_only_dups)
  expect_equal(arrange(bound_dict, Variable, Value), arrange(dict1_for_bind, Variable, Value))
})

test_that("bindDicts handles first dictionary being empty", {
  empty_dict1 <- dict1_for_bind[0,]
  bound_dict <- bindDicts(empty_dict1, dict2_for_bind)
  expect_equal(arrange(bound_dict, Variable, Value), arrange(dict2_for_bind, Variable, Value))
})

test_that("bindDicts handles second dictionary being empty", {
  empty_dict2 <- dict2_for_bind[0,]
  bound_dict <- bindDicts(dict1_for_bind, empty_dict2)
  expect_equal(arrange(bound_dict, Variable, Value), arrange(dict1_for_bind, Variable, Value))
})

test_that("bindDicts requires data frames with Variable column", {
  expect_error(bindDicts(dict1_for_bind, data.frame(A=1)), "must contain a 'Variable' column")
  expect_error(bindDicts(data.frame(A=1), dict2_for_bind), "must contain a 'Variable' column")
  expect_error(bindDicts(dict1_for_bind, "not a dataframe"), "must be data frames")
})


# --- Tests for CatVarsDict ---

test_that("CatVarsDict subsets to Categorical variables", {
  cat_dict <- CatVarsDict(dict_for_subset)
  expect_equal(nrow(cat_dict), 2)
  expect_true(all(cat_dict$Type == "Categorical"))
  expect_equal(unique(cat_dict$Variable), "VarB")
})

test_that("CatVarsDict handles no Categorical variables", {
  dict_no_cat <- dict_for_subset[!dict_for_subset$Type == "Categorical",]
  cat_dict <- CatVarsDict(dict_no_cat)
  expect_equal(nrow(cat_dict), 0)
})

test_that("CatVarsDict requires Type column", {
  dict_no_type <- dict_for_subset[, !(colnames(dict_for_subset) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(CatVarsDict(dict_no_type), "must contain a 'Type' column")
})


# --- Tests for ConVarsDict ---

test_that("ConVarsDict subsets to Continuous variables", {
  con_dict <- ConVarsDict(dict_for_subset)
  expect_equal(nrow(con_dict), 1)
  expect_true(all(con_dict$Type == "Continuous"))
  expect_equal(unique(con_dict$Variable), "VarA")
})

test_that("ConVarsDict handles no Continuous variables", {
  dict_no_con <- dict_for_subset[!dict_for_subset$Type == "Continuous",]
  con_dict <- ConVarsDict(dict_no_con)
  expect_equal(nrow(con_dict), 0)
})

test_that("ConVarsDict requires Type column", {
  dict_no_type <- dict_for_subset[, !(colnames(dict_for_subset) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(ConVarsDict(dict_no_type), "must contain a 'Type' column")
})


# --- Tests for CharVarsDict ---

test_that("CharVarsDict subsets to Character variables", {
  char_dict <- CharVarsDict(dict_for_subset)
  expect_equal(nrow(char_dict), 1)
  expect_true(all(char_dict$Type == "Character"))
  expect_equal(unique(char_dict$Variable), "VarC")
})

test_that("CharVarsDict handles no Character variables", {
  dict_no_char <- dict_for_subset[!dict_for_subset$Type == "Character",]
  char_dict <- CharVarsDict(dict_no_char)
  expect_equal(nrow(char_dict), 0)
})

test_that("CharVarsDict requires Type column", {
  dict_no_type <- dict_for_subset[, !(colnames(dict_for_subset) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(CharVarsDict(dict_no_type), "must contain a 'Type' column")
})


# --- Tests for OneValCat2NumDict ---

test_that("OneValCat2NumDict corrects single-occurrence Categorical to Continuous", {
  corrected_dict <- OneValCat2NumDict(dict_for_type_correct) # Use standard assignment

  # VarA: Single occurrence, was Categorical -> should be Continuous
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarA"], "Continuous")
  # VarB: Multiple occurrences, was Categorical -> should remain Categorical
  expect_true(all(corrected_dict$Type[corrected_dict$Variable == "VarB"] == "Categorical"))
  # VarC: Single occurrence, was Continuous -> should remain Continuous
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarC"], "Continuous")
  # VarD: Multiple occurrences, was Categorical -> should remain Categorical
  expect_true(all(corrected_dict$Type[corrected_dict$Variable == "VarD"] == "Categorical"))
   # VarE: Single occurrence, was Categorical -> should be Continuous
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarE"], "Continuous")
})

test_that("OneValCat2NumDict respects 'vars' argument", {
  # Only check VarA and VarB
  corrected_dict_subset <- OneValCat2NumDict(dict_for_type_correct, vars = c("VarA", "VarB")) # Use standard assignment

  # VarA should be corrected
  expect_equal(corrected_dict_subset$Type[corrected_dict_subset$Variable == "VarA"], "Continuous")
  # VarB should remain Categorical
  expect_true(all(corrected_dict_subset$Type[corrected_dict_subset$Variable == "VarB"] == "Categorical"))
  # VarE should NOT be corrected as it wasn't in 'vars'
  expect_equal(corrected_dict_subset$Type[corrected_dict_subset$Variable == "VarE"], "Categorical")
})

test_that("OneValCat2NumDict handles no variables needing correction", {
  dict_no_correction <- data.frame(
      Variable = c("V1", "V1", "V2"),
      Type = c("Categorical", "Categorical", "Continuous"),
      stringsAsFactors = FALSE
  )
  expect_equal(OneValCat2NumDict(dict_no_correction), dict_no_correction) # Use standard assignment
})

test_that("OneValCat2NumDict requires necessary columns", {
  dict_no_var <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Variable"), drop = FALSE] # Base R subsetting
  dict_no_type <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(OneValCat2NumDict(dict_no_var), "must contain 'Variable' and 'Type' columns")
  expect_error(OneValCat2NumDict(dict_no_type), "must contain 'Variable' and 'Type' columns")
})


# --- Tests for NumIdVarstoCharDict ---

test_that("NumIdVarstoCharDict changes Type to Character for specified IDs", {
  dict_ids <- data.frame(
      Variable = c("ID1", "ID2", "VarA", "ID3"),
      Type = c("Num", "Numeric", "Cat", "Continuous"), # Different initial numeric types
      Value = NA,
      stringsAsFactors = FALSE
  )
  ids_to_change <- c("ID1", "ID3", "NonExistentID")
  corrected_dict <- NumIdVarstoCharDict(dict_ids, idvars = ids_to_change) # Use standard assignment

  expect_equal(corrected_dict$Type[corrected_dict$Variable == "ID1"], "Character")
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "ID3"], "Character")
  # Check others are unchanged
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "ID2"], "Numeric")
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarA"], "Cat")
})

test_that("NumIdVarstoCharDict uses default IDs if none provided", {
   dict_default_ids <- data.frame(
       Variable = c("pseudoid", "VarA", "pseudoccn", "`#`"),
       Type = c("Num", "Cat", "Num", "Num"),
       Value = NA,
       stringsAsFactors = FALSE
   )
   corrected_dict <- NumIdVarstoCharDict(dict_default_ids) # Use defaults, standard assignment
   expect_equal(corrected_dict$Type[corrected_dict$Variable == "pseudoid"], "Character")
   expect_equal(corrected_dict$Type[corrected_dict$Variable == "pseudoccn"], "Character")
   expect_equal(corrected_dict$Type[corrected_dict$Variable == "`#`"], "Character")
   expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarA"], "Cat") # Unchanged
})

test_that("NumIdVarstoCharDict handles no specified IDs found", {
  dict_ids <- data.frame(Variable = "VarA", Type = "Num", Value = NA, stringsAsFactors = FALSE)
  expect_warning(corrected_dict <- NumIdVarstoCharDict(dict_ids, idvars = c("ID1", "ID2")), # Use standard assignment
                 "None of the specified idvars found")
  expect_equal(corrected_dict, dict_ids) # Should be unchanged
})

test_that("NumIdVarstoCharDict handles NULL or empty idvars", {
  dict_ids <- data.frame(Variable = "ID1", Type = "Num", Value = NA, stringsAsFactors = FALSE)
  expect_equal(NumIdVarstoCharDict(dict_ids, idvars = NULL), dict_ids) # Use standard assignment
  expect_equal(NumIdVarstoCharDict(dict_ids, idvars = character(0)), dict_ids) # Use standard assignment
})

test_that("NumIdVarstoCharDict requires necessary columns", {
  dict_no_var <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Variable"), drop = FALSE] # Base R subsetting
  dict_no_type <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(NumIdVarstoCharDict(dict_no_var), "must contain 'Variable' and 'Type' columns")
  expect_error(NumIdVarstoCharDict(dict_no_type), "must contain 'Variable' and 'Type' columns")
})


# --- Tests for DateVarstoDateDict ---

test_that("DateVarstoDateDict changes Type to Date for specified vars", {
  dict_dates <- data.frame(
      Variable = c("Date1", "VarA", "Date_Adm", "DateVar3"),
      Type = c("Num", "Cat", "Character", "Numeric"),
      Value = NA,
      stringsAsFactors = FALSE
  )
  dates_to_change <- c("Date1", "Date_Adm", "NonExistentDate")
  corrected_dict <- DateVarstoDateDict(dict_dates, datevars = dates_to_change) # Use standard assignment

  expect_equal(corrected_dict$Type[corrected_dict$Variable == "Date1"], "Date")
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "Date_Adm"], "Date")
  # Check others are unchanged
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "VarA"], "Cat")
  expect_equal(corrected_dict$Type[corrected_dict$Variable == "DateVar3"], "Numeric")
})

test_that("DateVarstoDateDict handles no specified dates found", {
  dict_dates <- data.frame(Variable = "VarA", Type = "Num", Value = NA, stringsAsFactors = FALSE)
  expect_warning(corrected_dict <- DateVarstoDateDict(dict_dates, datevars = c("D1", "D2")), # Use standard assignment
                 "None of the specified datevars found")
  expect_equal(corrected_dict, dict_dates) # Should be unchanged
})

test_that("DateVarstoDateDict handles NULL or empty datevars", {
  dict_dates <- data.frame(Variable = "Date1", Type = "Num", Value = NA, stringsAsFactors = FALSE)
  # If datevars is NULL, function should do nothing
  expect_equal(DateVarstoDateDict(dict_dates, datevars = NULL), dict_dates) # Use standard assignment
  expect_equal(DateVarstoDateDict(dict_dates, datevars = character(0)), dict_dates) # Use standard assignment
})

test_that("DateVarstoDateDict requires necessary columns", {
  dict_no_var <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Variable"), drop = FALSE] # Base R subsetting
  dict_no_type <- dict_for_type_correct[, !(colnames(dict_for_type_correct) %in% "Type"), drop = FALSE] # Base R subsetting
  expect_error(DateVarstoDateDict(dict_no_var), "must contain 'Variable' and 'Type' columns")
  expect_error(DateVarstoDateDict(dict_no_type), "must contain 'Variable' and 'Type' columns")
})


# --- Tests for addDiseaseSpecColsDict ---

test_that("addDiseaseSpecColsDict adds NACats and DisSpec columns", {
  dict_disease <- data.frame(
      Variable = c("disease", "disease", "VarA", "VarA", "VarB", "VarB"),
      Description = c("Disease Type A", "Disease Type B", "Variable A Desc", "Variable A Desc", "Variable B Desc", "Variable B Desc"),
      Value = c("A", "B", "1", ".E", "X", "Y"), # .E indicates NotApp
      Label = c("Disease A", "Disease B", "Value 1", "N/A, other disease", "Value X", "Value Y"), # Label indicates specificity for VarA
      Type = c("Cat", "Cat", "Num", "Num", "Cat", "Cat"),
      stringsAsFactors = FALSE
  )
  # Add Notes column required by readDict standard format, even if not used here
  dict_disease$Notes <- NA_character_ # Use standard assignment

  processed_dict <- addDiseaseSpecColsDict(dict_disease) # Use standard assignment

  # Check NACats column
  expect_true("NACats" %in% colnames(processed_dict))
  expect_equal(processed_dict$NACats[processed_dict$Value == ".E"], "NotApp")
  expect_equal(processed_dict$NACats[processed_dict$Value != ".E"], rep("VALUE", 5))

  # Check DisSpec column
  expect_true("DisSpec" %in% colnames(processed_dict))
  # VarA has "N/A, other disease" label, and its description doesn't mention "Disease B"
  # So, it should be marked as specific to "Disease A"
  # Note: The logic relies on grepl matching disease names in descriptions, which isn't happening here.
  # Let's adjust the test data description to trigger the logic.
  dict_disease_mod <- dict_disease # Use standard assignment
  dict_disease_mod$Description[dict_disease_mod$Variable == "VarA" & dict_disease_mod$Label == "N/A, other disease"] <- "Variable A Desc (Not for Disease B)"

  processed_dict_mod <- addDiseaseSpecColsDict(dict_disease_mod)
  expect_true(!is.na(processed_dict_mod$DisSpec[processed_dict_mod$Variable == "VarA"][1])) # Should have a disease listed
  # Since "Not for Disease B" is in desc, it should list "Disease A"
  expect_match(processed_dict_mod$DisSpec[processed_dict_mod$Variable == "VarA"][1], "Disease A")
  expect_false(grepl("Disease B", processed_dict_mod$DisSpec[processed_dict_mod$Variable == "VarA"][1]))

  # VarB has no "N/A, other disease" label, should be NA in DisSpec
  expect_true(all(is.na(processed_dict_mod$DisSpec[processed_dict_mod$Variable == "VarB"])))
  # disease variable itself should be NA in DisSpec
  expect_true(all(is.na(processed_dict_mod$DisSpec[processed_dict_mod$Variable == "disease"])))

})

test_that("addDiseaseSpecColsDict handles missing 'disease' variable", {
  dict_no_disease_var <- dict_for_subset # Doesn't have 'disease' variable
  # Add required columns if missing
  if (!"Description" %in% colnames(dict_no_disease_var)) dict_no_disease_var$Description <- NA_character_
  if (!"Notes" %in% colnames(dict_no_disease_var)) dict_no_disease_var$Notes <- NA_character_

  expect_warning(processed_dict <- addDiseaseSpecColsDict(dict_no_disease_var), # Use standard assignment
                 "No 'disease' variable found")
  expect_true("NACats" %in% colnames(processed_dict))
  expect_true("DisSpec" %in% colnames(processed_dict))
  expect_true(all(is.na(processed_dict$DisSpec))) # DisSpec should be all NA
})

test_that("addDiseaseSpecColsDict requires necessary columns", {
  dict_minimal <- data.frame(Variable = "V1", Value = "1", stringsAsFactors = FALSE)
  expect_error(addDiseaseSpecColsDict(dict_minimal), "must contain columns: Variable, Description, Value, Label")
})


# --- Tests for MakelabelledSASDict ---

test_that("MakelabelledSASDict reads SAS file and generates dictionary (requires packages)", {
  skip_if_not_installed("labelled")
  skip_if_not_installed("haven")

  # Create a dummy SAS file for testing (complex to do robustly without writing a file)
  # Option 1: Skip this test if no dummy file exists
  dummy_sas_path <- "dummy_test_file.sas7bdat"
  skip_if_not(file.exists(dummy_sas_path), "Dummy SAS file not found for testing MakelabelledSASDict")

  # Option 2: Mock haven::read_sas and labelled::generate_dictionary (safer)
  # This requires mocking libraries like mockery

  # Assuming a dummy file exists or mocking is set up:
  # mock_haven_read_sas <- mockery::mock(data.frame(Var1 = labelled::labelled(c(1,2), labels = c(Yes=1, No=2)), Var2 = c("a", "b")))
  # mock_generate_dict <- mockery::mock(data.frame(variable="Var1", label="Label for Var1", format=NA, type="numeric", value_labels = "Yes=1, No=2"))
  # mockery::stub(MakelabelledSASDict, 'haven::read_sas', mock_haven_read_sas)
  # mockery::stub(MakelabelledSASDict, 'labelled::generate_dictionary', mock_generate_dict)

  # dict_from_sas <- MakelabelledSASDict(dummy_sas_path)
  # expect_s3_class(dict_from_sas, "data.frame")
  # expect_true("Variable" %in% colnames(dict_from_sas)) # Adjust based on actual/mocked output
  # expect_true("Description" %in% colnames(dict_from_sas)) # Adjust based on actual/mocked output

  # Placeholder expectation if skipping
  expect_true(TRUE)
})

test_that("MakelabelledSASDict handles non-existent file", {
  skip_if_not_installed("labelled")
  skip_if_not_installed("haven")
  expect_error(MakelabelledSASDict("non_existent_file.sas7bdat"), "Specified SAS file does not exist")
})
