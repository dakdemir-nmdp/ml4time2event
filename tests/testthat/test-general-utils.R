library(testthat)
library(here)
# library(data.table) # Removed
library(partykit) # For listrules tests

# Assuming the functions are available in the environment
# Use here::here for robustness
source(here::here("R/general_utils.R"))

context("Testing general_utils functions")

# --- Test Data Setup ---
test_data_profile <- data.frame(
  num_var = c(1, 2, 3, 4, 5, NA),
  fact_var = factor(c("A", "B", "A", "C", "B", NA), levels = c("A", "B", "C")),
  char_var = c("X", "Y", "X", "Z", "Y", "X"),
  logi_var = c(TRUE, FALSE, TRUE, TRUE, FALSE, NA), # Unsupported type for explicit handling
  const_num = rep(5, 6),
  const_fact = factor(rep("K", 6)),
  stringsAsFactors = FALSE
)
expvars_profile <- c("num_var", "fact_var", "char_var", "logi_var", "const_num", "const_fact", "non_existent_var")


# --- Tests for VariableProfile ---

test_that("VariableProfile returns a list named by expvars", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  expect_type(profile, "list")
  expect_named(profile, expvars_profile)
  expect_length(profile, length(expvars_profile))
})

test_that("VariableProfile profiles numeric variables correctly", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  # num_var
  expect_equal(profile$num_var, c(min = 1, max = 5))
  # const_num
  expect_equal(profile$const_num, c(min = 5, max = 5))
})

test_that("VariableProfile profiles factor variables correctly (with NAs)", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  # fact_var
  expected_fact_table <- table(test_data_profile$fact_var, useNA = "ifany")
  expect_equivalent(profile$fact_var, expected_fact_table)
  # Check NA level is present by checking if NA is among the names
  expect_true(any(is.na(names(profile$fact_var))))


  # const_fact
  expected_const_fact_table <- table(test_data_profile$const_fact, useNA = "ifany")
  expect_equivalent(profile$const_fact, expected_const_fact_table)
})

test_that("VariableProfile profiles character variables correctly (with NAs)", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  # char_var (no NAs in this test data)
  expected_char_table <- table(test_data_profile$char_var, useNA = "ifany")
  expect_equivalent(profile$char_var, expected_char_table)
})

test_that("VariableProfile handles unsupported types", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  expect_equal(profile$logi_var, "Unsupported type: logical")
})

test_that("VariableProfile handles variables not found in data", {
  profile <- VariableProfile(test_data_profile, expvars_profile)
  expect_equal(profile$non_existent_var, "Variable not found in data")
})

test_that("VariableProfile handles empty data or expvars", {
  empty_data <- data.frame() # Replaced data.table()
  expect_equal(VariableProfile(empty_data, "A"), list(A = "Variable not found in data"))

  # Check that providing an empty character vector results in an empty list-like object
  profile_empty_vars <- VariableProfile(test_data_profile, character(0))
  expect_type(profile_empty_vars, "list")
  expect_length(profile_empty_vars, 0)
})


# --- Tests for listrules ---

# Setup data and tree for listrules tests
listrules_data <- data.frame(
    y = factor(c(rep("A", 5), rep("B", 5), rep("C", 5))),
    x1 = c(1:10, 20:24),
    x2 = factor(rep(c("G1", "G2", "G3"), each = 5)),
    stringsAsFactors = FALSE
)

# Build a simple tree using ctree
# Ensure partykit is installed
if (requireNamespace("partykit", quietly = TRUE)) {
    listrules_tree <- partykit::ctree(y ~ x1 + x2, data = listrules_data,
                                      control = partykit::ctree_control(minsplit = 4, minbucket = 2))

    test_that("listrules extracts rules for terminal nodes", {
      skip_if_not_installed("partykit")
      rules <- listrules(listrules_tree)

      # Check output type and names - sapply might simplify to vector
      expect_type(rules, "character")
      terminal_nodes <- partykit::nodeids(listrules_tree, terminal = TRUE)
      expect_named(rules, as.character(terminal_nodes))

      # Check that rules are non-empty strings (basic check)
      expect_true(all(sapply(rules, function(r) is.character(r) && nchar(r) > 0)))
    })

    test_that("listrules extracts rules for specific node ID", {
      skip_if_not_installed("partykit")
      rule_node5 <- listrules(listrules_tree, i = 5) # Get rule for node 5
      expect_type(rule_node5, "character")
      expect_length(rule_node5, 1)
      expect_true(nchar(rule_node5) > 0) # Basic check
    })

    test_that("listrules handles non-terminal nodes (returns empty string)", {
       skip_if_not_installed("partykit")
       # Rule for root node (ID 1) should be empty
       expect_equal(listrules(listrules_tree, i = 1), "")
       # Rule for internal node 2 - Should involve x2 split for G1/G2 vs G3
       rule_node2 <- listrules(listrules_tree, i = 2)
       expect_true(nchar(rule_node2) > 0) # Basic check
    })

} else {
    message("Skipping listrules tests because partykit is not installed.")
}
