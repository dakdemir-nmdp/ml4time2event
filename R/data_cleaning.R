
#' @title RemoveMonoVarsData
#'
#' @description Remove variables (columns) from a data frame that have only one unique non-NA value.
#'
#' @param data Data frame.
#' @return Data frame with single-value columns removed.
#' @export
RemoveMonoVarsData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")

  # Handle empty data frame
  if (ncol(data) == 0 || nrow(data) == 0) {
    return(data)
  }

  num_unique_vals <- sapply(data, function(x) length(unique(x[!is.na(x)])))
  # Ensure it's a numeric vector (sapply can return a list for empty data)
  num_unique_vals <- as.numeric(num_unique_vals)

  # Only remove columns with exactly 1 unique non-NA value (true monotonous)
  # Keep all-NA columns (0 unique values) for RemoveAllNAVars to handle
  cols_to_keep <- num_unique_vals != 1
  if (sum(!cols_to_keep) > 0) {
      cat("Removing", sum(!cols_to_keep), "columns with only one unique non-NA value:", paste(colnames(data)[!cols_to_keep], collapse=", "), "\n")
  }
  return(data[, cols_to_keep, drop = FALSE])
}

#' @title RemoveAllNAVars
#'
#' @description Remove variables (columns) from a data frame that consist entirely of NA values.
#'
#' @param data Data frame.
#' @return Data frame with all-NA columns removed.
#' @export
RemoveAllNAVars <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Handle empty data frame
  if (ncol(data) == 0 || nrow(data) == 0) {
    return(data)
  }
  
  all_na_cols <- sapply(data, function(x) all(is.na(x)))
  # Ensure it's a logical vector (sapply can return a list for empty data)
  all_na_cols <- as.logical(all_na_cols)
  
  if (sum(all_na_cols) > 0) {
      cat("Removing", sum(all_na_cols), "columns containing only NA values:", paste(colnames(data)[all_na_cols], collapse=", "), "\n")
  }
  return(data[, !all_na_cols, drop = FALSE])
}



#' RemoveMissinVarsData
#' @title RemoveMissinVarsData
#'
#' @description Remove variables (columns) with a proportion of missing values exceeding a threshold.
#'
#' @param data Data frame.
#' @param maxprop Numeric threshold (0-1) for the maximum allowed proportion of NA values (default: 0.2).
#' @return Data frame with high-missingness columns removed.
#' @export
RemoveMissinVarsData <- function(data, maxprop = .2) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!is.numeric(maxprop) || maxprop < 0 || maxprop > 1) stop("'maxprop' must be between 0 and 1.")

  na_props <- sapply(data, function(x) sum(is.na(x)) / length(x))
  cols_to_keep <- na_props <= maxprop
  cols_removed <- colnames(data)[!cols_to_keep]

  if (length(cols_removed) > 0) {
      cat("Removing", length(cols_removed), "columns with >=", maxprop*100, "% missing values:", paste(cols_removed, collapse=", "), "\n")
  }
  return(data[, cols_to_keep, drop = FALSE])
}

#' @title RemoveMissinRecordsData
#'
#' @description Remove records (rows) with a proportion of missing values exceeding a threshold.
#'
#' @param data Data frame.
#' @param maxprop Numeric threshold (0-1) for the maximum allowed proportion of NA values per row (default: 0.2).
#' @return Data frame with high-missingness rows removed.
#' @export
RemoveMissinRecordsData <- function(data, maxprop = .2) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!is.numeric(maxprop) || maxprop < 0 || maxprop > 1) stop("'maxprop' must be between 0 and 1.")

  na_props_row <- apply(data, 1, function(x) sum(is.na(x)) / length(x))
  rows_to_keep <- na_props_row <= maxprop
  rows_removed_count <- sum(!rows_to_keep)

  if (rows_removed_count > 0) {
      cat("Removing", rows_removed_count, "rows with >=", maxprop*100, "% missing values.\n")
  }
  return(data[rows_to_keep, , drop = FALSE])
}



#' @title getcharcolsData
#'
#' @description Utility function to get the names of character columns in a data frame.
#'
#' @param data Data frame.
#' @return Character vector with names of character columns.
#' @export
getcharcolsData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Handle empty data frame
  if (ncol(data) == 0) {
    return(character(0))
  }
  
  charcols_logical <- sapply(data, is.character)
  # Ensure it's a logical vector (sapply can return a list for empty data)
  charcols_logical <- as.logical(charcols_logical)
  
  return(colnames(data)[charcols_logical])
}


#' @title ImputeMissinRecordsData
#'
#' @description Impute missing values in a data frame using missRanger.
#' Excludes specified columns (e.g., character IDs, outcome variables) from the imputation process.
#'
#' @param data Data frame to be imputed.
#' @param dontuse Character vector of column names to exclude from imputation (neither used for prediction nor imputed).
#' @param ... Additional arguments passed to `missRanger::missRanger`.
#' @return Data frame with missing values imputed (excluding 'dontuse' columns).
#' @seealso [missRanger::missRanger()]
#' @examples
#' \dontrun{
#' if (requireNamespace("missRanger", quietly = TRUE)) {
#'   iris_na <- iris
#'   iris_na[1:3, 1] <- NA
#'   iris_imputed <- ImputeMissinRecordsData(iris_na, dontuse = "Species")
#' }
#' }
#' @importFrom missRanger missRanger
#' @export
ImputeMissinRecordsData <- function(data, dontuse = NULL, ...) {
  if (!requireNamespace("missRanger", quietly = TRUE)) {
    stop("Package 'missRanger' needed for this function. Please install it.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("'data' must be a data frame.")

  original_colnames <- colnames(data)
  cols_to_impute <- setdiff(original_colnames, dontuse)

  if (length(cols_to_impute) == 0) {
      warning("No columns selected for imputation after excluding 'dontuse' columns.")
      return(data)
  }

  data_to_impute <- data[, cols_to_impute, drop = FALSE]

  # Check for columns with zero variance or all NA after subsetting
  all_na <- sapply(data_to_impute, function(x) all(is.na(x)))
  if(any(all_na)) {
      warning("Columns contain only NA values and cannot be imputed: ", paste(colnames(data_to_impute)[all_na], collapse=", "))
      # Optionally remove them before imputation
      data_to_impute <- data_to_impute[, !all_na, drop = FALSE]
      if (ncol(data_to_impute) == 0) {
          warning("No columns left to impute after removing all-NA columns.")
          return(data)
      }
  }
  # missRanger might handle zero variance, but checking can be useful
  # zero_var <- sapply(data_to_impute, function(x) length(unique(x[!is.na(x)])) <= 1)
  # if(any(zero_var)) {
  #     warning("Columns contain zero variance (excluding NAs): ", paste(colnames(data_to_impute)[zero_var], collapse=", "))
  # }


  # Make column names syntactically valid if needed (missRanger might require this)
  impute_colnames_orig <- colnames(data_to_impute)
  impute_colnames_valid <- make.names(impute_colnames_orig, unique = TRUE)
  colnames(data_to_impute) <- impute_colnames_valid

  cat("Starting imputation with missRanger for", ncol(data_to_impute), "columns...\n")
  data_imputed_valid_names <- tryCatch(
      missRanger::missRanger(data_to_impute, verbose = 1, ...), # Set verbose=1 for basic progress
      error = function(e) {
          stop("Error during missRanger imputation: ", e$message)
      }
  )
  cat("Imputation finished.\n")

  # Restore original column names
  colnames(data_imputed_valid_names) <- impute_colnames_orig

  # Combine imputed data with excluded columns
  data_out <- data
  data_out[, impute_colnames_orig] <- data_imputed_valid_names

  return(data_out)
}




#' @title RemoveRareCategoriesData
#'
#' @description Replace rare categories in factor variables with NA.
#'
#' @param data Data frame.
#' @param minfreq Numeric threshold (0-1). Categories with a frequency below this proportion will be set to NA (default: 0.01).
#' @return Data frame with rare factor levels replaced by NA and unused levels dropped.
#' @export
RemoveRareCategoriesData <- function(data, minfreq = .01) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!is.numeric(minfreq) || minfreq < 0 || minfreq > 1) stop("'minfreq' must be between 0 and 1.")

  data_out <- data
  for (vari in colnames(data_out)) {
    col_data <- data_out[[vari]]
    if (is.factor(col_data)) {
      counts <- table(col_data)
      props <- prop.table(counts)
      rare_levels <- names(props)[props < minfreq]

      if (length(rare_levels) > 0) {
          cat("Variable '", vari, "': Replacing rare levels (freq < ", minfreq, ") with NA: ", paste(rare_levels, collapse=", "), "\n")
          # Convert to character, replace, then back to factor (to handle level changes)
          col_data_char <- as.character(col_data)
          col_data_char[col_data_char %in% rare_levels] <- NA
          # Keep original levels minus the rare ones
          remaining_levels <- levels(col_data)[!levels(col_data) %in% rare_levels]
          data_out[[vari]] <- factor(col_data_char, levels = remaining_levels)
      }
    }
  }
  return(data_out)
}


#' @title RemoveRareBinaryVarsData
#'
#' @description For binary variables (numeric 0/1, logical, or 2-level factors), replace the rarer category with NA if its frequency is below a threshold.
#'
#' @param data Data frame.
#' @param minfreq Numeric threshold (0-1). The rarer category is set to NA if its frequency is below this (default: 0.01).
#' @return Data frame with rare binary categories replaced by NA.
#' @export
RemoveRareBinaryVarsData <- function(data, minfreq = .01) {
   if (!is.data.frame(data)) stop("'data' must be a data frame.")
   if (!is.numeric(minfreq) || minfreq < 0 || minfreq > 1) stop("'minfreq' must be between 0 and 1.")

   data_out <- data
   for (vari in colnames(data_out)) {
       col_data <- data_out[[vari]]
       unique_non_na <- unique(col_data[!is.na(col_data)])

       # Only process binary variables: numeric, logical, or factors (not character)
       if (length(unique_non_na) == 2 && !is.character(col_data)) { # Check if binary (ignoring NAs)
           counts <- table(col_data, useNA = "no") # Exclude NA from table
           props <- prop.table(counts)
           min_prop <- min(props)
           rare_level <- names(props)[which.min(props)]

           if (min_prop < minfreq) {
               cat("Variable '", vari, "': Replacing rare level '", rare_level, "' (freq = ", round(min_prop, 4), ") with NA.\n")
               # Need to handle different types (factor, numeric, logical)
               if (is.factor(col_data)) {
                   col_data_char <- as.character(col_data)
                   col_data_char[col_data_char == rare_level] <- NA
                   # Keep only the non-rare level
                   remaining_level <- names(props)[which.max(props)]
                   data_out[[vari]] <- factor(col_data_char, levels = remaining_level)
               } else { # Numeric or Logical
                   # Convert rare_level back to original type for comparison
                   original_type <- class(unique_non_na)
                   rare_level_typed <- methods::as(rare_level, original_type)
                   col_data[col_data == rare_level_typed & !is.na(col_data)] <- NA
                   data_out[[vari]] <- col_data
               }
           }
       }
   }
   return(data_out)
}


#' @title CollapseRareCategoriestoOtherData
#'
#' @description Collapse rare categories in factor or character variables into a new level called "Other".
#' Only applies to variables with more than 3 unique levels initially.
#'
#' @param data Data frame.
#' @param minfreq Numeric threshold (0-1). Categories with frequency below this proportion are collapsed (default: 0.01).
#' @return Data frame with rare categories collapsed into "Other".
#' @export
CollapseRareCategoriestoOtherData <- function(data, minfreq = .01) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!is.numeric(minfreq) || minfreq < 0 || minfreq > 1) stop("'minfreq' must be between 0 and 1.")

  data_out <- data
  for (vari in colnames(data_out)) {
    col_data <- data_out[[vari]]
    if (is.factor(col_data) || is.character(col_data)) {
      col_data_char <- as.character(col_data)
      counts <- table(col_data_char, useNA = "no") # Exclude NA

      if (length(counts) > 3) { # Only collapse if more than 3 levels
          props <- prop.table(counts)
          levels_to_collapse <- names(props)[props < minfreq]

          if (length(levels_to_collapse) > 0) {
              cat("Variable '", vari, "': Collapsing levels (freq < ", minfreq, ") to 'Other': ", paste(levels_to_collapse, collapse=", "), "\n")
              col_data_char[col_data_char %in% levels_to_collapse] <- "Other"
              # Convert back to factor, ensuring "Other" is a level
              data_out[[vari]] <- factor(col_data_char)
          }
      }
    }
  }
  return(data_out)
}



#' @title droplevelsoffactorsData
#'
#' @description Drop unused factor levels from all factor columns in a data frame.
#'
#' @param data Data frame.
#' @return Data frame with unused factor levels dropped.
#' @export
droplevelsoffactorsData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  data_out <- data
  factor_cols <- sapply(data_out, is.factor)
  if (any(factor_cols)) {
      data_out[factor_cols] <- lapply(data_out[factor_cols], droplevels)
      cat("Dropped unused levels for factor columns:", paste(colnames(data_out)[factor_cols], collapse=", "), "\n")
  }
  return(data_out)
}

#' @title findvarsnamesthatrepeatData
#'
#' @description Find variable names where one name is a substring of another (potential redundancy).
#' This is a simple substring check, not checking for semantic similarity.
#'
#' @param data Data frame.
#' @return A named character vector where names are the shorter variable names and
#'   values are comma-separated lists of longer variable names containing the shorter name.
#' @export
findvarsnamesthatrepeatData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  col_names <- colnames(data)
  
  # Handle empty data frame
  if (length(col_names) == 0) {
    cat("No simple substring relationships found between variable names.\n")
    return(character(0))
  }
  
  out_list <- list()

  for (i in seq_along(col_names)) {
    vari <- col_names[i]
    # Find other column names that contain vari as a substring
    # Ensure it's not just matching itself and use word boundaries for better matching?
    # Simple substring match as in original:
    matches <- grep(vari, col_names[-i], fixed = TRUE, value = TRUE) # Use fixed=TRUE for literal match

    # Alternative: word boundary match (might be too strict)
    # pattern <- paste0("\\b", vari, "\\b")
    # matches <- grep(pattern, col_names[-i], value = TRUE)

    if (length(matches) > 0) {
      out_list[[vari]] <- paste(matches, collapse = ", ")
    }
  }

  if (length(out_list) > 0) {
      cat("Found potential substring relationships between variable names:\n")
      print(utils::str(out_list))
  } else {
      cat("No simple substring relationships found between variable names.\n")
  }
  # Convert list to named character vector for original return type
  result <- unlist(out_list)
  if (is.null(result)) {
    return(character(0))
  }
  return(result)
}


#' ReplaceOutlierNumValsData
#'
#' @title ReplaceOutlierNumValsData
#'
#' @description Replace outlier values in numeric columns using the IQR method.
#' Outliers are defined as values below Q1 - multIQR * IQR or above Q3 + multIQR * IQR.
#' Replacement is done by capping at these lower and upper bounds.
#' Only applies to numeric columns with at least `minnumgroup` unique non-NA values.
#'
#' @param data Data frame.
#' @param multIQR Multiplier for the Interquartile Range (IQR) to define outlier bounds (default: 1.5).
#' @param minnumgroup Minimum number of unique non-NA values required for a numeric column to be processed (default: 10).
#' @return Data frame with outliers in numeric columns capped.
#' @importFrom stats quantile IQR
#' @export
ReplaceOutlierNumValsData<-function(data, multIQR=1.5, minnumgroup=10){
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Handle empty data frame
  if (ncol(data) == 0 || nrow(data) == 0) {
    return(data)
  }
  
  data_out <- data

  for (i in seq_len(ncol(data_out))){
    col_data <- data_out[[i]]
    col_name <- colnames(data_out)[i]

    if (is.numeric(col_data)){
      unique_vals <- unique(col_data[!is.na(col_data)])
      if (length(unique_vals) >= minnumgroup){
          qnt <- stats::quantile(col_data, probs=c(.25, .75), na.rm = TRUE)
          H <- multIQR * stats::IQR(col_data, na.rm = TRUE)
          lower_bound <- unname(qnt[1] - H)
          upper_bound <- unname(qnt[2] + H)

          outliers_low <- col_data < lower_bound & !is.na(col_data)
          outliers_high <- col_data > upper_bound & !is.na(col_data)

          if(any(outliers_low) || any(outliers_high)) {
              cat("Variable '", col_name, "': Capping", sum(outliers_low), "low and", sum(outliers_high), "high outliers.\n")
              col_data[outliers_low] <- unname(lower_bound)
              col_data[outliers_high] <- unname(upper_bound)
              data_out[[i]] <- col_data
          }
      }
    }
  }
  data_out
}


#' @title MakeTestDataConfWithTrainData
#'
#' @description Ensure test data conforms to training data structure and factor levels.
#' 1. Selects only columns present in training data.
#' 2. For factor columns, replaces values not seen in training with NA.
#' 3. Re-applies factor levels from training data to ensure consistency.
#'
#' @param traindata Data frame used for training.
#' @param testdata Data frame to be conformed.
#' @return Conformed test data frame.
#' @export
MakeTestDataConfWithTrainData<-function(traindata, testdata){
  if (!is.data.frame(traindata) || !is.data.frame(testdata)) {
      stop("Both 'traindata' and 'testdata' must be data frames.")
  }

  # Handle empty train data
  if (ncol(traindata) == 0) {
    return(data.frame())
  }

  train_cols <- colnames(traindata)
  test_cols <- colnames(testdata)
  common_cols <- intersect(train_cols, test_cols)
  missing_in_test <- setdiff(train_cols, test_cols)
  extra_in_test <- setdiff(test_cols, train_cols)

   if (length(extra_in_test) > 0) {
      cat("Removing columns from test data not present in training data:", paste(extra_in_test, collapse=", "), "\n")
  }

  # Subset test data to common columns, maintaining order of training data
  testdata_out <- testdata[, intersect(train_cols, common_cols), drop = FALSE]

  # Add missing columns filled with NA
  if (length(missing_in_test) > 0) {
      cat("Adding columns from training data missing in test data:", paste(missing_in_test, collapse=", "), "\n")
      for (col in missing_in_test) {
        # Create NA values with the same type as training column
        if (is.factor(traindata[[col]])) {
          testdata_out[[col]] <- factor(rep(NA, nrow(testdata_out)), levels = levels(traindata[[col]]))
        } else {
          testdata_out[[col]] <- rep(NA, nrow(testdata_out))
          # Ensure same type as training data
          if (is.numeric(traindata[[col]])) {
            testdata_out[[col]] <- as.numeric(testdata_out[[col]])
          } else if (is.character(traindata[[col]])) {
            testdata_out[[col]] <- as.character(testdata_out[[col]])
          } else if (is.logical(traindata[[col]])) {
            testdata_out[[col]] <- as.logical(testdata_out[[col]])
          }
        }
      }
  }

  # Ensure columns are in the same order as training data
  testdata_out <- testdata_out[, train_cols, drop = FALSE]
  for (vari in colnames(testdata_out)){
    if (is.factor(traindata[[vari]])){ # Check type in training data
      train_levels <- levels(traindata[[vari]])
      test_col_data <- testdata_out[[vari]]

      if (!is.factor(test_col_data)) {
          # If test column is not factor, try converting using train levels
          warning("Column '", vari, "' is factor in train data but not in test data. Attempting conversion.")
          test_col_data_char <- as.character(test_col_data)
          test_col_data_char[!(test_col_data_char %in% train_levels)] <- NA
          testdata_out[[vari]] <- factor(test_col_data_char, levels = train_levels)
      } else {
          # If test column is already factor, conform levels
          current_test_levels <- levels(test_col_data)
          new_vals_in_test <- setdiff(current_test_levels, train_levels)

          test_col_data_char <- as.character(test_col_data)
          # Set values with levels not in training data to NA
          test_col_data_char[test_col_data_char %in% new_vals_in_test] <- NA
          # Re-apply factor with training levels
          testdata_out[[vari]] <- factor(test_col_data_char, levels = train_levels)

          if (length(new_vals_in_test) > 0) {
              cat("Variable '", vari, "': Replaced levels not seen in training data with NA:", paste(new_vals_in_test, collapse=", "), "\n")
          }
      }
    }
    # Optional: Check numeric/character type consistency?
    # else if (class(traindata[[vari]]) != class(testdata_out[[vari]])) {
    #    warning("Type mismatch for column '", vari, "': Train='", class(traindata[[vari]]), "', Test='", class(testdata_out[[vari]]), "'")
    # }
  }
  testdata_out
}



#' @title RemoveEmptySpacesData
#'
#' @description Remove leading/trailing whitespace from character columns and factor levels.
#'
#' @param DATA Data frame.
#' @return Data frame with whitespace trimmed.
#' @export
RemoveEmptySpacesData<-function(DATA){
  if (!is.data.frame(DATA)) stop("'DATA' must be a data frame.")
  DATA_out <- DATA
  for (vari in colnames(DATA_out)){
    col_data <- DATA_out[[vari]]
    if (is.character(col_data)){
      DATA_out[[vari]] <- trimws(col_data, "both")
    } else if (is.factor(col_data)) {
      # Trim existing levels
      current_levels <- levels(col_data)
      trimmed_levels <- trimws(current_levels, "both")
      # Check if trimming creates duplicate levels
      if (any(duplicated(trimmed_levels))) {
          warning("Trimming whitespace for factor '", vari, "' creates duplicate levels. Consider collapsing levels first.")
          # Simple approach: just trim values and refactor (might merge levels)
          DATA_out[[vari]] <- factor(trimws(as.character(col_data), "both"))
      } else {
          levels(DATA_out[[vari]]) <- trimmed_levels
      }
    }
  }
  DATA_out
}
