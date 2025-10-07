#' ApplyDictToData
#'
#' Apply dictionary specifications (types, labels, NA values) to a data frame.
#' Converts columns to specified types (Factor, Character, Numeric, Date).
#' Applies factor levels and labels based on dictionary entries.
#' Sets values identified as NA in the dictionary (using 'NACats' column) to NA in the data.
#'
#' @param dict Data frame dictionary with columns: Variable, Type, Value, Label, NACats (optional, created by DetectNAStringsDict).
#' @param data Data frame to process. Column names should match 'Variable' in dict (case-insensitive matching attempted).
#' @return Data frame with types converted, factors labelled, and NAs applied according to the dictionary.
#' @importFrom lubridate parse_date_time
#' @export
ApplyDictToData <- function(dict, data) {
  if (!is.data.frame(dict) || !is.data.frame(data)) {
      stop("Both 'dict' and 'data' must be data frames.")
  }
  required_dict_cols <- c("Variable", "Type", "Value", "Label")
  if (!all(required_dict_cols %in% colnames(dict))) {
      stop("Dictionary missing required columns: ", paste(setdiff(required_dict_cols, colnames(dict)), collapse=", "))
  }
  # Ensure NACats exists, default to "VALUE" if not
  if (!"NACats" %in% colnames(dict)) {
      warning("Dictionary missing 'NACats' column. Assuming all values are valid (not NA). Consider running DetectNAStringsDict first.")
      dict$NACats <- "VALUE"
  }

  data_out <- as.data.frame(data) # Work on a copy
  dict_vars <- unique(dict$Variable)
  data_cols_lower <- tolower(colnames(data_out))

  for (vari in dict_vars) {
    dictx <- dict[dict$Variable == vari, ]
    # Find matching column in data (case-insensitive)
    match_idx <- which(data_cols_lower == tolower(vari))

    if (length(match_idx) == 1) {
      VarName <- colnames(data_out)[match_idx] # Get original case name
      col_data <- data_out[[VarName]] # Extract column data

      # --- 1. Apply NA values based on dictionary ---
      na_values <- dictx$Value[dictx$NACats == "NAVAL"]
      if (length(na_values) > 0) {
          # Handle potential type mismatches before comparison
          tryCatch({
              if (is.numeric(col_data)) {
                  na_values_num <- suppressWarnings(as.numeric(na_values))
                  col_data[col_data %in% na_values_num[!is.na(na_values_num)]] <- NA
              } else { # Character or Factor
                  col_data[trimws(as.character(col_data)) %in% trimws(na_values)] <- NA
              }
          }, warning = function(w) {
              warning("Potential type mismatch when applying NA values for variable '", VarName, "': ", w$message)
          })
      }

      # --- 2. Convert to specified Type ---
      target_type <- dictx$Type[[1]]
      current_type <- class(col_data)

      tryCatch({
          if (target_type == "Categorical") {
              valid_entries <- dictx[dictx$NACats != "NAVAL", ]
              Values <- trimws(valid_entries$Value)
              Labels <- trimws(valid_entries$Label)
              # Ensure unique values map to unique labels if possible
              if (length(Values) != length(Labels)) {
                  warning("Mismatch between number of Values and Labels for Categorical variable '", VarName, "'. Using Values as Labels.")
                  Labels <- Values
              }
              # Remove duplicates carefully, keeping first occurrence
              unique_idx <- !duplicated(Values)
              Values <- Values[unique_idx]
              Labels <- Labels[unique_idx]

              if (length(Values) > 0) {
                  # Convert data column to character before applying factor levels
                  col_data_char <- trimws(as.character(col_data))
                  # Apply factor levels and labels
                  col_data <- factor(col_data_char, levels = Values, labels = Labels)
              } else {
                  warning("No valid levels defined for Categorical variable '", VarName, "'. Converting to character.")
                  col_data <- as.character(col_data)
              }
          } else if (target_type == "Character") {
              col_data <- as.character(col_data)
          } else if (target_type == "Continuous" || target_type == "Numeric") {
              # Attempt numeric conversion, warn if NAs introduced
              original_na_count <- sum(is.na(col_data))
              col_data_num <- suppressWarnings(as.numeric(as.character(col_data))) # Convert via character to handle factors
              new_na_count <- sum(is.na(col_data_num))
              if (new_na_count > original_na_count) {
                  warning("NAs introduced by coercion to numeric for variable '", VarName, "'")
              }
              col_data <- col_data_num
          } else if (target_type == "Date") {
              # Attempt date parsing using lubridate (more flexible)
              if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package 'lubridate' needed for Date conversion.")
              # Try common formats; add more as needed
              col_data <- lubridate::parse_date_time(as.character(col_data), orders = c("Ymd", "mdY", "dmY", "Y-m-d", "m/d/Y", "d/m/Y", "b d Y", "d b Y"))
              if (all(is.na(col_data[!is.na(data_out[[VarName]])]))) { # Check if parsing failed for non-NA original values
                  warning("Date parsing failed for variable '", VarName, "'. Check format or dictionary type.")
              } else {
                   col_data <- as.Date(col_data) # Convert to Date object
              }
          } else {
              warning("Unsupported dictionary Type '", target_type, "' for variable '", VarName, "'. Column not converted.")
          }
      }, error = function(e) {
          warning("Error processing variable '", VarName, "' with type '", target_type, "': ", e$message)
      })


      # --- 3. Update data frame and add comment ---
      data_out[[VarName]] <- col_data
      comment(data_out[[VarName]]) <- dictx$Description[[1]]

    } else if (length(match_idx) > 1) {
        warning("Variable '", vari, "' matched multiple columns in data (case-insensitive). Skipping.")
    } else {
       # warning("Variable '", vari, "' from dictionary not found in data.") # Reduce noise?
    }
  }
  return(data_out)
}


#' RemoveMonoVarsData
#'
#' Remove variables (columns) from a data frame that have only one unique non-NA value.
#'
#' @param data Data frame.
#' @return Data frame with single-value columns removed.
#' @export
RemoveMonoVarsData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  num_unique_vals <- sapply(data, function(x) length(unique(x[!is.na(x)])))
  cols_to_keep <- num_unique_vals > 1
  if (sum(!cols_to_keep) > 0) {
      cat("Removing", sum(!cols_to_keep), "columns with only one unique non-NA value:", paste(colnames(data)[!cols_to_keep], collapse=", "), "\n")
  }
  return(data[, cols_to_keep, drop = FALSE])
}

#' RemoveAllNAVars
#'
#' Remove variables (columns) from a data frame that consist entirely of NA values.
#'
#' @param data Data frame.
#' @return Data frame with all-NA columns removed.
#' @export
RemoveAllNAVars <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  all_na_cols <- sapply(data, function(x) all(is.na(x)))
  if (sum(all_na_cols) > 0) {
      cat("Removing", sum(all_na_cols), "columns containing only NA values:", paste(colnames(data)[all_na_cols], collapse=", "), "\n")
  }
  return(data[, !all_na_cols, drop = FALSE])
}



#' findDiseaseSpecVarsData
#'
#' Utility function to heuristically detect disease-specific variables from data.
#' A variable is considered specific to a disease if it's mostly non-NA within
#' that disease group and mostly NA outside that group.
#' Assumes a 'disease' column exists in the data.
#'
#' @param data Data frame containing a 'disease' column and potential disease-specific variables.
#' @param pin Numeric threshold (0-1) for the minimum proportion of non-NA values within the disease group (default: 0.5).
#' @param pout Numeric threshold (0-1) for the minimum proportion of NA values outside the disease group (default: 0.9).
#' @return A named list where names are disease levels and values are character vectors of variables deemed specific to that disease.
#' @export
findDiseaseSpecVarsData <- function(data, pin = .5, pout = .9) {
  if (!"disease" %in% colnames(data)) {
      stop("Data frame must contain a 'disease' column.")
  }
  varstocheck <- setdiff(colnames(data), "disease")
  diseases <- unique(data$disease[!is.na(data$disease)])
  dislist <- vector(mode = "list")

  for (disi in diseases) {
    varlist <- c()
    rows_disease <- which(data$disease == disi)
    rows_other <- which(data$disease != disi)
    n_disease <- length(rows_disease)
    n_other <- length(rows_other)

    if (n_disease == 0 || n_other == 0) next # Skip if no obs in/out of group

    for (vari in varstocheck) {
      col_data <- data[[vari]]
      na_indices <- is.na(col_data)
      non_na_indices <- !na_indices

      # Proportion non-NA within the disease group
      prop_non_na_in <- sum(non_na_indices[rows_disease]) / n_disease
      # Proportion NA outside the disease group
      prop_na_out <- sum(na_indices[rows_other]) / n_other

      # Check thresholds
      if (!is.nan(prop_non_na_in) && !is.nan(prop_na_out) && prop_non_na_in > pin && prop_na_out > pout) {
        varlist <- c(varlist, vari)
      }
    }

    if (length(varlist) > 0) {
      dislist[[as.character(disi)]] <- varlist # Use disease level as name
    }
  }
  if (length(dislist) > 0) {
      cat("Found potentially disease-specific variables:\n")
      print(utils::str(dislist))
  } else {
      cat("No potentially disease-specific variables found based on thresholds.\n")
  }
  dislist
}




#' processDiseaseSpecVarsData
#'
#' Utility function to process potentially disease-specific variables.
#' For variables identified as specific to a disease (via `dislist`),
#' replaces NA values outside that disease group with "NotApp" and converts the column to factor.
#'
#' @param data Data frame to process.
#' @param dislist Named list (output of `findDiseaseSpecVarsData`) mapping diseases to specific variable names.
#' @return Data frame with NAs replaced by "NotApp" for specific variables/disease combinations.
#' @export
processDiseaseSpecVarsData <- function(data, dislist) {
  if (!is.list(dislist) || is.null(names(dislist))) {
      stop("'dislist' must be a named list (e.g., output from findDiseaseSpecVarsData).")
  }
  if (!"disease" %in% colnames(data)) {
      stop("Data frame must contain a 'disease' column.")
  }

  data_out <- data # Work on a copy

  for (disease in names(dislist)) {
    vars <- dislist[[disease]]
    rows_not_disease <- which(data_out$disease != disease)

    for (vari in vars) {
      if (vari %in% colnames(data_out)) {
          col_data <- data_out[[vari]]
          # Convert to character to allow adding "NotApp" level
          col_data_char <- as.character(col_data)
          # Replace NA for rows *not* in the specific disease group
          na_in_other_rows <- is.na(col_data_char[rows_not_disease])
          col_data_char[rows_not_disease][na_in_other_rows] <- "NotApp"
          # Convert back to factor, ensuring "NotApp" is a level
          data_out[[vari]] <- factor(col_data_char, levels = unique(c(levels(col_data), "NotApp")))
          cat("Processed variable '", vari, "' for disease '", disease, "', replacing NAs outside group with 'NotApp'.\n")
      } else {
          warning("Variable '", vari, "' listed for disease '", disease, "' not found in data.")
      }
    }
  }
  return(data_out)
}




#' RemoveMissinVarsData
#'
#' Remove variables (columns) with a proportion of missing values exceeding a threshold.
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
      cat("Removing", length(cols_removed), "columns with >", maxprop*100, "% missing values:", paste(cols_removed, collapse=", "), "\n")
  }
  return(data[, cols_to_keep, drop = FALSE])
}

#' RemoveMissinRecordsData
#'
#' Remove records (rows) with a proportion of missing values exceeding a threshold.
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
      cat("Removing", rows_removed_count, "rows with >", maxprop*100, "% missing values.\n")
  }
  return(data[rows_to_keep, , drop = FALSE])
}



#' getcharcolsData
#'
#' Utility function to get the names of character columns in a data frame.
#'
#' @param data Data frame.
#' @return Character vector with names of character columns.
#' @export
getcharcolsData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  charcols_logical <- sapply(data, is.character)
  return(names(which(charcols_logical)))
}


#' ImputeMissinRecordsData
#'
#' Impute missing values in a data frame using missRanger.
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




#' RemoveRareCategoriesData
#'
#' Replace rare categories in factor variables with NA.
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


#' RemoveRareBinaryVarsData
#'
#' For binary variables (numeric 0/1, logical, or 2-level factors), replace the rarer category with NA if its frequency is below a threshold.
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

       if (length(unique_non_na) == 2) { # Check if binary (ignoring NAs)
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


#' CollapseRareCategoriestoOtherData
#'
#' Collapse rare categories in factor or character variables into a new level called "Other".
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



#' droplevelsoffactorsData
#'
#' Drop unused factor levels from all factor columns in a data frame.
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

#' findvarsnamesthatrepeatData
#'
#' Find variable names where one name is a substring of another (potential redundancy).
#' This is a simple substring check, not checking for semantic similarity.
#'
#' @param data Data frame.
#' @return A named character vector where names are the shorter variable names and
#'   values are comma-separated lists of longer variable names containing the shorter name.
#' @export
findvarsnamesthatrepeatData <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  col_names <- colnames(data)
  out_list <- list()

  for (i in 1:length(col_names)) {
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
  unlist(out_list)
}


#' ReplaceOutlierNumValsData
#'
#' Replace outlier values in numeric columns using the IQR method.
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
  data_out <- data

  for (i in 1:ncol(data_out)){
    col_data <- data_out[[i]]
    col_name <- colnames(data_out)[i]

    if (is.numeric(col_data)){
      unique_vals <- unique(col_data[!is.na(col_data)])
      if (length(unique_vals) >= minnumgroup){
          qnt <- stats::quantile(col_data, probs=c(.25, .75), na.rm = TRUE)
          H <- multIQR * stats::IQR(col_data, na.rm = TRUE)
          lower_bound <- qnt[1] - H
          upper_bound <- qnt[2] + H

          outliers_low <- col_data < lower_bound & !is.na(col_data)
          outliers_high <- col_data > upper_bound & !is.na(col_data)

          if(any(outliers_low) || any(outliers_high)) {
              cat("Variable '", col_name, "': Capping", sum(outliers_low), "low and", sum(outliers_high), "high outliers.\n")
              col_data[outliers_low] <- lower_bound
              col_data[outliers_high] <- upper_bound
              data_out[[i]] <- col_data
          }
      }
    }
  }
  data_out
}


#' MakeTestDataConfWithTrainData
#'
#' Ensure test data conforms to training data structure and factor levels.
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

  train_cols <- colnames(traindata)
  test_cols <- colnames(testdata)
  common_cols <- intersect(train_cols, test_cols)
  missing_in_test <- setdiff(train_cols, test_cols)
  extra_in_test <- setdiff(test_cols, train_cols)

  if (length(missing_in_test) > 0) {
      warning("Columns from training data missing in test data: ", paste(missing_in_test, collapse=", "), ". They cannot be included.")
      # Or potentially add them as NA columns? Depends on requirement.
  }
   if (length(extra_in_test) > 0) {
      cat("Removing columns from test data not present in training data:", paste(extra_in_test, collapse=", "), "\n")
  }

  # Subset test data to common columns, maintaining order of training data
  testdata_out <- testdata[, intersect(train_cols, common_cols), drop = FALSE]

  # Conform factor levels
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



#' RemoveEmptySpacesData
#'
#' Remove leading/trailing whitespace from character columns and factor levels.
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
