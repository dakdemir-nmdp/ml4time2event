#' readDict
#'
#' Import a dictionary that is in standard format.
#' The dictionary should be in a csv, xlsx, or sas7bdat file.
#' Standard columns: Variable, Description, Type, Value, Label, Notes.
#'
#' @param file path to dictionary file.
#' @param filetype type of dictionary file ("csv", "xlsx", "sas7bdat").
#' @param ... additional arguments to pass to read.csv, openxlsx::read.xlsx,
#' haven::read_sas.
#' @return data frame with dictionary information.
#' @examples
#' \dontrun{
#' dict <- readDict(file = "path/to/dictionary.csv", filetype = "csv")
#' }
#' @importFrom utils read.csv
#' @importFrom openxlsx read.xlsx
#' @importFrom haven read_sas
#' @export
readDict <- function(file = NULL, filetype="csv", ...) {
  if (is.null(file)) {
    stop("Provide file location using the 'file' argument.")
  }
  if (!file.exists(file)) {
      stop("Specified dictionary file does not exist: ", file)
  }

  DICDATA <- NULL
  tryCatch({
      if (tolower(filetype) == "csv"){
        DICDATA <- as.data.frame(utils::read.csv(file, stringsAsFactors = FALSE, ...)) # Avoid factors initially
      } else if (tolower(filetype) == "xlsx") {
        if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' needed for xlsx files.")
        DICDATA <- as.data.frame(openxlsx::read.xlsx(file, ...))
      } else if (tolower(filetype) == "sas7bdat"){
        if (!requireNamespace("haven", quietly = TRUE)) stop("Package 'haven' needed for sas7bdat files.")
       DICDATA <- as.data.frame(haven::read_sas(file, ...))
      } else {
        stop("Invalid filetype provided. Choose one of 'csv', 'xlsx', 'sas7bdat'.")
      }
  }, error = function(e) {
      stop("Error reading dictionary file '", file, "' with filetype '", filetype, "': ", e$message)
  })

  required_cols <- c("Variable", "Description", "Type", "Value", "Label", "Notes")
  missing_cols <- setdiff(required_cols, colnames(DICDATA))

  if (length(missing_cols) > 0) {
      stop(
        "Dictionary file is missing required columns: ", paste(missing_cols, collapse=", "), ".\n",
        "Required columns are: Variable, Description, Type, Value, Label, Notes."
      )
  }
  # Select and reorder columns to standard format
  DICDATA <- DICDATA[, required_cols]

  # Basic type conversion attempt (can be refined)
  DICDATA$Value <- as.character(DICDATA$Value)
  DICDATA$Label <- as.character(DICDATA$Label)
  DICDATA$Variable <- as.character(DICDATA$Variable)
  DICDATA$Description <- as.character(DICDATA$Description)
  DICDATA$Type <- as.character(DICDATA$Type)
  DICDATA$Notes <- as.character(DICDATA$Notes)


  return(DICDATA)
}



#' filllasttoDict
#'
#' Utility function to fill missing values (NA) in dictionary columns
#' by carrying forward the last non-NA value.
#'
#' @param dict data frame with dictionary information.
#' @return data frame with NA values filled forward column-wise.
#' @export
filllasttoDict <- function(dict) {
  # Helper function to fill NAs forward in a single vector
  filllastvector <- function(x) {
    y <- x
    filli <- NA # Stores the last non-NA value encountered
    for (i in 1:length(y)) {
      if (!is.na(y[i]) && y[i] != "") { # Consider empty strings as missing too? Optional.
        filli <- y[i]
      } else if (!is.na(filli)) { # Only fill if we have a value to carry forward
        y[i] <- filli
      }
    }
    return(y) # No need for unlist here
  }
  # Apply the function to each column
  DICT_filled <- as.data.frame(lapply(dict, filllastvector), stringsAsFactors = FALSE)
  colnames(DICT_filled) <- colnames(dict) # Preserve original column names
  return(DICT_filled)
}



#' bindDicts
#'
#' Utility function to combine two dictionaries. Rows from Dict2 for variables
#' already present in Dict1 are ignored.
#'
#' @param Dict1 data frame, the primary dictionary.
#' @param Dict2 data frame, the secondary dictionary to append from.
#' @return data frame containing combined dictionary information.
#' @importFrom dplyr bind_rows
#' @export
bindDicts <- function(Dict1, Dict2) {
  if (!is.data.frame(Dict1) || !is.data.frame(Dict2)) {
      stop("Both Dict1 and Dict2 must be data frames.")
  }
  if (!"Variable" %in% colnames(Dict1) || !"Variable" %in% colnames(Dict2)) {
      stop("Both dictionaries must contain a 'Variable' column.")
  }

  Names1 <- unique(as.character(Dict1$Variable))
  Names2 <- unique(as.character(Dict2$Variable))

  # Keep only rows from Dict2 where the variable is NOT in Dict1
  Dict2rem <- Dict2[!(Dict2$Variable %in% Names1), ]

  # Combine Dict1 and the remaining rows from Dict2
  DICT <- dplyr::bind_rows(Dict1, Dict2rem)
  # Reset row names for cleanliness
  rownames(DICT) <- NULL
  return(DICT)
}


#' CatVarsDict
#'
#' Utility function to subset dictionary for categorical variables.
#'
#' @param dict data frame with dictionary information (must have 'Type' column).
#' @return data frame containing only rows where Type is "Categorical".
#' @export
CatVarsDict <- function(dict) {
  if (!"Type" %in% colnames(dict)) stop("Dictionary must contain a 'Type' column.")
  return(dict[dict$Type == "Categorical", , drop = FALSE])
}


#' ConVarsDict
#'
#' Utility function to subset dictionary for continuous variables.
#' @param dict data frame with dictionary information (must have 'Type' column).
#' @return data frame containing only rows where Type is "Continuous".
#'
#' @export
ConVarsDict <- function(dict) {
  if (!"Type" %in% colnames(dict)) stop("Dictionary must contain a 'Type' column.")
  return(dict[dict$Type == "Continuous", , drop = FALSE])
}


#' CharVarsDict
#'
#' Utility function to subset dictionary for character variables.
#' @param dict data frame with dictionary information (must have 'Type' column).
#' @return data frame containing only rows where Type is "Character".
#'
#' @export
CharVarsDict <- function(dict) {
  if (!"Type" %in% colnames(dict)) stop("Dictionary must contain a 'Type' column.")
  return(dict[dict$Type == "Character", , drop = FALSE])
}


#' OneValCat2NumDict
#'
#' Utility function to correct the 'Type' for variables that appear only once
#' in the dictionary (implying continuous) but are marked as 'Categorical'.
#' This assumes single-row entries correspond to continuous variables.
#'
#' @param dict data frame with dictionary information.
#' @param vars character vector of variable names to potentially correct. If NULL, checks all variables.
#' @return data frame with corrected 'Type' information.
#' @export
OneValCat2NumDict <- function(dict, vars = NULL) {
  if (!"Variable" %in% colnames(dict) || !"Type" %in% colnames(dict)) {
      stop("Dictionary must contain 'Variable' and 'Type' columns.")
  }
  if (is.null(vars)) {
      vars_to_check <- unique(dict$Variable)
  } else {
      vars_to_check <- intersect(vars, unique(dict$Variable))
  }

  # Count occurrences of each variable name in the dictionary
  tabx <- table(dict$Variable)
  # Identify variables appearing only once
  single_occurrence_vars <- names(tabx)[tabx == 1]

  # Find variables within the target set ('vars_to_check') that appear only once
  # AND are currently marked as 'Categorical'
  vars_to_correct <- intersect(intersect(vars_to_check, single_occurrence_vars),
                               dict$Variable[dict$Type == "Categorical"])

  # Correct the type for these variables
  if (length(vars_to_correct) > 0) {
      dict[dict$Variable %in% vars_to_correct, "Type"] <- "Continuous"
      cat("Corrected type to 'Continuous' for single-occurrence 'Categorical' variables:", paste(vars_to_correct, collapse=", "), "\n")
  }
  return(dict)
}


#' NumIdVarstoCharDict
#'
#' Utility function to change the 'Type' of specified ID variables to "Character".
#'
#' @param dict data frame with dictionary information.
#' @param idvars character vector of variable names to treat as character IDs.
#' @return data frame with updated 'Type' for ID variables.
#' @export
NumIdVarstoCharDict <- function(dict, idvars = c("`#`", "pseudoid", "pseudoccn")) {
   if (!"Variable" %in% colnames(dict) || !"Type" %in% colnames(dict)) {
      stop("Dictionary must contain 'Variable' and 'Type' columns.")
  }
  if (!is.null(idvars) && length(idvars) > 0) {
      # Find which idvars are actually present in the dictionary
      vars_in_dict <- idvars[idvars %in% dict$Variable]
      if (length(vars_in_dict) > 0) {
          dict[dict$Variable %in% vars_in_dict, "Type"] <- "Character"
          cat("Set type to 'Character' for ID variables:", paste(vars_in_dict, collapse=", "), "\n")
      } else {
          warning("None of the specified idvars found in the dictionary.")
      }
  }
  return(dict)
}



#' DateVarstoDateDict
#'
#' Utility function to change the 'Type' of specified date variables to "Date".
#' @param dict data frame with dictionary information.
#' @param datevars character vector of variable names to mark as "Date".
#' @return data frame with updated 'Type' for date variables.
#'
#' @export
DateVarstoDateDict <- function(dict, datevars = NULL) {
  if (!"Variable" %in% colnames(dict) || !"Type" %in% colnames(dict)) {
      stop("Dictionary must contain 'Variable' and 'Type' columns.")
  }
  if (!is.null(datevars) && length(datevars) > 0) {
      vars_in_dict <- datevars[datevars %in% dict$Variable]
       if (length(vars_in_dict) > 0) {
            dict[dict$Variable %in% vars_in_dict, 'Type'] <- "Date"
            cat("Set type to 'Date' for variables:", paste(vars_in_dict, collapse=", "), "\n")
       } else {
           warning("None of the specified datevars found in the dictionary.")
       }
  }
  return(dict)
}


#' addDiseaseSpecColsDict
#'
#' Process the dictionary to add columns indicating disease-specific applicability.
#' Adds 'NACats' column marking 'NotApp' for specific values (e.g., ".E").
#' Adds 'DisSpec' column listing diseases for which a variable might be specific,
#' based on labels like "N/A, other disease".
#'
#' @param dict data frame with dictionary information (must include Value, Label, Description, Variable columns).
#' @return data frame with added 'NACats' and 'DisSpec' columns.
#' @export
addDiseaseSpecColsDict<-function(dict){
  required_cols <- c("Variable", "Description", "Value", "Label")
  if (!all(required_cols %in% colnames(dict))) {
      stop("Dictionary must contain columns: ", paste(required_cols, collapse=", "))
  }

  # Initialize NACats if it doesn't exist
  if (!"NACats" %in% colnames(dict)) {
      dict$NACats <- "VALUE" # Default value
  }
  # Mark specific values as Not Applicable
  dict$NACats[dict$Value %in% ".E"] <- "NotApp" # Example: ".E" means Not Applicable

  # Identify potentially disease-specific variables based on Label
  # This logic assumes labels like "N/A, other disease" indicate specificity
  na_other_disease_label <- "N/A, other disease" # Define the target label
  DiseaseSpecificVars <- unique(dict$Variable[dict$Label %in% na_other_disease_label])

  # Get the list of all diseases (assuming a variable named 'disease' exists with disease labels)
  if ("disease" %in% dict$Variable) {
      alldiseases <- unique(dict$Label[dict$Variable == "disease"])
      alldiseases <- alldiseases[!is.na(alldiseases) & alldiseases != ""] # Clean disease list
  } else {
      warning("No 'disease' variable found in dictionary to identify disease list for specificity check.")
      alldiseases <- character(0)
  }

  # Initialize DisSpec column
  dict$DisSpec <- NA_character_

  if (length(DiseaseSpecificVars) > 0 && length(alldiseases) > 0) {
      # Get descriptions associated with the "N/A, other disease" label for potentially specific vars
      DiseaseSpecificVarsDescr <- unique(dict$Description[dict$Variable %in% DiseaseSpecificVars & dict$Label %in% na_other_disease_label])

      # For each disease, find descriptions that mention it and update DisSpec
      for (i in 1:length(alldiseases)) {
          disease_name <- alldiseases[i]
          # Find descriptions containing the current disease name (case-insensitive search)
          matching_descr_indices <- which(grepl(disease_name, DiseaseSpecificVarsDescr, ignore.case = TRUE))

          if (length(matching_descr_indices) > 0) {
              matching_descriptions <- DiseaseSpecificVarsDescr[matching_descr_indices]
              # Find rows in the original dictionary matching these descriptions
              dict_rows_to_update <- which(dict$Description %in% matching_descriptions)

              # Update the DisSpec column for these rows
              for (row_idx in dict_rows_to_update) {
                  current_spec <- dict$DisSpec[row_idx]
                  # Append disease name if not already present
                  if (is.na(current_spec)) {
                      dict$DisSpec[row_idx] <- disease_name
                  } else if (!grepl(paste0("\\b", disease_name, "\\b"), current_spec)) { # Avoid duplicates
                      dict$DisSpec[row_idx] <- paste(current_spec, disease_name, sep = " ")
                  }
              }
          }
      }
  }
  dict
}


#' MakelabelledSASDict
#'
#' Create a dictionary data frame from a labelled SAS dataset file.
#' Uses the 'labelled' package to generate the dictionary.
#'
#' @param fileloc Path to the SAS dataset file (.sas7bdat).
#' @return A data frame representing the dictionary derived from SAS labels.
#' @importFrom labelled generate_dictionary
#' @importFrom haven read_sas
#' @export
MakelabelledSASDict<-function(fileloc){
  if (!requireNamespace("labelled", quietly = TRUE)) stop("Package 'labelled' needed for this function.")
  if (!requireNamespace("haven", quietly = TRUE)) stop("Package 'haven' needed for this function.")
  if (!file.exists(fileloc)) stop("Specified SAS file does not exist: ", fileloc)

  dataraw <- tryCatch(
      haven::read_sas(fileloc),
      error = function(e) stop("Error reading SAS file '", fileloc, "': ", e$message)
  )
  dictionary <- labelled::generate_dictionary(dataraw)
  # Convert to standard data frame if needed (generate_dictionary might return tibble)
  dictionary <- as.data.frame(dictionary)
  # Potentially rename columns to match standard format (Variable, Label, etc.) if needed
  # Example renaming (adjust based on actual output of generate_dictionary):
  # colnames(dictionary)[colnames(dictionary) == "variable"] <- "Variable"
  # colnames(dictionary)[colnames(dictionary) == "label"] <- "Description"
  # ... add more renames and potentially create missing standard columns ...
  return(dictionary)
}
