#' DetectIdVarsDict
#'
#' Utility function to detect potential ID variables in a dictionary based on common patterns in variable names.
#'
#' @param dict data frame with dictionary information (must have 'Variable' column).
#' @param idvars Optional character vector of known ID variable names to include.
#' @return character vector of potential ID variable names found in the dictionary.
#' @importFrom stringr str_detect
#' @export
DetectIdVarsDict <- function(dict, idvars = NULL) {
  if (!"Variable" %in% colnames(dict)) {
      stop("Dictionary must contain a 'Variable' column.")
  }
  variables <- unique(as.character(dict$Variable)) # Ensure character vector

  # Define patterns for common ID variable names (case-insensitive)
  id_patterns <- c("id$", "^id", "ccn", "pseudo", "team", "case", "identif") # Added end/start anchors and 'identifier' pattern

  # Detect variables matching the patterns
  detected_vars <- variables[sapply(id_patterns, function(pattern) {
      grep(pattern, variables, ignore.case = TRUE)
  }) |> unlist() |> unique()] # Combine results from all patterns

  # Combine detected variables with user-provided known ID variables
  potential_idvars <- union(detected_vars, idvars)

  # Return only those potential ID vars that are actually in the dictionary
  return(intersect(potential_idvars, variables))
}


#' DetectNAStringsDict
#'
#' Utility function to detect potential NA representations in dictionary 'Label' and 'Value' columns.
#' Adds an 'NACats' column to the dictionary, marking rows as "NAVAL" if suspected NA, otherwise "VALUE".
#'
#' @param dict data frame with dictionary information (must have 'Label' and 'Value' columns).
#' @param knownNaLabels Optional character vector of additional strings in 'Label' column to consider NA.
#' @param knownNaValues Optional character vector of additional strings in 'Value' column to consider NA.
#' @return data frame with added 'NACats' column ("NAVAL" or "VALUE").
#' @importFrom stringr str_detect
#' @export
DetectNAStringsDict <-
  function(dict,
           knownNaLabels = NULL,
           knownNaValues = NULL) {
    if (!"Label" %in% colnames(dict) || !"Value" %in% colnames(dict)) {
        stop("Dictionary must contain 'Label' and 'Value' columns.")
    }

    # Define default patterns/strings indicating NA (case-insensitive)
    default_nalabels <-
      c("missing", "unknown", "not answ", "not answered", "n/a", "unspecified", "nan", "declined", "refused")
    default_navalues <-
      c("^\\.$", "^\\.b$", "^\\.c$", "^\\.f$", "^\\.y$", "^\\.$", # Specific SAS missing codes (anchored)
        "<NA>", "99", "999", "88", "888", # Common numeric codes
        "null", "na", "" # Common string codes (empty string)
        )

    # Combine default and user-provided NA indicators
    nalabels <- union(default_nalabels, tolower(knownNaLabels))
    navalues <- union(default_navalues, tolower(knownNaValues))

    # Create regex patterns for matching (case-insensitive)
    # Match whole words/codes where appropriate
    label_pattern <- paste0("\\b(", paste(nalabels, collapse = "|"), ")\\b")
    # For values, match exact strings or specific patterns
    value_pattern <- paste0("^(", paste(navalues, collapse = "|"), ")$") # Anchor values

    # Find rows where Label or Value suggests NA
    label_na_indices <- which(stringr::str_detect(tolower(dict$Label), label_pattern))
    value_na_indices <- which(stringr::str_detect(tolower(dict$Value), value_pattern))

    # Combine indices
    na_indices <- unique(c(label_na_indices, value_na_indices))

    # Add or update the NACats column
    dict$NACats <- "VALUE" # Default to VALUE
    if (length(na_indices) > 0) {
        dict$NACats[na_indices] <- "NAVAL" # Mark suspected NAs
    }

    cat("Marked", length(na_indices), "rows as potential NA values ('NAVAL').\n")
    return(dict)
  }


#' finddatevars
#'
#' Utility function to find potential date variables based on common patterns in names.
#'
#' @param names character vector of variable names.
#' @param knowndatevars Optional character vector of known date variable names to include.
#' @return character vector of potential date variable names found in 'names'.
#' @export
finddatevars<-function(names, knowndatevars=NULL){
  # Define patterns for common date variable names (case-insensitive)
  date_patterns <- c("dte", "date", "_dt$", "^dt_") # Common substrings and endings/beginnings

  # Detect variables matching the patterns
  detected_vars <- names[sapply(date_patterns, function(pattern) {
      grep(pattern, names, ignore.case = TRUE)
  }) |> unlist() |> unique()]

  # Combine detected variables with user-provided known date variables
  potential_datevars <- union(detected_vars, knowndatevars)

  # Return only those potential date vars that are actually in the input 'names'
  return(intersect(potential_datevars, names))
}
