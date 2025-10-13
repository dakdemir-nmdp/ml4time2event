#' @title pairwiserelationshipsDataSummmary
#'
#' @description Calculate pairwise non-linear correlations among variables using the
#' bias-corrected distance correlation (bcdcor) from the 'energy' package.
#' Note: This can be computationally intensive for large datasets or many variables.
#'
#' @param data Data frame with variables to correlate. Factors/characters are converted to dummy variables.
#' @return A symmetric matrix of pairwise distance correlations.
#' @importFrom energy bcdcor
#' @importFrom stats model.matrix
#' @importFrom pbapply pblapply
#' @export
pairwiserelationshipsDataSummmary <- function(data) {
  if (!requireNamespace("energy", quietly = TRUE)) stop("Package 'energy' needed.")
  if (!requireNamespace("pbapply", quietly = TRUE)) stop("Package 'pbapply' needed for progress bar.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")

  n_vars <- ncol(data)
  if (n_vars < 2) stop("Need at least two columns to calculate pairwise correlations.")

  # Convert factors/characters to numeric matrix using model.matrix
  # Handle potential errors during conversion
  tryCatch({
      # Using formula ~ . - 1 to get dummy variables without intercept
      data_matrix <- stats::model.matrix(~ . - 1, data = data)
      # Store original-like column names (model.matrix can alter them)
      matrix_colnames <- colnames(data_matrix)
  }, error = function(e) {
      stop("Failed to convert data to numeric matrix using model.matrix: ", e$message)
  })

  n_matrix_vars <- ncol(data_matrix)
  Cormat <- matrix(0, nrow = n_matrix_vars, ncol = n_matrix_vars)
  colnames(Cormat) <- rownames(Cormat) <- matrix_colnames

  cat("Calculating pairwise distance correlations for", n_matrix_vars, "numeric features...\n")

  # Use pblapply for parallel processing with progress bar
  results <- pbapply::pblapply(1:(n_matrix_vars - 1), function(i) {
      cor_row <- numeric(n_matrix_vars) # Initialize row for correlations
      xi <- data_matrix[, i, drop = FALSE] # Keep as matrix
      for (j in (i + 1):n_matrix_vars) {
          xj <- data_matrix[, j, drop = FALSE] # Keep as matrix
          # Calculate bcdcor, handle potential errors
          dcorij <- tryCatch(
              energy::bcdcor(xi, xj),
              error = function(e) {
                  warning("bcdcor failed for columns '", matrix_colnames[i], "' and '", matrix_colnames[j], "': ", e$message)
                  return(NA) # Return NA on error
              }
          )
          cor_row[j] <- dcorij
      }
      return(cor_row) # Return correlations for row i (upper triangle)
  })

  # Fill the correlation matrix from the results list
  for (i in 1:(n_matrix_vars - 1)) {
      Cormat[i, ] <- results[[i]]
  }

  # Make the matrix symmetric and set diagonal to 1
  Cormat <- Cormat + t(Cormat)
  diag(Cormat) <- 1

  return(Cormat)
}




#' @title gethighcorvarsDataSummmary
#'
#' @description Identify pairs of variables from a correlation matrix with correlation exceeding a threshold.
#'
#' @param pmat A symmetric correlation matrix (e.g., output from `pairwiserelationshipsDataSummmary` or standard cor()).
#' @param corcutoff Numeric correlation cutoff threshold (default: 0.8).
#' @return A matrix with three columns: Var1, Var2, Correlation for pairs exceeding the cutoff. Returns NULL if no pairs found.
#' @export
gethighcorvarsDataSummmary <- function(pmat, corcutoff = .8) {
  if (!is.matrix(pmat) || !isSymmetric(pmat)) stop("'pmat' must be a symmetric matrix.")
  if (!is.numeric(corcutoff) || corcutoff < -1 || corcutoff > 1) stop("'corcutoff' must be between -1 and 1.")

  pmat[lower.tri(pmat, diag = TRUE)] <- NA # Set lower triangle and diagonal to NA to avoid duplicates
  high_cor_indices <- which(abs(pmat) > corcutoff, arr.ind = TRUE) # Use absolute correlation

  if (nrow(high_cor_indices) == 0) {
    cat("No variable pairs found with absolute correlation >", corcutoff, "\n")
    return(NULL)
  }

  Var1 <- rownames(pmat)[high_cor_indices[, 1]]
  Var2 <- colnames(pmat)[high_cor_indices[, 2]]
  Correlation <- pmat[high_cor_indices]

  corvarsmat <- data.frame(Var1, Var2, Correlation, stringsAsFactors = FALSE)
  cat("Found", nrow(corvarsmat), "variable pairs with absolute correlation >", corcutoff, "\n")
  return(corvarsmat) # Return as data.frame
}


#' OneAgainstRestCorDataSummmary
#'
#' @title OneAgainstRestCorDataSummmary
#'
#' @description Calculate the distance correlation of each variable against all other variables combined.
#' Note: This can be computationally intensive.
#'
#' @param data Data frame. Factors/characters are converted to dummy variables.
#' @return A named numeric vector where names are the original variable names (or derived numeric feature names)
#'   and values are the distance correlations of each against the rest.
#' @importFrom energy bcdcor
#' @importFrom stats model.matrix
#' @importFrom pbapply pbsapply
#' @export
OneAgainstRestCorDataSummmary <- function(data) {
  if (!requireNamespace("energy", quietly = TRUE)) stop("Package 'energy' needed.")
   if (!requireNamespace("pbapply", quietly = TRUE)) stop("Package 'pbapply' needed for progress bar.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (ncol(data) < 2) stop("Need at least two columns.")

  # Convert data to numeric matrix
   tryCatch({
      data_matrix <- stats::model.matrix(~ . - 1, data = data)
      matrix_colnames <- colnames(data_matrix)
  }, error = function(e) {
      stop("Failed to convert data to numeric matrix using model.matrix: ", e$message)
  })

  n_matrix_vars <- ncol(data_matrix)
  if (n_matrix_vars < 2) stop("Need at least two numeric features after conversion.")

  cat("Calculating distance correlation of each feature against the rest for", n_matrix_vars, "features...\n")

  Corvec <- pbapply::pbsapply(1:n_matrix_vars, function(i) {
      xi <- data_matrix[, i, drop = FALSE]
      xrest <- data_matrix[, -i, drop = FALSE]
      dcor_val <- tryCatch(
          energy::bcdcor(xi, xrest),
          error = function(e) {
               warning("bcdcor failed for column '", matrix_colnames[i], "' against rest: ", e$message)
               return(NA)
          }
      )
      return(dcor_val)
  })

  names(Corvec) <- matrix_colnames
  return(Corvec)
}



#' @title SummaryTableDataSummmary
#'
#' @description Create a summary table for selected variables using the 'gtsummary' package.
#'
#' @param data Data frame.
#' @param UseVars Character vector of variable names to include in the summary.
#' @return A gtsummary table object.
#' @examples
#' \dontrun{
#' if (requireNamespace("gtsummary", quietly = TRUE)) {
#'   my_summary <- SummaryTableDataSummmary(iris, UseVars = c("Sepal.Length", "Species"))
#'   print(my_summary)
#' }
#' }
#' @importFrom gtsummary tbl_summary modify_header bold_labels italicize_labels
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @export
SummaryTableDataSummmary <- function(data, UseVars) {
    if (!requireNamespace("gtsummary", quietly = TRUE)) {
        stop("Package 'gtsummary' needed for this function. Please install it.", call. = FALSE)
    }
     if (!is.data.frame(data)) stop("'data' must be a data frame.")
     missing_vars <- setdiff(UseVars, colnames(data))
     if (length(missing_vars) > 0) {
         stop("Variables specified in 'UseVars' not found in data: ", paste(missing_vars, collapse=", "))
     }

    # Select data and create summary table using dplyr::select with all_of()
    # Ensure data is a data.frame for dplyr compatibility if it's a data.table
    summary_data <- as.data.frame(data)
    table_out <- summary_data %>%
      dplyr::select(dplyr::all_of(UseVars)) %>%
      gtsummary::tbl_summary() %>%
      gtsummary::modify_header(label = "**Variable**") %>% # Update the column header
      gtsummary::bold_labels() %>%
      gtsummary::italicize_labels() # Italicize variable labels

    return(table_out)
  }
