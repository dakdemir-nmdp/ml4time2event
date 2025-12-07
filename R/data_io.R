#' Read Data from Various File Formats
#'
#' Reads data from CSV, XLSX, or SAS formats and converts column names to lowercase.
#'
#' @param file Character. Path to the data file to read.
#' @param filetype Character. File type - "csv", "xlsx", or "sas7bdat". Default is "csv".
#' @param ... Additional arguments passed to the underlying read function.
#'
#' @return A data frame with column names converted to lowercase.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Read CSV file
#'   data <- readData("mydata.csv")
#'   # Read XLSX file
#'   data <- readData("mydata.xlsx", filetype = "xlsx")
#' }
readData <- function(file = NULL, filetype="csv",...) {
  if (is.null(file)) {
    stop("Provide file location using the 'file' argument.")
  }
   if (!file.exists(file)) {
      stop("Specified data file does not exist: ", file)
  }

  DATA <- NULL
  tryCatch({
      if (tolower(filetype) == "csv"){
        DATA <- as.data.frame(utils::read.csv(file, stringsAsFactors = FALSE, ...)) # Avoid factors initially
      } else if (tolower(filetype) == "xlsx") {
        if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' needed for xlsx files.")
        DATA <- as.data.frame(openxlsx::read.xlsx(file, ...))
      } else if (tolower(filetype) == "sas7bdat"){
         if (!requireNamespace("haven", quietly = TRUE)) stop("Package 'haven' needed for sas7bdat files.")
        DATA <- as.data.frame(haven::read_sas(file, ...))
      } else {
        stop("Invalid filetype provided. Choose one of 'csv', 'xlsx', 'sas7bdat'.")
      }
  }, error = function(e) {
      stop("Error reading data file '", file, "' with filetype '", filetype, "': ", e$message)
  })

  # Convert column names to lowercase
  colnames(DATA) <- tolower(colnames(DATA))
  return(DATA)
}

# Consider adding writeData function here as well if needed.
