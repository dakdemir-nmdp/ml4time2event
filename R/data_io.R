#' @title readData
#'
#' @description Import data from csv, xlsx, or sas7bdat file.
#' Converts column names to lowercase.
#'
#' @param file Path to the data file.
#' @param filetype Type of data file ("csv", "xlsx", "sas7bdat").
#' @param ... Additional arguments passed to the respective reading function
#'   (read.csv, openxlsx::read.xlsx, haven::read_sas).
#' @return A data frame containing the imported data with lowercase column names.
#' @examples
#' \dontrun{
#' mydata <- readData(file = "path/to/data.csv", filetype = "csv")
#' }
#' @importFrom utils read.csv
#' @importFrom openxlsx read.xlsx
#' @importFrom haven read_sas
#' @export
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
