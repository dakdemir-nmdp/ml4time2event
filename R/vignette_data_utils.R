## Helper utilities to provide prepared datasets used in vignettes ---------

.vignette_extdata_path <- function(filename) {
  candidate_paths <- c(
    system.file("extdata", filename, package = "ml4time2event"),
    file.path("inst", "extdata", filename),
    file.path("..", "inst", "extdata", filename)
  )

  for (path in candidate_paths) {
    if (!identical(path, "") && file.exists(path)) {
      return(path)
    }
  }

  stop(
    "Unable to locate extdata file '", filename,
    "'. Consider ensuring the package is installed or run from the project root."
  )
}

.prepare_lung_survival_dataset <- function() {
  lung_path <- .vignette_extdata_path("lung_dataset.csv")
  lung_data <- readData(lung_path)

  if ("inst" %in% names(lung_data)) {
    lung_data <- dplyr::select(lung_data, -inst)
  }

  lung_data <- lung_data |>
    dplyr::mutate(
      sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
      ph.ecog = factor(ph.ecog)
    )

  lung_data <- ImputeMissinRecordsData(
    data = lung_data,
    dontuse = c("time", "status"),
    num.trees = 100
  )

  lung_data |>
    dplyr::mutate(
      age_group = cut(
        age,
        breaks = c(0, 50, 65, Inf),
        labels = c("Young", "Middle", "Elderly"),
        right = FALSE
      ),
      performance_good = as.integer(ph.karno >= 80)
    )
}

.prepare_bmt_competing_risks_dataset <- function() {
  bmt_path <- .vignette_extdata_path("bmtcrr_competing_risks.csv")
  bmt_data <- readData(bmt_path)

  bmt_data |>
    dplyr::mutate(
      sex = factor(sex),
      d = factor(d),
      phase = factor(phase),
      source = factor(source),
      status = as.integer(status),
      ftime = as.numeric(ftime)
    )
}

#' Retrieve the prepared lung cancer survival dataset
#'
#' @return A data frame containing the prepared NCCTG lung cancer data used in
#'   the survival vignette.
#' @examples
#' lung_df <- get_lung_survival_data()
#' @export
get_lung_survival_data <- function() {
  data_env <- new.env(parent = emptyenv())
  has_data <- tryCatch({
    utils::data(
      "lung_survival_prepared",
      package = "ml4time2event",
      envir = data_env
    )
    exists("lung_survival_prepared", envir = data_env, inherits = FALSE)
  }, warning = function(w) FALSE, error = function(e) FALSE)

  if (isTRUE(has_data)) {
    return(get("lung_survival_prepared", envir = data_env, inherits = FALSE))
  }

  .prepare_lung_survival_dataset()
}

#' Retrieve the prepared bone marrow transplant competing-risks dataset
#'
#' @return A data frame containing the prepared BMT competing-risks data used in
#'   the competing-risks vignette.
#' @examples
#' bmt_df <- get_bmt_competing_risks_data()
#' @export
get_bmt_competing_risks_data <- function() {
  data_env <- new.env(parent = emptyenv())
  has_data <- tryCatch({
    utils::data(
      "bmt_competing_risks_prepared",
      package = "ml4time2event",
      envir = data_env
    )
    exists("bmt_competing_risks_prepared", envir = data_env, inherits = FALSE)
  }, warning = function(w) FALSE, error = function(e) FALSE)

  if (isTRUE(has_data)) {
    return(get("bmt_competing_risks_prepared", envir = data_env, inherits = FALSE))
  }

  .prepare_bmt_competing_risks_dataset()
}
