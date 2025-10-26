# This script prepares cleaned datasets used in the package vignettes
# and stores them in the data/ directory for easy access via data().

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop("The 'usethis' package is required to store vignette datasets.")
}

suppressPackageStartupMessages({
  library(dplyr)
})

# Ensure we can source package functions when running from the package root
package_root <- normalizePath("..", mustWork = FALSE)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  package_root <- normalizePath(".", mustWork = TRUE)
}

source(file.path(package_root, "R", "data_io.R"))
source(file.path(package_root, "R", "data_cleaning.R"))
source(file.path(package_root, "R", "vignette_data_utils.R"))

# --- Prepare Lung Survival Dataset ---
lung_survival_prepared <- .prepare_lung_survival_dataset()

usethis::use_data(lung_survival_prepared, overwrite = TRUE)

# --- Prepare BMT Competing Risks Dataset ---
bmt_competing_risks_prepared <- .prepare_bmt_competing_risks_dataset()

usethis::use_data(bmt_competing_risks_prepared, overwrite = TRUE)

message("Prepared vignette datasets saved to data/.")
