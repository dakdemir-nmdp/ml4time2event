tryCatch(
    {
        devtools::document()
        devtools::install(upgrade = "never")
        message("Package installed successfully.")
    },
    error = function(e) {
        stop("Installation failed: ", e$message)
    }
)

# Run the integration test
tryCatch(
    {
        testthat::test_file("tests/testthat/test-integration-full.R")
    },
    error = function(e) {
        message("Test execution failed: ", e$message)
    }
)
