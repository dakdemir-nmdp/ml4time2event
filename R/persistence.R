#' Save a fitted model or pipeline
#'
#' Persists a `t2e_fit` or `ml4t2e_pipeline` object to disk, preserving
#' all metadata, trained models, and preprocessing recipes.
#'
#' @param object A `t2e_fit` or `ml4t2e_pipeline` object.
#' @param file Path to save the object.
#' @param compress Logical or character string specifying compression. Default is `TRUE`.
#' @return Invisible `NULL`.
#' @export
ml4t2e_save <- function(object, file, compress = TRUE) {
    if (!inherits(object, "t2e_fit") && !inherits(object, "ml4t2e_pipeline")) {
        rlang::abort("`object` must be a `t2e_fit` or `ml4t2e_pipeline`.")
    }

    model_data <- list(
        object = object,
        metadata = list(
            saved_date = Sys.time(),
            r_version = R.version.string,
            package_version = utils::packageVersion("ml4time2event")
        )
    )

    saveRDS(model_data, file = file, compress = compress)
    invisible(NULL)
}

#' Load a fitted model or pipeline
#'
#' Loads a previously saved `ml4time2event` object.
#'
#' @param file Path to the saved object.
#' @return The loaded `t2e_fit` or `ml4t2e_pipeline` object.
#' @export
ml4t2e_load <- function(file) {
    if (!file.exists(file)) {
        rlang::abort(glue::glue("File not found: {file}"))
    }

    model_data <- readRDS(file)

    # Check if it was saved with our wrapper
    if (is.list(model_data) && "metadata" %in% names(model_data) && "object" %in% names(model_data)) {
        meta <- model_data$metadata

        current_version <- utils::packageVersion("ml4time2event")
        if (!is.null(meta$package_version) && meta$package_version != current_version) {
            rlang::warn(glue::glue(
                "Object was saved with package version {meta$package_version}, ",
                "but you are using version {current_version}. Behavior may differ."
            ))
        }

        return(model_data$object)
    }

    # Fallback: maybe it was saved with plain saveRDS
    if (inherits(model_data, "t2e_fit") || inherits(model_data, "ml4t2e_pipeline")) {
        return(model_data)
    }

    rlang::abort("Loaded object is not a valid `ml4time2event` model or pipeline.")
}
