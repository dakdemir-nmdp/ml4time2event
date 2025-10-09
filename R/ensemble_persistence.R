# ==============================================================================
# Ensemble Model Persistence - Save and Load Functions
# ==============================================================================

#' @title Create SurvEnsemble S3 Class
#' @description Converts ensemble model output to S3 class for better organization
#' @param ensemble_list List returned by RunSurvModels
#' @return Object of class "SurvEnsemble"
#' @export
SurvEnsemble <- function(ensemble_list) {
  # Validate input
  if (!is.list(ensemble_list)) {
    stop("ensemble_list must be a list")
  }

  if (!"input" %in% names(ensemble_list)) {
    stop("ensemble_list must contain 'input' element")
  }

  # Add class
  class(ensemble_list) <- c("SurvEnsemble", "list")
  ensemble_list
}

#' @title Create CREnsemble S3 Class
#' @description Converts ensemble model output to S3 class for better organization
#' @param ensemble_list List returned by RunCRModels
#' @return Object of class "CREnsemble"
#' @export
CREnsemble <- function(ensemble_list) {
  # Validate input
  if (!is.list(ensemble_list)) {
    stop("ensemble_list must be a list")
  }

  if (!"input" %in% names(ensemble_list)) {
    stop("ensemble_list must contain 'input' element")
  }

  # Add class
  class(ensemble_list) <- c("CREnsemble", "list")
  ensemble_list
}

#' @title Check if Object is SurvEnsemble
#' @param x Object to check
#' @return Logical
#' @export
is.SurvEnsemble <- function(x) {
  inherits(x, "SurvEnsemble")
}

#' @title Check if Object is CREnsemble
#' @param x Object to check
#' @return Logical
#' @export
is.CREnsemble <- function(x) {
  inherits(x, "CREnsemble")
}

#' @title Print SurvEnsemble Object
#' @param x SurvEnsemble object
#' @param ... Additional arguments (ignored)
#' @export
print.SurvEnsemble <- function(x, ...) {
  cat("Survival Ensemble Model\n")
  cat("=======================\n\n")

  # Show input parameters
  cat("Variables:\n")
  cat("  Explanatory:", length(x$input$ExpVars), "variables\n")
  cat("  Selected:", length(x$input$ExpVars2), "variables\n")
  cat("  Time variable:", x$input$timevar, "\n")
  cat("  Event variable:", x$input$eventvar, "\n\n")

  # Show model status
  if (!is.null(x$model_status)) {
    cat("Models fitted:\n")
    successful <- names(x$model_status)[x$model_status]
    failed <- names(x$model_status)[!x$model_status]

    cat("  Successful (", length(successful), "):", paste(successful, collapse=", "), "\n")
    if (length(failed) > 0) {
      cat("  Failed (", length(failed), "):", paste(failed, collapse=", "), "\n")
    }
  } else {
    # Count non-NULL models
    model_names <- setdiff(names(x), c("input", "model_status"))
    n_models <- sum(sapply(model_names, function(m) !is.null(x[[m]])))
    cat("Models fitted:", n_models, "\n")
  }

  invisible(x)
}

#' @title Print CREnsemble Object
#' @param x CREnsemble object
#' @param ... Additional arguments (ignored)
#' @export
print.CREnsemble <- function(x, ...) {
  cat("Competing Risks Ensemble Model\n")
  cat("==============================\n\n")

  # Show input parameters
  cat("Variables:\n")
  cat("  Explanatory:", length(x$input$ExpVars), "variables\n")
  cat("  Selected:", length(x$input$ExpVars2), "variables\n")
  cat("  Time variable:", x$input$timevar, "\n")
  cat("  Event variable:", x$input$eventvar, "\n\n")

  # Show model status
  if (!is.null(x$model_status)) {
    cat("Models fitted:\n")
    successful <- names(x$model_status)[x$model_status]
    failed <- names(x$model_status)[!x$model_status]

    cat("  Successful (", length(successful), "):", paste(successful, collapse=", "), "\n")
    if (length(failed) > 0) {
      cat("  Failed (", length(failed), "):", paste(failed, collapse=", "), "\n")
    }
  } else {
    # Count non-NULL models
    model_names <- setdiff(names(x), c("input", "model_status"))
    n_models <- sum(sapply(model_names, function(m) !is.null(x[[m]])))
    cat("Models fitted:", n_models, "\n")
  }

  invisible(x)
}

#' @title Save Ensemble Model to File
#' @description Save a survival or competing risks ensemble model to disk
#' @param ensemble_model Ensemble model object (SurvEnsemble or CREnsemble)
#' @param file Path where the model should be saved (will use .rds extension)
#' @param compress Logical or character string specifying compression (default: TRUE)
#' @return Invisibly returns the file path
#' @export
SaveEnsemble <- function(ensemble_model, file, compress = TRUE) {
  # Validate input
  if (!is.SurvEnsemble(ensemble_model) && !is.CREnsemble(ensemble_model)) {
    # Try to convert if it's a plain list
    if (is.list(ensemble_model)) {
      if ("input" %in% names(ensemble_model)) {
        message("Converting plain list to ensemble object")
        # Guess the type based on structure
        if (any(grepl("CIF|FG_Model|Fine", names(ensemble_model), ignore.case = TRUE))) {
          ensemble_model <- CREnsemble(ensemble_model)
        } else {
          ensemble_model <- SurvEnsemble(ensemble_model)
        }
      } else {
        stop("ensemble_model must be a SurvEnsemble, CREnsemble, or list with 'input' element")
      }
    } else {
      stop("ensemble_model must be a SurvEnsemble or CREnsemble object")
    }
  }

  # Add .rds extension if not present
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".rds")
  }

  # Add metadata
  ensemble_model$.__metadata__ <- list(
    saved_date = Sys.time(),
    r_version = R.version.string,
    package_version = utils::packageVersion("ml4time2event")
  )

  # Save
  saveRDS(ensemble_model, file = file, compress = compress)

  message(sprintf("Ensemble model saved to: %s", file))
  invisible(file)
}

#' @title Load Ensemble Model from File
#' @description Load a survival or competing risks ensemble model from disk
#' @param file Path to the saved ensemble model (.rds file)
#' @return Ensemble model object (SurvEnsemble or CREnsemble)
#' @export
LoadEnsemble <- function(file) {
  # Check file exists
  if (!file.exists(file)) {
    stop(sprintf("File not found: %s", file))
  }

  # Load the object
  ensemble_model <- readRDS(file)

  # Validate it's an ensemble
  if (!is.SurvEnsemble(ensemble_model) && !is.CREnsemble(ensemble_model)) {
    # Try to detect and convert
    if (is.list(ensemble_model) && "input" %in% names(ensemble_model)) {
      message("Loaded model is a plain list, attempting to convert to ensemble class")
      # Guess the type
      if (any(grepl("CIF|FG_Model|Fine", names(ensemble_model), ignore.case = TRUE))) {
        ensemble_model <- CREnsemble(ensemble_model)
      } else {
        ensemble_model <- SurvEnsemble(ensemble_model)
      }
    } else {
      warning("Loaded object may not be a valid ensemble model")
    }
  }

  # Show metadata if available
  if (!is.null(ensemble_model$.__metadata__)) {
    meta <- ensemble_model$.__metadata__
    message(sprintf("Loaded ensemble saved on %s", meta$saved_date))
    if (!is.null(meta$package_version)) {
      current_version <- utils::packageVersion("ml4time2event")
      if (meta$package_version != current_version) {
        warning(sprintf("Model was saved with package version %s, currently using %s",
                        meta$package_version, current_version))
      }
    }
  }

  ensemble_model
}

#' @title Summary Method for SurvEnsemble
#' @param object SurvEnsemble object
#' @param ... Additional arguments (ignored)
#' @export
summary.SurvEnsemble <- function(object, ...) {
  cat("Survival Ensemble Model Summary\n")
  cat("===============================\n\n")

  # Input summary
  cat("Input Configuration:\n")
  cat("  Total variables:", length(object$input$ExpVars), "\n")
  cat("  Selected variables:", length(object$input$ExpVars2), "\n")
  cat("  Time variable:", object$input$timevar, "\n")
  cat("  Event variable:", object$input$eventvar, "\n\n")

  # Model status
  if (!is.null(object$model_status)) {
    cat("Model Status:\n")
    status_df <- data.frame(
      Model = names(object$model_status),
      Status = ifelse(object$model_status, "Success", "Failed"),
      stringsAsFactors = FALSE
    )
    print(status_df, row.names = FALSE)
  }

  # Metadata if available
  if (!is.null(object$.__metadata__)) {
    cat("\nMetadata:\n")
    cat("  Saved:", format(object$.__metadata__$saved_date), "\n")
    cat("  Package version:", as.character(object$.__metadata__$package_version), "\n")
  }

  invisible(object)
}

#' @title Summary Method for CREnsemble
#' @param object CREnsemble object
#' @param ... Additional arguments (ignored)
#' @export
summary.CREnsemble <- function(object, ...) {
  cat("Competing Risks Ensemble Model Summary\n")
  cat("======================================\n\n")

  # Input summary
  cat("Input Configuration:\n")
  cat("  Total variables:", length(object$input$ExpVars), "\n")
  cat("  Selected variables:", length(object$input$ExpVars2), "\n")
  cat("  Time variable:", object$input$timevar, "\n")
  cat("  Event variable:", object$input$eventvar, "\n\n")

  # Model status
  if (!is.null(object$model_status)) {
    cat("Model Status:\n")
    status_df <- data.frame(
      Model = names(object$model_status),
      Status = ifelse(object$model_status, "Success", "Failed"),
      stringsAsFactors = FALSE
    )
    print(status_df, row.names = FALSE)
  }

  # Metadata if available
  if (!is.null(object$.__metadata__)) {
    cat("\nMetadata:\n")
    cat("  Saved:", format(object$.__metadata__$saved_date), "\n")
    cat("  Package version:", as.character(object$.__metadata__$package_version), "\n")
  }

  invisible(object)
}
