#' Install Extra Dependencies
#'
#' Installs the optional heavy dependencies required for specific engines
#' (e.g., xgboost, gbm, glmnet, BART) and advanced explainability features.
#'
#' @param force Logical. If TRUE, force re-installation even if already present. Default is FALSE.
#' @return Invisible NULL.
#' @export
ml4t2e_install_extras <- function(force = FALSE) {
    extras <- c(
        "BART",
        "fastcmprsk",
        "gbm",
        "glmnet",
        "mgcv",
        "nnet",
        "partykit",
        "pseudo",
        "randomForestSRC",
        "rpart",
        "xgboost",
        "fastshap",
        "kernelshap",
        "shapviz",
        "mlr3",
        "mlr3learners"
    )

    installed <- rownames(utils::installed.packages())
    to_install <- if (force) extras else setdiff(extras, installed)

    if (length(to_install) == 0) {
        message("All extra dependencies are already installed.")
        return(invisible())
    }

    message("Installing extra dependencies: ", paste(to_install, collapse = ", "))
    utils::install.packages(to_install)

    invisible()
}
