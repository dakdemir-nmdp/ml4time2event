#' RuleFit Competing Risk Model (R6)
#'
#' Fits a RuleFit model for competing risks.
#' Uses pseudo-observations to train regression trees for rule generation,
#' then fits a Fine-Gray model on the augmented feature space.
#'
#' @keywords internal
#' @noRd
RulefitCompetingRiskModel <- R6::R6Class(
    classname = "RulefitCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL,
        time_grid = NULL,
        task = NULL,
        trees_list = NULL,
        rules_list = NULL,
        fg_model = NULL, # Instance of FineGrayCompetingRiskModel
        # Need to store sample for alignment
        sample_data = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            expvars <- task$features
            timevar <- task$time_col
            statusvar <- task$event_col # This is 0/1/2.. causes

            # Params
            spec <- self$spec
            ntree <- spec$ntree %||% 30
            nsample <- spec$nsample %||% min(300, nrow(data))
            keepvars <- spec$keepvars

            cause_map <- task$metadata$cause_map
            causes <- as.character(cause_map$code)

            models_list <- list()

            for (cause_str in causes) {
                cause <- as.numeric(cause_str)

                # 1. Tree Generation
                x_time <- data[[timevar]]
                cuttimesReg <- stats::quantile(x_time, c(.25, .50, .60, .70, .80), na.rm = TRUE)

                trees <- list()
                for (i in 1:ntree) {
                    # Sampling
                    sampcols <- union(keepvars, sample(expvars, min(length(expvars), sample(seq_len(min(10, length(expvars))), 1))))
                    samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)

                    datasampl <- data[samprows, c(timevar, statusvar, sampcols), drop = FALSE]

                    # Pseudo
                    tmax <- sample(cuttimesReg, 1)
                    pout <- tryCatch(
                        {
                            pseudo::pseudoci(time = datasampl[[timevar]], event = datasampl[[statusvar]], tmax = tmax)
                        },
                        error = function(e) NULL
                    )

                    if (is.null(pout)) next

                    cause_name <- paste0("cause", cause)
                    if (!cause_name %in% names(pout$pseudo)) next

                    pseudo_val <- pout$pseudo[[cause_name]]
                    if (length(pseudo_val) != nrow(datasampl)) next

                    # Standardize
                    reg_var <- (pseudo_val - mean(pseudo_val)) / (stats::sd(pseudo_val) + 1e-10)

                    datasampl$RegVar <- reg_var
                    datasampl <- datasampl[, c("RegVar", sampcols), drop = FALSE]

                    # Fit Regression Tree
                    ctrl <- rpart::rpart.control(
                        minsplit = stats::rpois(1, 2) + 1,
                        minbucket = stats::rpois(1, 3) + 1,
                        cp = 0.01 * stats::runif(1),
                        maxdepth = stats::rpois(1, 2) + 1,
                        usesurrogate = 0, xval = 0
                    )

                    tree_fit <- tryCatch(
                        {
                            rpart::rpart(RegVar ~ ., data = datasampl, control = ctrl, model = TRUE)
                        },
                        error = function(e) NULL
                    )

                    if (!is.null(tree_fit)) trees[[length(trees) + 1]] <- tree_fit
                }

                if (length(trees) == 0) {
                    rlang::warn(glue::glue("No trees generated for cause {cause}. Skipping."))
                    next
                }

                # 2. Extract Rules
                rules_list_c <- lapply(trees, function(x) {
                    tryCatch(.extract_rules_from_party(partykit::as.party(x)), error = function(e) list())
                })

                # 3. Create Augmented Data
                rules_mats <- list()
                for (k in seq_along(trees)) {
                    tr <- trees[[k]]
                    rls <- rules_list_c[[k]]
                    if (length(rls) <= 1) next

                    ptree <- partykit::as.party(tr)
                    nd <- predict(ptree, newdata = data, type = "node")
                    tf <- factor(as.character(nd), levels = names(rls))

                    if (length(levels(tf)) > 1) {
                        mm <- tryCatch(stats::model.matrix(~ -1 + tf), error = function(e) NULL)
                        if (!is.null(mm)) {
                            colnames(mm) <- paste(paste(k, seq_along(rls), sep = "_"), rls, sep = "***")
                            rules_mats[[length(rules_mats) + 1]] <- mm
                        }
                    }
                }

                if (length(rules_mats) > 0) {
                    mm_rules <- Reduce("cbind", rules_mats)
                    train_aug <- cbind(data, mm_rules)
                } else {
                    train_aug <- data
                }

                # 4. Fit Fine-Gray
                # features
                rule_vars <- if (length(rules_mats) > 0) colnames(mm_rules) else character(0)
                fg_features <- c(expvars, rule_vars)

                fg_task <- ml4t2e_task_cr(
                    train_aug,
                    time = timevar,
                    status = task$metadata$status_col,
                    cause = task$cause_col,
                    features = fg_features
                )

                fg_model <- FineGrayCompetingRiskModel$new(spec = list(engine = "cr_fine_gray"))

                # Target only current cause
                fg_task$metadata$cause_map <- fg_task$metadata$cause_map[fg_task$metadata$cause_map$code == cause, ]

                fg_model$fit(fg_task, time_grid = self$time_grid)

                models_list[[cause_str]] <- list(
                    trees = trees,
                    rules = rules_list_c,
                    fg_model = fg_model,
                    rule_vars = rule_vars
                )
            }

            self$model <- models_list
            self$sample_data <- data[seq_len(min(nrow(data), 10)), expvars, drop = FALSE]
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()
            complete_data <- .ensure_prediction_data(newdata, self$task)
            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

            preds_list <- list()
            id_col <- complete_data[[self$task$id_col]]

            # Prepare data factors once
            data_test <- complete_data[, self$task$features, drop = FALSE]
            for (v in self$task$features) {
                # Ensure levels match if needed
            }

            for (cause_str in names(self$model)) {
                comp <- self$model[[cause_str]]
                trees <- comp$trees
                rules <- comp$rules
                fg <- comp$fg_model

                # Generate rules matrix
                rules_mats <- list()
                for (k in seq_along(trees)) {
                    tr <- trees[[k]]
                    rls <- rules[[k]]
                    if (length(rls) <= 1) next

                    ptree <- partykit::as.party(tr)
                    nd <- predict(ptree, newdata = data_test, type = "node")
                    tf <- factor(as.character(nd), levels = names(rls))

                    if (length(levels(tf)) > 1) {
                        mm <- tryCatch(stats::model.matrix(~ -1 + tf), error = function(e) NULL)
                        if (!is.null(mm)) {
                            colnames(mm) <- paste(paste(k, seq_along(rls), sep = "_"), rls, sep = "***")
                            rules_mats[[length(rules_mats) + 1]] <- mm
                        }
                    }
                }

                if (length(rules_mats) > 0) {
                    mm_rules <- Reduce("cbind", rules_mats)
                    test_aug <- cbind(complete_data, mm_rules)
                } else {
                    test_aug <- complete_data
                }

                # Predict
                fg_preds <- fg$predict_cif(test_aug, times = req_times, set = set)
                preds_list[[length(preds_list) + 1]] <- fg_preds
            }

            dplyr::bind_rows(preds_list)
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "RuleFit Competing Risk"
            info
        },
        required_packages = function() {
            c("rpart", "partykit", "glmnet", "survival", "pseudo", "fastcmprsk")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
        }
    )
)

# Helper reuse
.extract_rules_from_party <- function(x) {
    if (is.null(x) || is.null(x$node)) {
        return(list())
    }
    if (exists(".list.rules.party", where = asNamespace("partykit"), mode = "function")) {
        rule_strings <- partykit:::.list.rules.party(x)
        node_ids <- names(partykit::nodeids(x, terminal = TRUE))
        names(rule_strings) <- node_ids
        return(rule_strings)
    }
    return(list())
}

.register_time_to_event_model(
    engine = "cr_rulefit",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        RulefitCompetingRiskModel$new(spec = modifyList(list(engine = "cr_rulefit"), spec))
    },
    packages = c("rpart", "partykit", "glmnet", "survival", "pseudo", "fastcmprsk"),
    tags = c("rules", "tree-based"),
    label = "RuleFit Competing Risk"
)
