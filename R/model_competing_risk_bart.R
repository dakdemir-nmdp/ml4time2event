#' BART Competing Risk Model (R6)
#'
#' Fits a cause-specific BART model for each competing risk (using One vs Rest-Competing strategy).
#'
#' @keywords internal
#' @noRd
BartCompetingRiskModel <- R6::R6Class(
    classname = "BartCompetingRiskModel",
    inherit = CompetingRiskModel,
    public = list(
        model = NULL, # list of models per cause
        time_grid = NULL,
        task = NULL,
        cause_codes = NULL,
        factor_levels = NULL,
        x_train = NULL, # Need to store for prediction (BART Requirement)
        feature_names = NULL,

        # Store training metadata per cause for prediction
        train_metadata = NULL,
        fit = function(task, time_grid, ...) {
            super$fit(task = task, ...)
            self$task <- task
            self$time_grid <- time_grid

            data <- as.data.frame(task$data)
            cause_map <- task$metadata$cause_map
            self$cause_codes <- setNames(cause_map$code, as.character(cause_map$cause))

            spec <- self$spec
            K <- spec$K %||% 10
            ntree <- spec$ntree %||% 50
            ndpost <- spec$ndpost %||% 200
            nskip <- spec$nskip %||% 50
            keepevery <- spec$keepevery %||% 10L
            numcut <- spec$numcut %||% 100
            mc.cores <- spec$mc.cores %||% 1L
            verbose <- spec$verbose %||% FALSE

            expvars <- task$features
            timevar <- task$time_col
            eventcodevar <- task$event_col # numeric codes

            # Drop rows with NA features before model.matrix so y vectors and
            # x_train stay synchronized (model.matrix otherwise omits rows).
            complete_rows <- stats::complete.cases(data[, expvars, drop = FALSE])
            if (!all(complete_rows)) {
                rlang::warn(sprintf(
                    "BART CR fit: dropping %d rows with missing feature values.",
                    sum(!complete_rows)
                ))
                data <- data[complete_rows, , drop = FALSE]
            }

            # Factor levels
            self$factor_levels <- lapply(data[, expvars, drop = FALSE], function(x) {
                if (is.factor(x)) levels(x) else NULL
            })

            # Prepare Matrix
            x_train_mat <- as.matrix(stats::model.matrix(~ -1 + ., data = data[, expvars, drop = FALSE]))
            self$x_train <- x_train_mat
            self$feature_names <- colnames(x_train_mat)

            y_times <- data[[timevar]]
            y_codes <- data[[eventcodevar]]

            models_list <- list()
            train_meta_list <- list()
            cause_labels <- names(self$cause_codes)

            for (cause_label in cause_labels) {
                cause_code <- self$cause_codes[[cause_label]]

                # Recode delta for this cause:
                # 1 = cause of interest
                # 2 = competing risk
                # 0 = censored

                delta_k <- integer(length(y_codes))
                delta_k[y_codes == 0] <- 0L
                delta_k[y_codes == cause_code] <- 1L
                delta_k[y_codes != 0 & y_codes != cause_code] <- 2L

                # If no competing risks (only 0 and 1)? crisk.bart might fail or reduce to surv.bart?
                # If no competing events observed, delta has 0 and 1. 2 is absent.
                # BART might complain.
                # But for CR task, improved logic:
                # If sum(delta==2) == 0, use surv.bart?
                # But we want CR interface.
                # Assuming proper CR data.

                if (sum(delta_k == 1) < 1) {
                    rlang::warn(glue::glue("No events for cause {cause_label}. Skipping."))
                    next
                }

                # Fit
                if (verbose) {
                    bart_fit <- suppressMessages(BART::mc.crisk.bart(
                        mc.cores = mc.cores,
                        x.train = x_train_mat,
                        times = y_times,
                        delta = delta_k,
                        x.test = x_train_mat,
                        K = K,
                        ntree = ntree,
                        ndpost = ndpost,
                        nskip = nskip,
                        keepevery = keepevery,
                        numcut = numcut
                    ))
                } else {
                    bart_fit <- NULL
                    invisible(capture.output({
                        bart_fit <- suppressMessages(BART::mc.crisk.bart(
                            mc.cores = mc.cores,
                            x.train = x_train_mat,
                            times = y_times,
                            delta = delta_k,
                            x.test = x_train_mat,
                            K = K,
                            ntree = ntree,
                            ndpost = ndpost,
                            nskip = nskip,
                            keepevery = keepevery,
                            numcut = numcut
                        ))
                    }))
                }

                models_list[[cause_label]] <- bart_fit
                train_meta_list[[cause_label]] <- list(
                    times_train = y_times,
                    delta_train = delta_k
                )
            }

            self$model <- models_list
            self$train_metadata <- train_meta_list
            invisible(self)
        },
        predict_cif = function(newdata, times, set = "test", ...) {
            private$ensure_fitted()

            # Prepare newdata (factor levels)
            newdata_df <- newdata
            for (vari in self$task$features) {
                levels_train <- self$factor_levels[[vari]]
                if (!is.null(levels_train)) {
                    if (vari %in% colnames(newdata_df)) {
                        newdata_df[[vari]] <- factor(as.character(newdata_df[[vari]]), levels = levels_train)
                    }
                }
            }

            # Matrix fast path for numerical data
            feat_data <- newdata_df[, self$task$features, drop = FALSE]
            if (all(vapply(feat_data, is.numeric, FUN.VALUE = logical(1)))) {
                x_test_mat <- as.matrix(feat_data)
            } else {
                x_test_mat <- stats::model.matrix(~ -1 + ., data = feat_data)
            }

            # Align test columns to the exact training model.matrix columns.
            # This prevents fold-specific dummy level differences from causing
            # mismatches when generating BART prediction design matrices.
            missing_train_cols <- setdiff(self$feature_names, colnames(x_test_mat))
            if (length(missing_train_cols) > 0) {
                add_mat <- matrix(
                    0,
                    nrow = nrow(x_test_mat),
                    ncol = length(missing_train_cols),
                    dimnames = list(NULL, missing_train_cols)
                )
                x_test_mat <- cbind(x_test_mat, add_mat)
            }
            x_test_mat <- x_test_mat[, self$feature_names, drop = FALSE]

            preds_list <- list()
            cause_labels <- names(self$cause_codes)
            id_complete <- newdata[[self$task$id_col]]
            if (is.null(id_complete)) id_complete <- seq_len(nrow(x_test_mat))

            req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))
            n_obs <- nrow(x_test_mat)

            for (cause in cause_labels) {
                bart_model <- self$model[[cause]]
                train_meta <- self$train_metadata[[cause]]

                if (is.null(bart_model)) {
                    # NA predictions
                    cif_vec <- rep(NA_real_, n_obs * length(req_times))
                } else {
                    # Pre
                    pre <- BART::crisk.pre.bart(
                        time = train_meta$times_train,
                        delta = train_meta$delta_train,
                        x.train = self$x_train,
                        x.test = x_test_mat,
                        x.train2 = self$x_train,
                        x.test2 = x_test_mat,
                        K = bart_model$K
                    )

                    tx_test <- pre$tx.test
                    tx_test2 <- pre$tx.test2

                    # crisk.bart may drop constant columns during fit; apply the
                    # same column selections to prediction matrices.
                    if (!is.null(bart_model$rm.const)) {
                        tx_test <- tx_test[, bart_model$rm.const, drop = FALSE]
                    }
                    if (!is.null(bart_model$rm.const2)) {
                        tx_test2 <- tx_test2[, bart_model$rm.const2, drop = FALSE]
                    }

                    # Predict
                    mc.cores <- self$spec$mc.cores %||% 1L
                    bOff1 <- if (is.null(bart_model$binaryOffset)) 0 else bart_model$binaryOffset
                    bOff2 <- if (is.null(bart_model$binaryOffset2)) 0 else bart_model$binaryOffset2

                    pred <- BART::mc.crisk.pwbart(
                        x.test = tx_test,
                        x.test2 = tx_test2,
                        treedraws = bart_model$treedraws,
                        treedraws2 = bart_model$treedraws2,
                        binaryOffset = bOff1,
                        binaryOffset2 = bOff2,
                        mc.cores = mc.cores
                    )

                    # Reshape
                    # pred$cif.test.mean [N*K]
                    n_times_bart <- bart_model$K
                    cif_matrix <- matrix(pred$cif.test.mean, nrow = n_obs, ncol = n_times_bart, byrow = TRUE)

                    # Sort times
                    time_order <- order(bart_model$times)
                    sorted_times <- bart_model$times[time_order]
                    sorted_cif <- cif_matrix[, time_order, drop = FALSE]

                    # Add t=0
                    cif_with_0 <- cbind(0, sorted_cif)
                    times_with_0 <- c(0, sorted_times)

                    # Interpolate
                    # (rows=times)
                    cif_t <- t(cif_with_0)
                    interp_cif <- cifMatInterpolator(cif_t, times_with_0, req_times)

                    cif_vec <- as.vector(interp_cif)
                }

                new_cif_prediction(
                    id = rep(id_complete, each = length(req_times)),
                    time = rep(req_times, times = length(id_complete)),
                    cause = rep(cause, length(id_complete) * length(req_times)),
                    cif = cif_vec,
                    model = rep("bart", length(id_complete) * length(req_times)),
                    ensemble = FALSE,
                    set = set
                ) -> pred_obj
                preds_list[[length(preds_list) + 1]] <- pred_obj
            }

            result <- dplyr::bind_rows(preds_list)
            result$cif <- pmin(pmax(result$cif, 0), 1)
            result
        },
        model_info = function() {
            info <- super$model_info()
            info$label <- "BART Competing Risk"
            info
        },
        required_packages = function() {
            c("BART")
        }
    ),
    private = list(
        ensure_fitted = function() {
            if (!isTRUE(self$fitted)) {
                rlang::abort("Model not fitted")
            }
        }
    )
)

.register_time_to_event_model(
    engine = "cr_bart",
    outcome = "competing_risk",
    constructor = function(spec = list()) {
        BartCompetingRiskModel$new(spec = modifyList(list(engine = "cr_bart"), spec))
    },
    packages = "BART",
    tags = c("bart", "competing-risk"),
    label = "BART Competing Risk"
)
