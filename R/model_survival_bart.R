#' Bayesian Additive Regression Trees Survival Model (R6)
#'
#' Fits a Survival BART model.
#'
#' @keywords internal
#' @noRd
BartSurvivalModel <- R6::R6Class(
  classname = "BartSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    factor_levels = NULL,
    x_train = NULL,
    times_train = NULL,
    delta_train = NULL,
    varprof = NULL,
    feature_names = NULL,
    bart_feature_names = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      self$varprof <- VariableProfile(data, task$features)

      # Params
      spec <- self$spec
      K <- spec$K %||% 10
      ntree <- spec$ntree %||% 50
      ndpost <- spec$ndpost %||% 200
      nskip <- spec$nskip %||% 50
      keepevery <- spec$keepevery %||% 10L
      mc.cores <- spec$mc.cores %||% 1L
      verbose <- spec$verbose %||% FALSE

      # Data Prep
      expvars <- task$features
      timevar <- task$time_col
      eventvar <- task$event_col

      # Store factor levels
      self$factor_levels <- lapply(data[, expvars, drop = FALSE], function(x) {
        if (is.factor(x)) levels(x) else NULL
      })

      # Drop rows with NA features BEFORE model.matrix so that times_vec /
      # x_train_mat stay in sync.  (model.matrix uses na.action=na.omit by
      # default and silently shortens the matrix, which makes
      # length(times) != nrow(x.train) and breaks surv.pre.bart.)
      complete_rows <- stats::complete.cases(data[, expvars, drop = FALSE])
      if (!all(complete_rows)) {
        rlang::warn(sprintf(
          "BART fit: dropping %d rows with missing feature values.",
          sum(!complete_rows)
        ))
        data <- data[complete_rows, , drop = FALSE]
      }

      times_vec <- data[[timevar]]
      delta_vec <- as.integer(data[[eventvar]] == 1)

      # Model Matrix
      # BART requires matrix
      x_train_mat <- as.matrix(stats::model.matrix(~ -1 + ., data = data[, expvars, drop = FALSE]))
      self$feature_names <- colnames(x_train_mat)

      # Capture which columns survive bartModelMatrix (rm.const / rm.dup) so
      # we can align the test matrix to the exact same columns during predict.
      # surv.bart calls bartModelMatrix internally; if we don't mirror that here
      # the test matrix ends up with a different column count and predict fails.
      x_train_bart <- BART::bartModelMatrix(x_train_mat)
      self$bart_feature_names <- colnames(x_train_bart)

      # Fit
      bart_fit <- NULL
      if (verbose) {
        bart_fit <- suppressMessages(BART::mc.surv.bart(
          mc.cores = mc.cores,
          x.train = x_train_mat,
          times = times_vec,
          delta = delta_vec,
          x.test = x_train_mat, # minimal dummy
          K = K,
          ntree = ntree,
          ndpost = ndpost,
          nskip = nskip,
          keepevery = keepevery
        ))
      } else {
        invisible(capture.output({
          bart_fit <- suppressMessages(BART::mc.surv.bart(
            mc.cores = mc.cores,
            x.train = x_train_mat,
            times = times_vec,
            delta = delta_vec,
            x.test = x_train_mat,
            K = K,
            ntree = ntree,
            ndpost = ndpost,
            nskip = nskip,
            keepevery = keepevery
          ))
        }))
      }

      self$model <- bart_fit
      self$x_train <- x_train_mat
      self$times_train <- times_vec
      self$delta_train <- delta_vec

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()

      # Prepare newdata
      if (!self$task$id_col %in% colnames(newdata)) {
        # Warning or error handled by .ensure?
        # .ensure_prediction_data handles missing ID
      }

      # We need original columns for factor level enforcement
      newdata_df <- newdata

      # Enforce factor levels
      for (vari in self$task$features) {
        levels_train <- self$factor_levels[[vari]]
        if (!is.null(levels_train)) {
          if (vari %in% colnames(newdata_df)) {
            newdata_df[[vari]] <- factor(as.character(newdata_df[[vari]]), levels = levels_train)
          }
        }
      }

      # Model Matrix
      # Fast-path for numeric-only features to avoid slow model.matrix calls at inference
      feat_data <- newdata_df[, self$task$features, drop = FALSE]
      if (all(vapply(feat_data, is.numeric, FUN.VALUE = logical(1)))) {
        x_test_mat <- as.matrix(feat_data)
      } else {
        x_test_mat <- stats::model.matrix(~ -1 + ., data = feat_data)
      }

      # Align test columns to the exact training model.matrix columns used in
      # fit().  This is required so surv.pre.bart() can build tx.test with the
      # same structure expected by predict.survbart().
      missing_train_cols <- setdiff(self$feature_names, colnames(x_test_mat))
      if (length(missing_train_cols) > 0) {
        add_mat <- matrix(0,
          nrow = nrow(x_test_mat), ncol = length(missing_train_cols),
          dimnames = list(NULL, missing_train_cols)
        )
        x_test_mat <- cbind(x_test_mat, add_mat)
      }
      x_test_mat <- x_test_mat[, self$feature_names, drop = FALSE]

      # Prepare BART prediction matrix. x.train is required here; omitting it
      # can produce a malformed tx.test and trigger a column-mismatch error in
      # predict.survbart().  Pass events= to reuse the fitted model time grid.
      pre <- BART::surv.pre.bart(
        times = self$times_train,
        delta = self$delta_train,
        x.train = self$x_train,
        x.test = x_test_mat,
        events = self$model$times
      )

      # Predict using optimized C++ multi-threaded backend
      mc.cores <- self$spec$mc.cores %||% 1L
      pred <- BART::mc.surv.pwbart(
        x.test = pre$tx.test,
        treedraws = self$model$treedraws,
        binaryOffset = self$model$binaryOffset,
        mc.cores = mc.cores
      )

      # Extract mean survival
      # pred$surv.test.mean is vector [patient1_times, patient2_times...]
      n_times_bart <- self$model$K
      n_obs <- nrow(x_test_mat)

      # Matrix: rows=obs, cols=times
      pred_mat <- matrix(pred$surv.test.mean, ncol = n_times_bart, byrow = TRUE)

      # Sort times from fit
      time_order <- order(self$model$times)
      sorted_times <- self$model$times[time_order]
      sorted_probs <- pred_mat[, time_order, drop = FALSE]

      # Add t=0, S=1
      probs_with_0 <- cbind(1, sorted_probs)
      times_with_0 <- c(0, sorted_times)

      # Interpolate
      req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Interpolator needs (rows=times, cols=obs)
      probs_t <- t(probs_with_0)

      interp_probs <- survprobMatInterpolator(probs_t, times_with_0, req_times)

      # Flatten
      id_col <- newdata[[self$task$id_col]]
      # Use generated IDs if missing from newdata logic?
      if (is.null(id_col)) id_col <- seq_len(n_obs)

      new_survival_prediction(
        id = rep(id_col, each = length(req_times)),
        time = rep(req_times, times = length(id_col)),
        surv = as.vector(interp_probs),
        model = rep("bart", length(id_col) * length(req_times)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "BART Survival"
      info
    },
    required_packages = function() {
      c("BART")
    }
  ),
  private = list(
    ensure_fitted = function() {
      if (!isTRUE(self$fitted)) {
        rlang::abort("Model must be fitted before predictions can be generated.")
      }
    }
  )
)

# Register
.register_time_to_event_model(
  engine = "bart",
  outcome = "survival",
  constructor = function(spec = list()) {
    BartSurvivalModel$new(spec = modifyList(list(engine = "bart"), spec))
  },
  packages = "BART",
  tags = c("bart", "bayesian"),
  label = "BART survival"
)
