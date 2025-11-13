#' Penalised Cox model via glmnet (R6 implementation)
#'
#' Provides an elastic-net penalised Cox proportional hazards model that follows
#' the shared time-to-event interface.
#'
#' @keywords internal
#' @noRd
GlmnetSurvivalModel <- R6::R6Class(
  classname = "GlmnetSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    train_matrix = NULL,
    train_response = NULL,
    feature_levels = NULL,
    training_features = NULL,

    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      data <- task$data
      required <- c(task$time_col, task$event_col, task$features)
      data <- data[, required, drop = FALSE]
      data <- stats::na.omit(data)

      self$training_features <- data[, task$features, drop = FALSE]
      self$feature_levels <- lapply(self$training_features, function(col) {
        if (is.factor(col)) {
          levels(col)
        } else {
          NULL
        }
      })

      x <- stats::model.matrix(~ . - 1, data = self$training_features)
      y <- survival::Surv(
        time = data[[task$time_col]],
        event = data[[task$event_col]]
      )

      spec <- self$spec
      alpha <- spec$alpha
      if (is.null(alpha)) {
        alpha <- 0.5
      }
      nfolds <- spec$nfolds
      if (is.null(nfolds)) {
        nfolds <- min(10, nrow(x))
      }
      maxit <- spec$maxit
      if (is.null(maxit)) {
        maxit <- 1000L
      }

      self$model <- glmnet::cv.glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = alpha,
        nfolds = nfolds,
        maxit = maxit
      )
      self$train_matrix <- x
      self$train_response <- y
      self$time_grid <- time_grid
      self$task <- task
      invisible(self)
    },

    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      newdata <- .ensure_prediction_data(newdata, self$task)
      model_name <- self$spec$engine
      if (is.null(model_name)) {
        model_name <- "glmnet"
      }
      if (is.null(times)) {
        target_times <- self$time_grid
      } else {
        target_times <- sort(unique(as.numeric(times)))
      }
      if (length(target_times) == 0) {
        rlang::abort("`times` must be numeric and non-empty.")
      }

      complete_idx <- stats::complete.cases(newdata)
      id_values <- newdata[[self$task$id_col]]

      preds_complete <- NULL
      if (any(complete_idx)) {
        newdata_complete <- newdata[complete_idx, , drop = FALSE]
        id_complete <- newdata_complete[[self$task$id_col]]
        feature_frame <- newdata_complete[, self$task$features, drop = FALSE]
        mat <- private$make_model_matrix(feature_frame)
        sfit <- survival::survfit(
          self$model,
          s = "lambda.min",
          x = self$train_matrix,
          y = self$train_response,
          newx = mat
        )
        base_times <- sfit$time
        surv_matrix <- t(sfit$surv)
        aligned <- .glmnet_align_survival(surv_matrix, base_times, target_times)

        preds_complete <- new_survival_prediction(
          id = rep(id_complete, each = length(target_times)),
          time = rep(target_times, times = length(id_complete)),
          surv = as.vector(aligned),
          model = rep(model_name, length(id_complete) * length(target_times)),
          ensemble = FALSE,
          set = set
        )
      }

      preds_missing <- NULL
      if (!all(complete_idx)) {
        missing_ids <- id_values[!complete_idx]
        rlang::warn(glue::glue(
          "Omitting {length(missing_ids)} rows with missing predictors for engine '{model_name}'."
        ))
        preds_missing <- new_survival_prediction(
          id = rep(missing_ids, each = length(target_times)),
          time = rep(target_times, times = length(missing_ids)),
          surv = rep(NA_real_, length(missing_ids) * length(target_times)),
          model = rep(model_name, length(missing_ids) * length(target_times)),
          ensemble = FALSE,
          set = set
        )
      }

      pieces <- list(preds_complete, preds_missing)
      pieces <- pieces[!vapply(pieces, is.null, logical(1))]
      if (length(pieces) == 0) {
        return(new_survival_prediction(
          id = integer(0),
          time = numeric(0),
          surv = numeric(0),
          model = character(0),
          ensemble = logical(0),
          set = character(0)
        ))
      }
      dplyr::bind_rows(pieces)
    },

    predict_risk = function(newdata, times = NULL, set = "test", ...) {
      target_times <- times
      if (is.null(target_times)) {
        target_times <- self$time_grid
      }
      survival_tbl <- self$predict_survival(newdata = newdata, times = target_times, set = set, ...)
      if (nrow(survival_tbl) == 0) {
        return(new_risk_prediction(
          id = integer(0),
          risk = numeric(0),
          model = character(0),
          time = numeric(0),
          ensemble = logical(0),
          set = character(0)
        ))
      }
      last_time <- max(target_times)
      last_slice <- survival_tbl[survival_tbl$time == last_time, , drop = FALSE]
      new_risk_prediction(
        id = last_slice$id,
        risk = 1 - last_slice$surv,
        model = last_slice$model,
        time = last_time,
        ensemble = FALSE,
        set = set
      )
    },

    model_info = function() {
      info <- super$model_info()
      info$label <- "Penalised Cox (glmnet)"
      info
    },

    required_packages = function() {
      c("glmnet", "survival")
    }
  ),
  private = list(
    ensure_fitted = function() {
      if (!isTRUE(self$fitted)) {
        rlang::abort("Model must be fitted before predictions can be generated.")
      }
    },
    make_model_matrix = function(newdata) {
      features <- self$task$features
      for (feature in features) {
        levels <- self$feature_levels[[feature]]
        if (!is.null(levels)) {
          newdata[[feature]] <- factor(newdata[[feature]], levels = levels)
        }
      }
      combined <- rbind(self$training_features, newdata)
      combined_mat <- stats::model.matrix(~ . - 1, data = combined)
      n_train <- nrow(self$training_features)
      new_mat <- combined_mat[-seq_len(n_train), , drop = FALSE]
      train_cols <- colnames(self$train_matrix)
      missing_cols <- setdiff(train_cols, colnames(new_mat))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {
          new_mat <- cbind(new_mat, 0)
          colnames(new_mat)[ncol(new_mat)] <- col
        }
      }
      extra_cols <- setdiff(colnames(new_mat), train_cols)
      if (length(extra_cols) > 0) {
        new_mat <- new_mat[, setdiff(colnames(new_mat), extra_cols), drop = FALSE]
      }
      new_mat <- new_mat[, train_cols, drop = FALSE]
      new_mat
    }
  )
)

.glmnet_align_survival <- function(surv_matrix, base_times, target_times) {
  n_obs <- nrow(surv_matrix)
  out <- matrix(NA_real_, nrow = n_obs, ncol = length(target_times))
  for (i in seq_len(n_obs)) {
    surv_vec <- surv_matrix[i, ]
    interpolator <- stats::approxfun(
      x = base_times,
      y = surv_vec,
      method = "linear",
      yleft = 1,
      yright = tail(surv_vec, 1)
    )
    out[i, ] <- interpolator(target_times)
  }
  out
}

.register_time_to_event_model(
  engine = "glmnet",
  outcome = "survival",
  constructor = function(spec = list()) {
    defaults <- list(engine = "glmnet", package = "glmnet")
    GlmnetSurvivalModel$new(spec = modifyList(defaults, spec))
  },
  packages = "glmnet",
  tags = c("penalised", "elastic-net"),
  label = "Penalised Cox (glmnet)"
)
