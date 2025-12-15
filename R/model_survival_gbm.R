#' Gradient Boosted Survival Model (R6)
#'
#' Fits a gradient boosting machine (GBM) with Cox proportional hazards loss.
#'
#' @keywords internal
#' @noRd
GbmSurvivalModel <- R6::R6Class(
  classname = "GbmSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    varprof = NULL,
    factor_levels = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      expvars <- task$features
      timevar <- task$time_col
      eventvar <- task$event_col

      spec <- self$spec
      ntree <- spec$ntree %||% 200
      interaction.depth <- spec$max.depth %||% 3
      shrinkage <- spec$learning_rate %||% spec$shrinkage %||% 0.01
      bag.fraction <- spec$bag.fraction %||% 0.5
      train.fraction <- spec$train.fraction %||% 1.0 # Default to use all data if no validation set provided
      n.minobsinnode <- spec$n.minobsinnode %||% 10

      # 1. Profile and Factor Levels
      self$varprof <- VariableProfile(data, expvars)

      self$factor_levels <- list()
      for (v in expvars) {
        if (is.factor(data[[v]]) || is.character(data[[v]])) {
          self$factor_levels[[v]] <- levels(as.factor(data[[v]]))
        }
      }

      # Prepare data for GBM
      # GBM handles factors internally? Yes.
      # Event must be 0/1 numeric? "coxph" requires numeric event.
      data[[eventvar]] <- as.numeric(data[[eventvar]] == 1)

      # Adaptive parameters (from legacy logic)
      n_obs <- nrow(data)
      cv_folds <- 0 # Default disabled

      if (n_obs >= 100) {
        # Ensure bag.fraction * n_train > 2 * minobs + 1
        # If user didn't specify strict minobs, adapt it
        if (is.null(spec$n.minobsinnode)) {
          safe_bag <- min(bag.fraction, 0.8)
          max_safe_minobs <- floor((n_obs * safe_bag - 1) / 2)
          n.minobsinnode <- max(1, min(5, max_safe_minobs))
          bag.fraction <- safe_bag
        }
      } else {
        n.minobsinnode <- 1
        bag.fraction <- 1.0
      }

      # Formula
      form <- as.formula(paste("Surv(", timevar, ",", eventvar, ") ~ ."))
      # Subset data columns
      train_data <- data[, c(timevar, eventvar, expvars), drop = FALSE]

      message("Fitting GBM Survival...")

      fit_res <- tryCatch(
        {
          suppressMessages(gbm::gbm(
            formula = form,
            data = train_data,
            distribution = "coxph",
            n.trees = ntree,
            interaction.depth = interaction.depth,
            n.minobsinnode = n.minobsinnode,
            shrinkage = shrinkage,
            bag.fraction = bag.fraction,
            train.fraction = train.fraction,
            cv.folds = cv_folds,
            keep.data = TRUE,
            verbose = FALSE
          ))
        },
        error = function(e) {
          warning("GBM failed with default settings. Retrying with minimal params.")
          suppressMessages(gbm::gbm(
            formula = form,
            data = train_data,
            distribution = "coxph",
            n.trees = min(50, ntree),
            interaction.depth = 1,
            n.minobsinnode = 1,
            shrinkage = 0.01,
            bag.fraction = 1.0,
            train.fraction = 1.0,
            cv.folds = 0,
            keep.data = TRUE,
            verbose = FALSE
          ))
        }
      )

      # Determine best trees
      best_iter <- fit_res$n.trees
      if (cv_folds > 0 && n_obs >= 50) { # Only use perf if we did CV
        try({
          best_iter <- gbm::gbm.perf(fit_res, method = "cv", plot.it = FALSE)
        })
      }
      # If bag fraction < 1, could use OOB?
      if (bag.fraction < 1 && cv_folds == 0 && n_obs >= 50) {
        try({
          best_iter <- gbm::gbm.perf(fit_res, method = "OOB", plot.it = FALSE)
        })
      }

      fit_res$n.trees.best <- best_iter

      # Baseline Hazard (Breslow via gbm::basehaz.gbm)
      pred_train <- predict(fit_res, train_data, n.trees = best_iter, type = "link")

      u_times <- sort(unique(data[[timevar]][data[[eventvar]] == 1]))
      base_cum <- gbm::basehaz.gbm(
        t = data[[timevar]],
        delta = data[[eventvar]],
        f.x = pred_train,
        t.eval = u_times,
        cumulative = TRUE
      )

      self$model <- list(
        gbm = fit_res,
        baseline_hazard = data.frame(time = u_times, cumhaz = base_cum),
        best_iter = best_iter
      )

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)
      expvars <- self$task$features

      # Factor alignment
      for (v in expvars) {
        if (v %in% names(self$factor_levels)) {
          complete_data[[v]] <- factor(complete_data[[v]], levels = self$factor_levels[[v]])
        }
      }

      preds_lp <- predict(self$model$gbm, newdata = complete_data, n.trees = self$model$best_iter, type = "link")

      # Surv = exp(- cumhaz * exp(lp))
      bh_times <- c(0, self$model$baseline_hazard$time)
      bh_cum <- c(0, self$model$baseline_hazard$cumhaz)

      target_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Intepolate CumHaz
      bh_interp <- stats::approx(bh_times, bh_cum, xout = target_times, rule = 2)$y

      # Matrix: [ntimes x nobs]
      Lambda <- outer(bh_interp, exp(preds_lp))
      Surv <- exp(-Lambda)

      id_col <- complete_data[[self$task$id_col]]
      flat_surv <- as.vector(Surv)

      new_survival_prediction(
        id = rep(id_col, each = length(target_times)),
        time = rep(target_times, times = length(id_col)),
        surv = flat_surv,
        model = rep("gbm", length(flat_surv)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "Gradient Boosted Survival"
      info
    },
    required_packages = function() {
      c("gbm")
    }
  ),
  private = list(
    ensure_fitted = function() {
      if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
    }
  )
)

.register_time_to_event_model(
  engine = "gbm",
  outcome = "survival",
  constructor = function(spec = list()) {
    GbmSurvivalModel$new(spec = modifyList(list(engine = "gbm"), spec))
  },
  packages = "gbm",
  tags = c("boosting", "tree-based"),
  label = "Gradient Boosted Survival"
)
