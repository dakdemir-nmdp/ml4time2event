#' XGBoost Survival Model (R6)
#'
#' Fits a Cox Proportional Hazards model using XGBoost.
#'
#' @keywords internal
#' @noRd
XGBoostSurvivalModel <- R6::R6Class(
  classname = "XGBoostSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL, # Stores the xgb.Booster
    baseline_hazard = NULL, # Stores the estimated baseline hazard
    varprof = NULL,
    feature_names = NULL, # Store for alignment
    task = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task

      data <- as.data.frame(task$data)
      self$varprof <- VariableProfile(data, task$features)

      # Extract parameters
      spec <- self$spec
      nrounds <- spec$nrounds %||% 100
      eta <- spec$eta %||% 0.01
      max_depth <- spec$max_depth %||% 5
      verbose <- spec$verbose %||% 0

      # Prepare data
      # Convert data to DMatrix format
      # Label: negative time for censored, positive for event
      time_col <- task$time_col
      event_col <- task$event_col

      # Model matrix transformation for factors
      # We need to save the formula/structure for prediction alignment
      # Actually, xgboost just needs a numeric matrix.
      # We can use model.matrix.

      X <- stats::model.matrix(~ -1 + ., data[, task$features, drop = FALSE])
      self$feature_names <- colnames(X)

      y_times <- data[[time_col]]
      y_events <- data[[event_col]] # 0/1

      # XGBoost Cox objective expects:
      # label = time (if event), -time (if censored)
      # But strictly: time > 0 is event, time < 0 is censored?
      # Wait, standard xgb Cox:
      # "The label should be the time to event, with negative values indicating censoring."
      # So if t=10, event=1 -> label=10
      # if t=10, event=0 -> label=-10

      y_label <- ifelse(y_events == 1, y_times, -y_times)

      dtrain <- xgboost::xgb.DMatrix(data = X, label = y_label)

      params <- list(
        objective = "survival:cox",
        eval_metric = "cox-nloglik",
        eta = eta,
        max_depth = max_depth
      )

      # Train
      bst <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = verbose
      )

      self$model <- bst

      # Estimate Baseline Hazard (Breslow)
      # Re-use logic from legacy xgb.train.surv

      # Predicted linear predictor (log hazard ratio)
      lp_train <- predict(bst, X)

      df_train <- data.frame(
        time = y_times,
        status = y_events,
        exp_lp = exp(lp_train)
      )

      unique_event_times <- sort(unique(df_train$time[df_train$status == 1]))

      if (length(unique_event_times) == 0) {
        # Minimal fallback
        self$baseline_hazard <- data.frame(time = c(0, max(y_times)), hazard = c(0, 0.01))
      } else {
        # Breslow
        # h0(t) = d_t / sum_{j \in R_t} exp(lp_j)

        # Optimize loop?
        # Sort by time to make risk sets easier
        df_train <- df_train[order(df_train$time), ]

        # This can be vectorized or computed efficiently:
        # Denominator is sum of exp_lp for all having t_j >= t
        # We can compute reverse cumulative sum of exp_lp (if sorted descending)?
        # Actually risk set R_t is t_j >= t.
        # So it's suffix sums.

        # Sort data descending by time for cumsum
        df_desc <- df_train[order(df_train$time, decreasing = TRUE), ]
        total_risk <- cumsum(df_desc$exp_lp)
        # Map back to increasing time order for unique event times?
        # Let's stick strictly to unique_event_times.

        # Pre-calc risk sums at each unique event time
        # This is strictly: sum(exp_lp where time >= t)

        # For simplicity and correctness with vectorization:
        # unique_event_times are sorted ascending.

        # Let's map unique times to indices in the sorted descending array?
        # Or just use the robust loop if N is small?
        # Original code used a loop 'sum(data$time >= t)'.
        # Let's optimize slightly:

        risk_denominators <- vapply(unique_event_times, function(t) {
          sum(df_train$exp_lp[df_train$time >= t])
        }, numeric(1))

        # Number of events at each time
        # This is `table` or `tapply`
        events_at_t <- as.vector(table(factor(df_train$time[df_train$status == 1], levels = unique_event_times)))

        h0_vals <- ifelse(risk_denominators > 0, events_at_t / risk_denominators, 0)

        bh <- data.frame(
          time = unique_event_times,
          hazard = h0_vals
        )

        # Ensure 0,0 start
        if (bh$time[1] != 0) {
          bh <- rbind(data.frame(time = 0, hazard = 0), bh)
        }

        self$baseline_hazard <- bh
      }

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()

      # Prepare newdata
      # Must match training model matrix columns
      # Be careful with factors: missing levels logic from legacy

      complete_data <- .ensure_prediction_data(newdata, self$task)

      # Model Matrix
      # Handling new levels or missing columns:
      # We construct model matrix on newdata.
      # Then align columns with self$feature_names.

      # Basic model matrix on new data
      # Note: this might fail if newdata has new factor levels not seen in training
      # if model.matrix tries to encode them.
      # Or produce different columns.
      # Robust way: bind with dummy row?
      # Or just let model.matrix run and drop/add cols.
      # For simplicity:

      X_new <- stats::model.matrix(~ -1 + ., complete_data[, self$task$features, drop = FALSE])

      # Align columns
      missing_cols <- setdiff(self$feature_names, colnames(X_new))
      if (length(missing_cols) > 0) {
        # Add missing as 0
        matrix_add <- matrix(0, nrow = nrow(X_new), ncol = length(missing_cols))
        colnames(matrix_add) <- missing_cols
        X_new <- cbind(X_new, matrix_add)
      }

      # Drop extra cols
      extra_cols <- setdiff(colnames(X_new), self$feature_names)
      if (length(extra_cols) > 0) {
        X_new <- X_new[, !colnames(X_new) %in% extra_cols, drop = FALSE]
      }

      # Reorder
      X_new <- X_new[, self$feature_names, drop = FALSE]

      # Predict LP
      dtest <- xgboost::xgb.DMatrix(data = X_new)
      lp <- predict(self$model, dtest)

      # Compute Survival
      # S(t) = exp(-H0(t) * exp(lp))

      # Eval times
      if (is.null(times)) {
        # Use training event times
        eval_times <- self$baseline_hazard$time
      } else {
        eval_times <- sort(unique(c(0, as.numeric(times))))
      }

      # Compute Cumulative Hazard H0(t) at eval_times
      # sum of hazard increments <= t
      # Vectorize this lookup

      bh <- self$baseline_hazard
      # Cumulative sum of hazards in baseline reference
      # wait, bh$hazard IS the increment h0_t.
      # So H0(t_k) = sum_{j: t_j <= t_k} h0(t_j)

      bh$cumhaz <- cumsum(bh$hazard)

      # Step function lookup for H0 at eval_times
      # approx(..., method="constant", f=0 or 1)
      # method="constant", f=1 implies right-continuous (value at step)
      # We want value at t. If t matches step, take it.
      # stepfun logic.

      H0_t <- stats::approx(bh$time, bh$cumhaz, xout = eval_times, method = "constant", rule = 2, f = 1)$y

      # S matrix: rows=obs, cols=times
      exp_lp <- exp(lp)

      # Outer: exp_lp (n) * H0_t (m)
      cumhaz_matrix <- outer(exp_lp, H0_t)
      surv_matrix <- exp(-cumhaz_matrix)

      # Interpolate to requested times if we added 0 or sorted differently?
      # eval_times matches columns of surv_matrix.
      # If `times` was passed, we specifically computed for those (plus maybe 0).
      # If `times` was NULL, we used internal grid.

      if (!is.null(times)) {
        # We computed at eval_times = sort(unique(c(0, times)))
        # But we need to return exactly `times` requested?
        # No, usually predict_survival returns "tidy" predictions for requested times.

        # If `times` passed, we just selecting columns corresponding to `times`.
        # Note: we added 0 to eval_times.

        # Use utility interpolator if needed, but we essentially evaluated exactly at points.
        # Just subset columns?
        # But outer product assumes exact match?
        # Yes.
        # If user passed distinct times, we have them in eval_times.

        # This logic is simpler:
        # result times = unique(sort(times)).
        # But predict_survival interface usually allows general `times`.

        # Let's follow standard pattern: return grid result or specific
        # prediction object format asks for long format.
      } else {
        # Return all grid points
      }

      # Convert to matrix for gather
      # We used eval_times.

      # Use interpolator to map from eval_times to exactly `times` (or default grid)
      req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Transpose for interpolator: rows=times, cols=obs
      surv_t <- t(surv_matrix)

      interp_t <- survprobMatInterpolator(surv_t, eval_times, req_times)

      # Flatten
      id_col <- complete_data[[self$task$id_col]]
      n_samples <- nrow(complete_data)
      n_times <- length(req_times)

      new_survival_prediction(
        id = rep(id_col, each = n_times),
        time = rep(req_times, times = n_samples),
        surv = as.vector(interp_t),
        model = rep("xgboost", n_samples * n_times),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "XGBoost Cox"
      info
    },
    required_packages = function() {
      c("xgboost")
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
  engine = "xgboost",
  outcome = "survival",
  constructor = function(spec = list()) {
    XGBoostSurvivalModel$new(spec = modifyList(list(engine = "xgboost"), spec))
  },
  packages = "xgboost",
  label = "XGBoost Cox",
  tags = c("xgboost", "cox")
)
