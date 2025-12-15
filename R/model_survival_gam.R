#' Generalized Additive Survival Model (R6)
#'
#' Fits a GAM with Cox PH family using `mgcv`.
#' Includes logic for formula construction and baseline hazard estimation.
#'
#' @keywords internal
#' @noRd
GamSurvivalModel <- R6::R6Class(
  classname = "GamSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL, # The fitted GAM object
    time_grid = NULL,
    task = NULL,
    base_model = NULL, # The internal Cox model for baseline hazard
    varprof = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      self$varprof <- VariableProfile(data, task$features)

      expvars <- task$features
      timevar <- task$time_col
      eventvar <- task$event_col

      # Params
      spec <- self$spec
      shrinkTreshold <- spec$shrinkTreshold %||% 10

      # Ensure numeric event
      data[[eventvar]] <- as.numeric(data[[eventvar]] == 1)

      # Formula Construction
      # Identify variable types
      df_vars <- data[, expvars, drop = FALSE]
      is_num <- sapply(df_vars, is.numeric)
      is_fac <- sapply(df_vars, function(x) is.factor(x) || is.character(x))

      numvars <- expvars[is_num]
      fctvars <- expvars[is_fac]

      # Convert char to factor
      for (v in fctvars) {
        if (is.character(data[[v]])) data[[v]] <- as.factor(data[[v]])
      }

      # Splitting based on threshold
      cat_shrink <- character(0)
      cat_noshrink <- character(0)
      if (length(fctvars) > 0) {
        lvls <- sapply(data[, fctvars, drop = FALSE], function(x) length(levels(as.factor(x))))
        cat_shrink <- fctvars[lvls > shrinkTreshold]
        cat_noshrink <- fctvars[lvls <= shrinkTreshold]
      }

      num_smooth <- character(0)
      num_linear <- character(0)
      if (length(numvars) > 0) {
        uniques <- sapply(data[, numvars, drop = FALSE], function(x) length(unique(x[!is.na(x)])))
        num_smooth <- numvars[uniques > shrinkTreshold]
        num_linear <- numvars[uniques <= shrinkTreshold]
      }

      terms <- c()
      # Random effect for high-cardinality factors
      if (length(cat_shrink) > 0) terms <- c(terms, paste0("s(", cat_shrink, ", bs='re')"))
      # Smooth for numerics
      if (length(num_smooth) > 0) terms <- c(terms, paste0("s(", num_smooth, ")"))
      # Linear terms
      terms <- c(terms, num_linear, cat_noshrink)

      rhs <- if (length(terms) > 0) paste(terms, collapse = "+") else "1"
      formula <- stats::as.formula(paste(timevar, "~", rhs))

      # Fit GAM
      b <- tryCatch(
        {
          mgcv::gam(
            formula,
            family = mgcv::cox.ph(),
            data = data,
            weights = data[[eventvar]], # Weights are status? standard mgcv::cox.ph usage
            select = TRUE
          )
        },
        error = function(e) {
          rlang::abort(glue::glue("GAM fit failed: {e$message}"))
        }
      )

      # Baseline Hazard Estimation
      # Predict link (linear predictor)
      lp <- stats::predict(b, newdata = data, type = "link")

      # Fit simple Cox on LP to get baseline
      # Datasurv
      surv_df <- data.frame(time = data[[timevar]], event = data[[eventvar]], score = lp)

      # score2proba logic inline
      # Fit coxph(Surv(time, event) ~ score, init=1, iter.max=0)
      # This fixes beta=1 for score, estimating baseline hazard H0(t) for exp(1*score)
      # Actually, mgcv::cox.ph provides martingale residuals etc, but getting S(t) requires this trick often.

      base_cox <- survival::coxph(
        survival::Surv(time, event) ~ score,
        data = surv_df,
        init = 1,
        control = survival::coxph.control(iter.max = 0)
      )

      self$model <- b
      self$base_model <- base_cox

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)

      # Predict LP
      lp <- stats::predict(self$model, newdata = complete_data, type = "link", se.fit = FALSE)

      # Get Survival using baseline model
      # survfit(base_cox, newdata=data.frame(score=lp))
      # This returns S(t) for each obs

      sf <- survival::survfit(
        self$base_model,
        newdata = data.frame(score = lp),
        conf.int = 0.95
      )

      # sf$surv is matrix [times, obs] usually (if newdata has multiple rows)
      # Check dim
      surv_mat <- sf$surv
      if (is.null(dim(surv_mat))) {
        # If 1 obs, it is vector [times]
        surv_mat <- matrix(surv_mat, ncol = 1)
      }

      # Times from sf
      base_times <- sf$time

      # Add 0 if needed
      if (!0 %in% base_times) {
        base_times <- c(0, base_times)
        surv_mat <- rbind(rep(1, ncol(surv_mat)), surv_mat)
      }

      # Interpolate
      target_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Matrix is [times, obs]. We need transpose for Interpolator if it expects [obs, times]??
      # Checked earlier: survprobMatInterpolator expects [times, obs]? (rows=times)
      # Yes, existing code SurvModel_GAM used estSURVTest which was sf$surv
      # Let's verify interpolator arg.
      # But previous R6 implementations used transposes carefully.
      # survprobMatInterpolator: "probsMat: matrix of probabilities (rows = times, cols = observations)"
      # So sf$surv (rows=times) is correct.

      interp_mat <- survprobMatInterpolator(surv_mat, base_times, target_times)

      # Flatten
      # [times, obs] -> we want [id1_t1, id1_t2...] which is column-major flatten of [times, obs].
      # So as.vector(interp_mat) -> t1_obs1, t2_obs1... Correct.

      id_col <- complete_data[[self$task$id_col]]

      new_survival_prediction(
        id = rep(id_col, each = length(target_times)),
        time = rep(target_times, times = length(id_col)),
        surv = as.vector(interp_mat),
        model = rep("gam", length(id_col) * length(target_times)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "Generalized Additive Model"
      info
    },
    required_packages = function() {
      c("mgcv", "survival")
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
  engine = "gam",
  outcome = "survival",
  constructor = function(spec = list()) {
    GamSurvivalModel$new(spec = modifyList(list(engine = "gam"), spec))
  },
  packages = "mgcv",
  tags = c("smooth", "semiparametric"),
  label = "Generalised additive model"
)
