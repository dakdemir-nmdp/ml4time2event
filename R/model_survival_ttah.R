#' TTAH Survival Model (R6)
#'
#' Fits a discrete-time hazard model with time-varying coefficients and
#' latent feature interactions (Time-to-Event Additive Hazards?).
#' Uses ridge-penalized logistic regression on person-time data.
#'
#' @keywords internal
#' @noRd
TtahSurvivalModel <- R6::R6Class(
  classname = "TtahSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    varprof = NULL,
    factor_levels = NULL,
    basis_specs = NULL,
    latent_projection = NULL,
    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      expvars <- task$features
      timevar <- task$time_col
      eventvar <- task$event_col

      # Params
      spec <- self$spec
      n_time <- spec$n_time %||% 50
      spline_knots <- spec$spline_knots %||% 5
      latent_dim <- spec$latent_dim %||% 8
      time_basis_df <- spec$time_basis_df %||% 4
      lambda <- spec$lambda %||% 1e-3

      # 1. Profile and Prep
      self$varprof <- VariableProfile(data, expvars)

      self$factor_levels <- list()
      for (v in expvars) {
        if (is.factor(data[[v]]) || is.character(data[[v]])) {
          self$factor_levels[[v]] <- levels(as.factor(data[[v]]))
        }
      }

      obs_times <- data[[timevar]]
      obs_events <- as.numeric(data[[eventvar]] == 1)

      # Grid - Filter positive values if time_grid provided
      tg_pass <- if (!is.null(time_grid)) time_grid[time_grid > 0] else NULL
      if (!is.null(tg_pass) && length(tg_pass) == 0) tg_pass <- NULL

      grid <- ttah_build_time_grid(obs_times, time_grid = tg_pass, n_time = n_time)

      # Features
      prep <- ttah_prepare_features(
        data = data[, expvars, drop = FALSE],
        expvars = expvars,
        factor_levels = self$factor_levels,
        basis_specs = NULL,
        spline_knots = spline_knots
      )
      phi <- prep$phi
      self$basis_specs <- prep$basis_specs
      self$factor_levels <- prep$factor_levels # Updated with any new knowledge?

      # Latent
      self$latent_projection <- ttah_compute_latent_projection(phi, latent_dim = latent_dim)
      phi_latent <- if (ncol(self$latent_projection) > 0) {
        tmp <- phi %*% self$latent_projection
        colnames(tmp) <- colnames(self$latent_projection)
        tmp
      } else {
        NULL
      }

      # Time Basis
      time_basis <- ttah_time_basis(grid, df = time_basis_df)
      time_basis_matrix <- time_basis$matrix
      if (ncol(time_basis_matrix) > 0) {
        colnames(time_basis_matrix) <- paste0("time_basis", seq_len(ncol(time_basis_matrix)))
      }
      time_basis_names <- colnames(time_basis_matrix)

      # Long Design
      interval_index <- ttah_assign_intervals(obs_times, grid)

      long_design <- ttah_build_long_design(
        phi = phi,
        phi_latent = phi_latent,
        time_basis = time_basis_matrix,
        interval_index = interval_index,
        event = obs_events
      )

      X <- cbind(Intercept = 1, long_design$X)

      fit_res <- ttah_ridge_glm(
        X = X,
        y = long_design$y,
        lambda = lambda
      )

      # Extract coefs
      coef <- fit_res$coefficients
      coef[is.na(coef)] <- 0

      offset <- 1
      coef_add <- coef[offset + seq_len(ncol(phi))]
      offset <- offset + ncol(phi)

      coef_time <- coef[offset + seq_len(ncol(time_basis_matrix))]
      offset <- offset + ncol(time_basis_matrix)

      latent_coef <- numeric(0)
      if (ncol(self$latent_projection) > 0) {
        latent_len <- ncol(self$latent_projection) * ncol(time_basis_matrix)
        if (latent_len > 0) {
          latent_coef <- coef[offset + seq_len(latent_len)]
        }
      }

      self$model <- list(
        intercept = coef[1],
        coef_add = coef_add,
        coef_time = coef_time,
        coef_inter = latent_coef,
        feature_names = colnames(phi),
        time_basis_names = time_basis_names,
        time_basis_specs = time_basis$specs,
        grid = grid
      )

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)
      expvars <- self$task$features

      # 1. Factor Alignment
      for (v in expvars) {
        if (v %in% names(self$factor_levels)) {
          complete_data[[v]] <- factor(complete_data[[v]], levels = self$factor_levels[[v]])
        }
      }

      # 2. Features
      prep_new <- ttah_prepare_features(
        data = complete_data,
        expvars = expvars,
        factor_levels = self$factor_levels,
        basis_specs = self$basis_specs
      )
      phi_new <- prep_new$phi
      # Ensure cols match
      phi_new <- phi_new[, self$model$feature_names, drop = FALSE]

      # 3. Latent
      latent_dim <- if (is.null(self$latent_projection)) 0 else ncol(self$latent_projection)
      phi_latent_new <- if (latent_dim > 0) {
        phi_new %*% self$latent_projection
      } else {
        NULL
      }

      # 4. Time Basis (on Model Grid)
      time_grid <- self$model$grid
      time_basis <- ttah_eval_time_basis(time_grid, self$model$time_basis_specs)

      n_obs <- nrow(phi_new)
      K <- length(time_grid)

      # 5. Compute Linear Predictor (Eta)
      # Eta is K x N (Time x Obs). Wait, code below constructs K x N.

      # Intercept
      eta_matrix <- matrix(self$model$intercept, nrow = K, ncol = n_obs)

      # Main Effects (Phi * Beta_add). N x 1 -> 1 x N. Add to all times.
      main_eff <- as.vector(phi_new %*% self$model$coef_add)
      eta_matrix <- eta_matrix + matrix(main_eff, nrow = K, ncol = n_obs, byrow = TRUE)

      # Time Effects (TimeBasis * Beta_time). K x 1. Add to all Obs.
      time_lin <- drop(time_basis %*% self$model$coef_time)
      eta_matrix <- eta_matrix + matrix(time_lin, nrow = K, ncol = n_obs, byrow = FALSE)

      # Interaction (Latent_i * Theta * TimeBasis_k)
      if (!is.null(phi_latent_new) && length(self$model$coef_inter) > 0) {
        # Reshape Coef_inter
        theta <- matrix(self$model$coef_inter, nrow = ncol(phi_latent_new), ncol = ncol(time_basis))
        # (N x d) * (d x dt) * (dt x K)^T
        # N x K
        interaction <- phi_latent_new %*% theta %*% t(time_basis)
        eta_matrix <- eta_matrix + t(interaction)
      }

      # 6. Hazard & Survival
      hazard_matrix <- ttah_sigmoid(eta_matrix)
      survival_body <- apply(1 - hazard_matrix, 2, cumprod)
      if (!is.matrix(survival_body)) survival_body <- matrix(survival_body, ncol = n_obs)

      Probs <- rbind(rep(1, n_obs), survival_body) # times: 0, t1, t2...
      Times <- c(0, time_grid)

      req_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # Interpolate
      interp_probs <- survprobMatInterpolator(Probs, Times, req_times)

      id_col <- complete_data[[self$task$id_col]]
      flat_surv <- as.vector(interp_probs)

      new_survival_prediction(
        id = rep(id_col, each = length(req_times)),
        time = rep(req_times, times = length(id_col)),
        surv = flat_surv,
        model = rep("ttah", length(flat_surv)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "TTAH Survival"
      info
    },
    required_packages = function() {
      character()
    } # Uses ttah_utils (internal)
  ),
  private = list(
    ensure_fitted = function() {
      if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
    }
  )
)

.register_time_to_event_model(
  engine = "ttah",
  outcome = "survival",
  constructor = function(spec = list()) {
    TtahSurvivalModel$new(spec = modifyList(list(engine = "ttah"), spec))
  },
  packages = character(),
  tags = c("semiparametric", "discrete-time"),
  label = "TTAH Survival"
)
