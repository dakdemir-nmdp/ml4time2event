#' @title Sigmoid helper for logistic link
#' @param x numeric vector or matrix
#' @return element-wise logistic transform
#' @noRd
ttah_sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

#' @title Ridge-stabilised logistic regression via IRLS
#' @param X design matrix with intercept column
#' @param y binary response vector
#' @param lambda ridge penalty applied to coefficients excluding intercept
#' @param max_iter maximum number of IRLS iterations
#' @param tol convergence tolerance on coefficient changes
#' @return list with coefficients, converged flag, iterations, deviance
#' @noRd
ttah_ridge_glm <- function(X, y, lambda = 1e-3, max_iter = 100, tol = 1e-6) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(y) || any(!y %in% c(0, 1))) {
    stop("Logistic regression expects binary outcomes in {0, 1}.")
  }

  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  penalty_matrix <- diag(p)
  penalty_matrix[1, 1] <- 0  # do not penalise intercept

  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    eta <- as.vector(X %*% beta)
    mu <- ttah_sigmoid(pmin(pmax(eta, -20), 20))  # clamp to avoid overflow in exp
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)          # stability guard
    W <- mu * (1 - mu)
    z <- eta + (y - mu) / W

    WX <- X * sqrt(W)
    z_weighted <- sqrt(W) * z

    system_matrix <- crossprod(WX) + lambda * penalty_matrix
    rhs <- crossprod(WX, z_weighted)

    beta_new <- tryCatch(
      solve(system_matrix, rhs),
      error = function(e) {
        # Pseudo-inverse fallback for near-singular systems
        qr.solve(system_matrix, rhs)
      }
    )

    if (max(abs(beta_new - beta)) < tol) {
      beta <- as.vector(beta_new)
      converged <- TRUE
      break
    }

    beta <- as.vector(beta_new)
  }

  eta_final <- as.vector(X %*% beta)
  mu_final <- ttah_sigmoid(pmin(pmax(eta_final, -20), 20))
  mu_final <- pmin(pmax(mu_final, 1e-12), 1 - 1e-12)
  loglik <- sum(y * log(mu_final) + (1 - y) * log(1 - mu_final))
  deviance <- -2 * loglik

  list(
    coefficients = beta,
    converged = converged,
    iterations = iter,
    deviance = deviance
  )
}

#' @title Construct discrete time grid for TTAH models
#' @param time_values numeric vector of observed times
#' @param time_grid optional user-supplied grid (sorted within)
#' @param n_time desired number of intervals when automatically constructing grid
#' @return sorted numeric vector of time grid points (length >= 1)
#' @noRd
ttah_build_time_grid <- function(time_values, time_grid = NULL, n_time = 50) {
  if (!is.null(time_grid)) {
    if (!is.numeric(time_grid) || any(time_grid <= 0)) {
      stop("'time_grid' must be a numeric vector of positive values.")
    }
    return(sort(unique(time_grid)))
  }

  time_values <- time_values[is.finite(time_values) & time_values > 0]
  if (length(time_values) == 0) {
    stop("Cannot infer time grid: no positive finite time values supplied.")
  }

  probs <- seq(0.05, 0.95, length.out = max(3, min(n_time, length(time_values))))
  grid <- unique(stats::quantile(time_values, probs = probs, type = 7, na.rm = TRUE))
  grid <- grid[grid > 0]
  grid <- sort(unique(c(grid, max(time_values))))
  grid
}

#' @title Assign observations to discrete intervals
#' @param times numeric vector of observed times
#' @param grid numeric time grid (ascending)
#' @return integer vector of interval indices (>= 1)
#' @noRd
ttah_assign_intervals <- function(times, grid) {
  if (length(grid) == 0) {
    stop("Time grid must contain at least one value.")
  }
  idx <- findInterval(times, grid, rightmost.closed = TRUE)
  idx[idx < 1] <- 1
  idx[idx > length(grid)] <- length(grid)
  idx
}

#' @title Build B-spline time basis on grid points
#' @param grid numeric time grid
#' @param df desired degrees of freedom (excluding intercept)
#' @param degree spline degree
#' @return list with basis matrix and specifications for reuse
#' @noRd
ttah_time_basis <- function(grid, df = 4, degree = 3) {
  df_use <- min(df, length(grid))
  if (df_use < 1) {
    df_use <- 1
  }
  if (df_use == 1) {
    basis <- matrix(1, nrow = length(grid), ncol = 1)
    specs <- list(
      df = 1,
      degree = 0,
      knots = numeric(0),
      boundary = range(grid)
    )
  } else {
    basis <- splines::bs(grid, df = df_use, degree = degree, intercept = FALSE)
    specs <- list(
      df = df_use,
      degree = degree,
      knots = attr(basis, "knots"),
      boundary = attr(basis, "Boundary.knots")
    )
  }
  list(matrix = basis, specs = specs)
}

#' @title Evaluate stored B-spline basis at new times
#' @param times numeric vector of times
#' @param specs list with knots/boundary/degree information
#' @return matrix with rows = length(times), columns = df
#' @noRd
ttah_eval_time_basis <- function(times, specs) {
  if (is.null(specs) || is.null(specs$df)) {
    stop("Time basis specifications are missing.")
  }

  if (specs$df <= 1) {
    return(matrix(1, nrow = length(times), ncol = 1))
  }

  splines::bs(
    times,
    df = specs$df,
    degree = specs$degree,
    knots = specs$knots,
    Boundary.knots = specs$boundary,
    intercept = FALSE
  )
}

#' @title Prepare feature basis expansions
#' @param data data.frame with predictor variables
#' @param expvars character vector of predictors
#' @param factor_levels optional named list of factor levels for prediction-time alignment
#' @param basis_specs optional list of training-time basis parameters
#' @param spline_knots number of knots/df for numeric splines
#' @return list with matrix `phi`, updated `basis_specs`, and `factor_levels`
#' @noRd
ttah_prepare_features <- function(data, expvars, factor_levels = NULL, basis_specs = NULL,
                                  spline_knots = 5) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  if (length(expvars) == 0) {
    stop("'expvars' must contain at least one variable.")
  }

  phi_list <- list()
  spec_list <- list()
  levels_list <- list()

  for (var in expvars) {
    col <- data[[var]]
    if (is.null(basis_specs)) {
      spec_list[[var]] <- list()
    } else if (!is.null(basis_specs[[var]])) {
      spec_list[[var]] <- basis_specs[[var]]
    } else {
      spec_list[[var]] <- list()
    }

    if (is.numeric(col)) {
      if (is.null(basis_specs) || is.null(spec_list[[var]]$type)) {
        mean_val <- mean(col, na.rm = TRUE)
        sd_val <- stats::sd(col, na.rm = TRUE)
        if (!is.finite(sd_val) || sd_val == 0) {
          sd_val <- 1
        }
        scaled <- (col - mean_val) / sd_val
        df_use <- min(max(1, spline_knots), length(unique(stats::na.omit(col))) - 1)
        if (!is.finite(df_use) || df_use < 1) {
          df_use <- 1
        }
        if (df_use == 1) {
          basis_mat <- matrix(scaled, ncol = 1)
        } else {
          basis_mat <- splines::ns(scaled, df = df_use, intercept = FALSE)
        }
        spec_list[[var]] <- list(
          type = "numeric",
          mean = mean_val,
          sd = sd_val,
          df = df_use,
          knots = attr(basis_mat, "knots"),
          boundary = attr(basis_mat, "Boundary.knots")
        )
      } else {
        scaled <- (col - spec_list[[var]]$mean) / spec_list[[var]]$sd
        df_use <- spec_list[[var]]$df
        if (df_use == 1) {
          basis_mat <- matrix(scaled, ncol = 1)
        } else {
          basis_mat <- splines::ns(
            scaled,
            df = df_use,
            knots = spec_list[[var]]$knots,
            Boundary.knots = spec_list[[var]]$boundary,
            intercept = FALSE
          )
        }
      }
      colnames(basis_mat) <- paste0(var, "_basis", seq_len(ncol(basis_mat)))
      phi_list[[var]] <- basis_mat
    } else if (is.factor(col) || is.character(col)) {
      if (!is.null(factor_levels) && !is.null(factor_levels[[var]])) {
        levels_list[[var]] <- factor_levels[[var]]
      } else if (is.factor(col)) {
        levels_list[[var]] <- levels(col)
      } else {
        levels_list[[var]] <- sort(unique(col))
      }
      col_factor <- factor(col, levels = levels_list[[var]])
      dummy_mat <- stats::model.matrix(~ col_factor - 1)
      dummy_mat[is.na(dummy_mat)] <- 0
      colnames(dummy_mat) <- gsub("^col_factor", var, colnames(dummy_mat))
      phi_list[[var]] <- dummy_mat
      spec_list[[var]]$type <- "factor"
    } else {
      stop(sprintf("Variable '%s' must be numeric, factor, or character.", var))
    }
  }

  phi <- do.call(cbind, phi_list)
  if (!is.matrix(phi)) {
    phi <- matrix(phi, ncol = length(phi))
  }
  colnames(phi) <- make.names(colnames(phi), unique = TRUE)

  list(
    phi = phi,
    basis_specs = spec_list,
    factor_levels = levels_list
  )
}

#' @title Compute latent projection matrix
#' @param phi feature basis matrix
#' @param latent_dim desired latent dimensionality
#' @return matrix with columns spanning latent space
#' @noRd
ttah_compute_latent_projection <- function(phi, latent_dim = 8) {
  if (latent_dim <= 0 || ncol(phi) == 0) {
    return(matrix(0, ncol = 0, nrow = ncol(phi)))
  }
  latent_dim <- min(latent_dim, ncol(phi))
  centered <- scale(phi, center = TRUE, scale = FALSE)
  svd_res <- svd(centered, nu = 0, nv = latent_dim)
  proj <- svd_res$v
  colnames(proj) <- paste0("latent", seq_len(ncol(proj)))
  proj
}

#' @title Construct long-format design matrix for discrete hazards
#' @param phi feature basis matrix
#' @param phi_latent latent features matrix
#' @param time_basis matrix of time basis evaluations (K x dt)
#' @param interval_index integer vector of interval assignments per observation
#' @param event indicator vector (0/1)
#' @return list with design matrix X (without intercept) and response vector y
#' @noRd
ttah_build_long_design <- function(phi, phi_latent, time_basis, interval_index, event) {
  n <- nrow(phi)
  p <- ncol(phi)
  d <- if (is.null(phi_latent)) 0 else ncol(phi_latent)
  dt <- ncol(time_basis)
  total_rows <- sum(interval_index)

  X <- matrix(0, nrow = total_rows, ncol = p + dt + d * dt)
  y <- numeric(total_rows)
  obs_id <- integer(total_rows)
  time_id <- integer(total_rows)

  cursor <- 0
  for (i in seq_len(n)) {
    k_max <- interval_index[i]
    for (k in seq_len(k_max)) {
      cursor <- cursor + 1
      X[cursor, 1:p] <- phi[i, ]
      X[cursor, p + seq_len(dt)] <- time_basis[k, ]
      if (d > 0) {
        interaction_vec <- as.vector(phi_latent[i, , drop = FALSE] %o% time_basis[k, ])
        X[cursor, (p + dt + 1):(p + dt + d * dt)] <- interaction_vec
      }
      y[cursor] <- as.numeric(event[i] == 1 && k == k_max)
      obs_id[cursor] <- i
      time_id[cursor] <- k
    }
  }

  time_names <- colnames(time_basis)
  if (is.null(time_names) || length(time_names) == 0) {
    time_names <- paste0("time_basis", seq_len(dt))
  }
  latent_names <- if (d > 0 && !is.null(colnames(phi_latent))) {
    colnames(phi_latent)
  } else if (d > 0) {
    paste0("latent", seq_len(d))
  } else {
    character(0)
  }
  interaction_names <- if (d > 0) {
    as.vector(outer(latent_names, time_names, function(a, b) paste0(a, ":", b)))
  } else {
    character(0)
  }

  colnames(X) <- c(colnames(phi), time_names, interaction_names)

  list(
    X = X,
    y = y,
    obs = obs_id,
    time = time_id
  )
}

#' @title Construct long-format design for competing risks multi-class target
#' @param phi feature basis matrix
#' @param phi_latent latent features matrix
#' @param time_basis matrix of time basis evaluations
#' @param interval_index integer vector of interval assignments per observation
#' @param event vector of original event codes (0 = censor)
#' @param cause_codes character vector of event codes to model
#' @return list with design matrix X and factor target
#' @noRd
ttah_build_multiclass_design <- function(phi, phi_latent, time_basis, interval_index,
                                         event, cause_codes) {
  n <- nrow(phi)
  p <- ncol(phi)
  d <- if (is.null(phi_latent)) 0 else ncol(phi_latent)
  dt <- ncol(time_basis)
  total_rows <- sum(interval_index)

  X <- matrix(0, nrow = total_rows, ncol = p + dt + d * dt)
  target <- character(total_rows)
  obs_id <- integer(total_rows)
  time_id <- integer(total_rows)

  baseline_label <- "no_event"
  cause_labels <- paste0("cause_", cause_codes)
  valid_labels <- c(baseline_label, cause_labels)

  cursor <- 0
  for (i in seq_len(n)) {
    k_max <- interval_index[i]
    event_code <- as.character(event[i])
    event_cause <- if (!is.na(event_code) && event_code %in% cause_codes) {
      paste0("cause_", event_code)
    } else {
      baseline_label
    }

    for (k in seq_len(k_max)) {
      cursor <- cursor + 1
      X[cursor, 1:p] <- phi[i, ]
      X[cursor, p + seq_len(dt)] <- time_basis[k, ]
      if (d > 0) {
        interaction_vec <- as.vector(phi_latent[i, , drop = FALSE] %o% time_basis[k, ])
        X[cursor, (p + dt + 1):(p + dt + d * dt)] <- interaction_vec
      }
      target[cursor] <- if (event[i] != 0 && k == k_max) event_cause else baseline_label
      obs_id[cursor] <- i
      time_id[cursor] <- k
    }
  }

  time_names <- colnames(time_basis)
  if (is.null(time_names) || length(time_names) == 0) {
    time_names <- paste0("time_basis", seq_len(dt))
  }
  latent_names <- if (d > 0 && !is.null(colnames(phi_latent))) {
    colnames(phi_latent)
  } else if (d > 0) {
    paste0("latent", seq_len(d))
  } else {
    character(0)
  }
  interaction_names <- if (d > 0) {
    as.vector(outer(latent_names, time_names, function(a, b) paste0(a, ":", b)))
  } else {
    character(0)
  }
  colnames(X) <- c(colnames(phi), time_names, interaction_names)

  list(
    X = X,
    target = factor(target, levels = valid_labels),
    obs = obs_id,
    time = time_id,
    levels = valid_labels
  )
}
