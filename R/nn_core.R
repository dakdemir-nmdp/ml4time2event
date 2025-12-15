# ============================================================================
# Weight Initialization and Manipulation
# ============================================================================

#' Initialize neural network weights
#'
#' Creates random initial weights for a shallow neural network with one
#' hidden layer using sigmoid activation and linear output.
#'
#' @param n_in Integer; number of input features
#' @param n_hidden Integer; number of hidden units
#' @return List with matrices W1, W2 and vectors b1, b2
#' @keywords internal
.nn_initialize_weights <- function(n_in, n_hidden) {
    r <- 0.7
    list(
        W1 = matrix(stats::runif(n_in * n_hidden, -r, r), nrow = n_in, ncol = n_hidden),
        b1 = stats::runif(n_hidden, -r, r),
        W2 = matrix(stats::runif(n_hidden, -r, r), nrow = n_hidden, ncol = 1),
        b2 = stats::runif(1, -r, r)
    )
}

#' Unpack flattened weight vector
#'
#' Converts a single vector of parameters back into structured weight matrices
#' and bias vectors for the neural network.
#'
#' @param w Numeric vector of flattened weights
#' @param n_in Integer; number of input features
#' @param n_hidden Integer; number of hidden units
#' @return List with matrices W1, W2 and vectors b1, b2
#' @keywords internal
.nn_unpack_weights <- function(w, n_in, n_hidden) {
    len_W1 <- n_in * n_hidden
    len_b1 <- n_hidden
    len_W2 <- n_hidden

    ptr <- 0
    W1 <- matrix(w[(ptr + 1):(ptr + len_W1)], n_in, n_hidden)
    ptr <- ptr + len_W1
    b1 <- w[(ptr + 1):(ptr + len_b1)]
    ptr <- ptr + len_b1
    W2 <- matrix(w[(ptr + 1):(ptr + len_W2)], n_hidden, 1)
    ptr <- ptr + len_W2
    b2 <- w[ptr + 1]

    list(W1 = W1, b1 = b1, W2 = W2, b2 = b2)
}

# ============================================================================
# Forward Propagation
# ============================================================================

#' Forward pass through shallow neural network
#'
#' Computes hidden layer activations and final output (log relative risk)
#' using sigmoid activation in hidden layer and linear output.
#'
#' @param X Numeric matrix; design matrix (n x p)
#' @param ws List; weight structure from .nn_unpack_weights()
#' @return List with 'output' (n x 1 matrix of log risk) and
#'   'hidden_activations' (n x h matrix)
#' @keywords internal
.nn_forward_pass <- function(X, ws) {
    sigmoid <- function(z) 1 / (1 + exp(-pmin(pmax(z, -500), 500))) # Clamp for numerical stability
    Z1 <- sweep(X %*% ws$W1, 2, ws$b1, "+")
    A1 <- sigmoid(Z1)
    Z2 <- sweep(A1 %*% ws$W2, 2, ws$b2, "+")
    list(output = Z2, hidden_activations = A1)
}

# ============================================================================
# Cox Loss Optimization
# ============================================================================

#' Fit neural network with Cox partial likelihood loss
#'
#' Optimizes neural network weights using BFGS to maximize the Cox
#' partial likelihood. Handles sorting and implements efficient
#' gradient computation using cumulative sums.
#'
#' @param X_train Numeric matrix; scaled design matrix (n x p)
#' @param event Numeric vector; binary event indicator (1=event, 0=censored)
#' @param time Numeric vector; observed times (will be sorted internally)
#' @param n_features Integer; number of input features
#' @param size Integer; number of hidden units
#' @param decay Numeric; L2 penalty on weights
#' @param maxit Integer; maximum BFGS iterations
#' @return List with optimized weights and optimization info
#' @keywords internal
.nn_fit_cox_loss <- function(X_train, event, time, n_features, size, decay, maxit) {
    # Sort data by descending time for efficient risk set calculations
    ord <- order(-time, -event)
    X_sorted <- X_train[ord, , drop = FALSE]
    time_sorted <- time[ord]
    event_sorted <- event[ord]

    # Initialize weights
    initial_weights <- .nn_initialize_weights(n_features, size)
    w_vec <- unlist(initial_weights)

    # Objective: Negative Cox partial log-likelihood + L2 penalty
    obj_fn <- function(w) {
        ws <- .nn_unpack_weights(w, n_features, size)
        res <- .nn_forward_pass(X_sorted, ws)
        log_risk <- res$output

        lr_exp <- exp(log_risk)
        risk_sum <- cumsum(lr_exp) # Cumulative risk set sum

        # Log-likelihood for events only
        ll <- sum(log_risk[event_sorted == 1] - log(risk_sum[event_sorted == 1] + 1e-10))
        pen <- decay * sum(w^2)
        -ll + pen
    }

    # Gradient: Efficient computation using reverse cumsum
    grad_fn <- function(w) {
        ws <- .nn_unpack_weights(w, n_features, size)
        pass <- .nn_forward_pass(X_sorted, ws)
        log_risk <- pass$output
        A1 <- pass$hidden_activations

        lr_exp <- exp(log_risk)
        risk_sum <- cumsum(lr_exp)

        # Gradient w.r.t. log-risk (eta)
        term <- event_sorted / (risk_sum + 1e-10)
        weighted_sum <- rev(cumsum(rev(term))) # Sum over future events in risk set
        d_eta <- -event_sorted + weighted_sum * lr_exp

        # Backpropagation through network
        dZ2 <- d_eta
        dW2 <- t(A1) %*% dZ2
        db2 <- sum(dZ2)
        dA1 <- dZ2 %*% t(ws$W2)
        dZ1 <- dA1 * (A1 * (1 - A1)) # Sigmoid derivative
        dW1 <- t(X_sorted) %*% dZ1
        db1 <- colSums(dZ1)

        grad <- c(as.vector(dW1), db1, as.vector(dW2), db2)
        grad + 2 * decay * w
    }

    # Optimize
    opt <- stats::optim(w_vec, obj_fn, grad_fn,
        method = "BFGS",
        control = list(maxit = maxit)
    )

    final_ws <- .nn_unpack_weights(opt$par, n_features, size)

    list(
        weights = final_ws,
        convergence = opt$convergence,
        value = opt$value,
        iterations = opt$counts,
        sorted_data = list(X = X_sorted, time = time_sorted, event = event_sorted)
    )
}

# ============================================================================
# Baseline Hazard Estimation
# ============================================================================

#' Compute Breslow baseline hazard (vectorized)
#'
#' Efficiently computes the cumulative baseline hazard using vectorized
#' operations instead of explicit loops. This is O(n + k) instead of O(n*k).
#'
#' @param risks Numeric vector; exp(log-risk) from forward pass (n x 1)
#' @param time_sorted Numeric vector; times in DESCENDING order
#' @param event_sorted Numeric vector; binary events corresponding to sorted times
#' @return Data frame with columns 'time' and 'hazard' (cumulative)
#' @keywords internal
.nn_compute_baseline_hazard_vectorized <- function(risks, time_sorted, event_sorted) {
    # Get unique event times (ascending for output)
    event_idx <- which(event_sorted == 1)
    if (length(event_idx) == 0) {
        return(data.frame(time = numeric(0), hazard = numeric(0)))
    }

    event_times <- sort(unique(time_sorted[event_idx]))
    n_times <- length(event_times)
    base_haz <- numeric(n_times)

    # Vectorized computation
    for (i in seq_along(event_times)) {
        t <- event_times[i]
        # Risk set: all observations with time >= t
        risk_set_idx <- time_sorted >= t
        denom <- sum(risks[risk_set_idx])

        # Numerator: count events at exactly time t
        n_events <- sum(time_sorted == t & event_sorted == 1)
        base_haz[i] <- n_events / (denom + 1e-10) # Avoid division by zero
    }

    cum_base_haz <- cumsum(base_haz)
    data.frame(time = event_times, hazard = cum_base_haz)
}

# ============================================================================
# Data Preparation
# ============================================================================

#' Prepare data for neural network training
#'
#' Handles factor encoding, scaling, and model matrix construction
#' with proper bookkeeping for prediction-time alignment.
#'
#' @param data Data frame; raw training data
#' @param expvars Character vector; feature column names
#' @return List with model matrix, scaling params, and factor levels
#' @keywords internal
.nn_prepare_data <- function(data, expvars) {
    # Extract factor levels
    factor_levels <- list()
    for (v in expvars) {
        if (is.factor(data[[v]]) || is.character(data[[v]])) {
            factor_levels[[v]] <- levels(as.factor(data[[v]]))
        }
    }

    # Identify numeric variables for scaling
    numeric_vars <- expvars[sapply(data[, expvars, drop = FALSE], is.numeric)]
    scaling_params <- list(
        mean = sapply(data[, numeric_vars, drop = FALSE], mean, na.rm = TRUE),
        sd = sapply(data[, numeric_vars, drop = FALSE], sd, na.rm = TRUE)
    )

    # Scale numeric variables
    data_scaled <- data
    if (length(numeric_vars) > 0) {
        data_scaled[, numeric_vars] <- scale(data[, numeric_vars, drop = FALSE])
    }

    # Create model matrix
    x_train <- stats::model.matrix(~ . - 1, data = data_scaled[, expvars, drop = FALSE])

    list(
        X = x_train,
        scaling_params = scaling_params,
        factor_levels = factor_levels,
        n_features = ncol(x_train)
    )
}

#' Apply training-time transformations to new data
#'
#' Ensures factor levels and scaling match training data.
#'
#' @param newdata Data frame; new observations
#' @param expvars Character vector; feature column names
#' @param scaling_params List; mean and sd from training
#' @param factor_levels List; factor levels from training
#' @return Numeric matrix; transformed design matrix
#' @keywords internal
.nn_transform_newdata <- function(newdata, expvars, scaling_params, factor_levels) {
    # Align factors
    for (v in expvars) {
        if (v %in% names(factor_levels)) {
            newdata[[v]] <- factor(newdata[[v]], levels = factor_levels[[v]])
        }
    }

    # Scale numerics
    num_vars <- names(scaling_params$mean)
    if (length(num_vars) > 0) {
        newdata[, num_vars] <- scale(
            newdata[, num_vars, drop = FALSE],
            center = scaling_params$mean,
            scale = scaling_params$sd
        )
    }

    # Create model matrix
    stats::model.matrix(~ . - 1, data = newdata[, expvars, drop = FALSE])
}

# ============================================================================
# Prediction Helpers
# ============================================================================

#' Predict survival probabilities from fitted NN model
#'
#' Combines neural network risk scores with baseline hazard to compute
#' survival probabilities at requested time points.
#'
#' @param model List; fitted model object with weights, baseline_hazard, etc.
#' @param newdata Data frame; new observations
#' @param times Numeric vector; time points for prediction
#' @param task Task object; contains metadata (id_col, etc.)
#' @param expvars Character vector; feature names
#' @param time_grid Numeric vector; default time grid if times=NULL
#' @param set Character; prediction set label ("train", "test", etc.)
#' @return Tidy data frame with survival predictions
#' @keywords internal
.nn_predict_survival <- function(model, newdata, times, task, expvars, time_grid, set = "test") {
    # Transform new data
    X_new <- .nn_transform_newdata(newdata, expvars, model$scaling_params, model$factor_levels)

    # Compute risk scores
    res <- .nn_forward_pass(X_new, model$weights)
    risks <- exp(res$output)

    # Interpolate baseline hazard
    base_times <- c(0, model$baseline_hazard$time)
    base_H0 <- c(0, model$baseline_hazard$hazard)
    target_times <- if (is.null(times)) time_grid else sort(unique(as.numeric(times)))
    H0_interp <- stats::approx(base_times, base_H0, xout = target_times, rule = 2)$y

    # Compute survival: S(t|x) = exp(-H0(t) * exp(beta'x))
    Lambda <- outer(H0_interp, as.vector(risks))
    Surv <- exp(-Lambda)

    # Format as tidy data frame
    id_col <- newdata[[task$id_col]]
    flat_surv <- as.vector(Surv)

    new_survival_prediction(
        id = rep(id_col, each = length(target_times)),
        time = rep(target_times, times = length(id_col)),
        surv = flat_surv,
        model = rep("shallownn", length(flat_surv)),
        ensemble = FALSE,
        set = set
    )
}

#' Predict cumulative incidence functions for competing risks
#'
#' Uses Aalen-Johansen estimator to combine cause-specific hazards
#' into cumulative incidence functions.
#'
#' @param models_list Named list; one fitted model per cause
#' @param newdata Data frame; new observations
#' @param times Numeric vector; time points for prediction
#' @param task Task object; contains metadata
#' @param expvars Character vector; feature names
#' @param time_grid Numeric vector; default time grid
#' @param set Character; prediction set label
#' @return Tidy data frame with CIF predictions
#' @keywords internal
.nn_predict_cif <- function(models_list, newdata, times, task, expvars, time_grid, set = "test") {
    # Extract scaling params from first model (same across causes)
    first_model <- models_list[[1]]
    X_new <- .nn_transform_newdata(newdata, expvars, first_model$scaling_params, first_model$factor_levels)

    n_obs <- nrow(X_new)
    target_times <- if (is.null(times)) time_grid else sort(unique(as.numeric(times)))
    n_times <- length(target_times)
    cause_labels <- names(models_list)

    # Compute cumulative hazards for each cause
    cum_hazards_array <- array(0, dim = c(n_times, n_obs, length(cause_labels)))
    dimnames(cum_hazards_array)[[3]] <- cause_labels

    for (cause in cause_labels) {
        m <- models_list[[cause]]
        res <- .nn_forward_pass(X_new, m$weights)
        risks <- exp(res$output)

        # Interpolate baseline
        bt <- c(0, m$baseline_hazard$time)
        bh <- c(0, m$baseline_hazard$hazard)
        bh_interp <- stats::approx(bt, bh, xout = target_times, rule = 2)$y

        # Cause-specific cumulative hazard
        Lambda <- outer(bh_interp, as.vector(risks))
        cum_hazards_array[, , cause] <- Lambda
    }

    # Total hazard and overall survival
    total_Lambda <- apply(cum_hazards_array, c(1, 2), sum)
    S_overall <- exp(-total_Lambda)

    # Aalen-Johansen: CIF_k(t) = integral_0^t S(u-) dLambda_k(u)
    preds_list <- list()
    id_col <- newdata[[task$id_col]]

    for (cause in cause_labels) {
        cif_mat <- matrix(0, nrow = n_times, ncol = n_obs)
        L_k <- cum_hazards_array[, , cause]

        # Discrete approximation
        cif_mat[1, ] <- L_k[1, ] # At t=1, S(0-) ≈ 1
        if (n_times > 1) {
            for (t in 2:n_times) {
                dL <- pmax(L_k[t, ] - L_k[t - 1, ], 0)
                cif_mat[t, ] <- cif_mat[t - 1, ] + S_overall[t - 1, ] * dL
            }
        }

        flat_cif <- as.vector(cif_mat)
        preds_list[[cause]] <- new_cif_prediction(
            id = rep(id_col, each = n_times),
            time = rep(target_times, times = n_obs),
            cause = rep(cause, length(flat_cif)),
            cif = flat_cif,
            model = rep("shallownn", length(flat_cif)),
            ensemble = FALSE,
            set = set
        )
    }

    dplyr::bind_rows(preds_list)
}
