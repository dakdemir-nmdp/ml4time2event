# --- Internal Helper Functions for Custom NNet Training -#' @title (Helper) Profile Variables
#'
#' @description A placeholder function to profile variables. In a real scenario, this would
#' analyze variable types (numeric, factor), etc.
#' @param data A data frame.
#' @param expvars A character vector of variable names.
#' @return A list summarizing variable types.
#' @keywords internal
VariableProfile <- function(data, expvars) {
  varprofile <- vector(mode="list", length=length(expvars))
  names(varprofile) <- expvars
  for (vari in expvars) {
    if (vari %in% colnames(data)) {
      col_data <- data[[vari]]
      if (is.factor(col_data)) {
        varprofile[[vari]] <- table(col_data, useNA = "ifany")
      } else if (is.numeric(col_data)) {
        varprofile[[vari]] <- c(min = min(col_data, na.rm = TRUE), max = max(col_data, na.rm = TRUE))
      } else if (is.character(col_data)) {
        varprofile[[vari]] <- table(col_data, useNA = "ifany")
      } else {
        varprofile[[vari]] <- paste("Unsupported type:", class(col_data))
      }
    } else {
      varprofile[[vari]] <- "Variable not found in data"
    }
  }
  varprofile
}

#' @title Initialize Neural Network Weights
#' @description Initializes weights and biases for a single-hidden-layer network.
#' @param n_in Number of input features.
#' @param n_hidden Number of hidden units.
#' @return A list of weight matrices and bias vectors.
#' @keywords internal
initialize_weights <- function(n_in, n_hidden) {
  # Use a common initialization scheme (e.g., random uniform)
  # The nnet package uses a range of [-0.7, 0.7]
  rand_range <- 0.7
  W1 <- matrix(runif(n_in * n_hidden, -rand_range, rand_range), nrow = n_in, ncol = n_hidden)
  b1 <- runif(n_hidden, -rand_range, rand_range)
  W2 <- matrix(runif(n_hidden, -rand_range, rand_range), nrow = n_hidden, ncol = 1)
  b2 <- runif(1, -rand_range, rand_range)

  list(W1 = W1, b1 = b1, W2 = W2, b2 = b2)
}

#' @title Forward Pass
#' @description Computes the output of the neural network.
#' @param X Input matrix.
#' @param weights A list of weights and biases.
#' @return A list containing the final output (log-risk) and hidden layer activations.
#' @keywords internal
forward_pass <- function(X, weights) {
  # Sigmoid activation function (as in nnet package)
  sigmoid <- function(z) 1 / (1 + exp(-z))

  # Z1 = X %*% W1 + b1 (broadcast b1)
  Z1 <- sweep(X %*% weights$W1, 2, weights$b1, "+")
  A1 <- sigmoid(Z1)
  # Z2 = A1 %*% W2 + b2 (linear output)
  Z2 <- sweep(A1 %*% weights$W2, 2, weights$b2, "+")

  list(output = Z2, hidden_activations = A1)
}

#' @title Unpack Weights
#' @description Unpacks a weight vector into weight matrices and bias vectors.
#' @param w_vec Weight vector.
#' @param n_in Number of input features.
#' @param n_hidden Number of hidden units.
#' @return A list of weight matrices and bias vectors.
#' @keywords internal
unpack_weights <- function(w_vec, n_in, n_hidden) {
  W1_end <- n_in * n_hidden
  b1_end <- W1_end + n_hidden
  W2_end <- b1_end + n_hidden
  b2_end <- W2_end + 1

  list(
    W1 = matrix(w_vec[1:W1_end], nrow = n_in, ncol = n_hidden),
    b1 = w_vec[(W1_end + 1):b1_end],
    W2 = matrix(w_vec[(b1_end + 1):W2_end], nrow = n_hidden, ncol = 1),
    b2 = w_vec[b2_end]
  )
}

#' @title CRModel_DeepSurv
#'
#' @description Fit a DeepSurv neural network model for competing risks outcomes using Fine-Gray style loss.
#'
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (coded 0=censored, 1=cause1, 2=cause2, etc.)
#' @param failcode integer, the code for the event of interest (default: 1)
#' @param size integer, number of units in the hidden layer (default: 5)
#' @param decay numeric, L2 regularization parameter (default: 0.01)
#' @param maxit integer, maximum iterations for optimization (default: 1000)
#' @param verbose logical, print progress messages (default: FALSE)
#'
#' @return a list with the following components:
#'   \item{model}{fitted neural network model (list with weights)}
#'   \item{times}{unique event times from training data for the event of interest}
#'   \item{varprof}{variable profile list}
#'   \item{expvars}{character vector of explanatory variables}
#'   \item{factor_levels}{list of factor levels for categorical variables}
#'   \item{failcode}{the event code for the outcome of interest}
#'   \item{model_type}{character string "cr_deepsurv"}
#'
#' @importFrom stats model.matrix as.formula
#' @export
CRModel_DeepSurv <- function(data, expvars, timevar, eventvar, failcode = 1,
                            size = 5, decay = 0.01, maxit = 1000, verbose = FALSE) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.character(expvars) || length(expvars) == 0) {
    stop("'expvars' must be a non-empty character vector")
  }
  if (!timevar %in% colnames(data)) {
    stop("'timevar' not found in data: ", timevar)
  }
  if (!eventvar %in% colnames(data)) {
    stop("'eventvar' not found in data: ", eventvar)
  }
  missing_vars <- setdiff(expvars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following expvars not found in data: ", paste(missing_vars, collapse=", "))
  }
  if (!is.numeric(failcode) || length(failcode) != 1 || failcode < 1) {
    stop("'failcode' must be a positive integer")
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================
  if (verbose) cat("Creating variable profile...\n")
  varprof <- VariableProfile(data, expvars)

  # Store factor levels for prediction
  factor_levels <- list()
  for (var in expvars) {
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      factor_levels[[var]] <- levels(as.factor(data[[var]]))
    }
  }

  # Ensure event variable is numeric
  data[[eventvar]] <- as.numeric(data[[eventvar]])

  # Prepare data subset with complete cases
  XYTrain <- data[, c(timevar, eventvar, expvars), drop = FALSE]
  n_original <- nrow(XYTrain)
  XYTrain <- XYTrain[complete.cases(XYTrain), , drop = FALSE]
  n_complete <- nrow(XYTrain)

  if (n_complete < n_original) {
    warning(sprintf("Removed %d rows with missing values. %d rows remaining.",
                    n_original - n_complete, n_complete))
  }
  if (n_complete < 10) {
    stop("Insufficient data after removing missing values. Need at least 10 observations.")
  }

  # Get unique event times for the event of interest
  event_times <- XYTrain[[timevar]][XYTrain[[eventvar]] == failcode]
  if (length(event_times) == 0) {
    stop("No events of type ", failcode, " in training data. Cannot fit competing risks model.")
  }
  times <- sort(unique(event_times))

  # ============================================================================
  # Model Fitting (Custom DeepSurv Implementation for Competing Risks)
  # ============================================================================
  if (verbose) cat("Fitting DeepSurv neural network for competing risks (event type", failcode, ")...\n")

  # Standardize numeric variables
  numeric_vars <- expvars[sapply(XYTrain[, expvars, drop=FALSE], is.numeric)]
  scaling_params <- list(
    mean = sapply(XYTrain[, numeric_vars, drop=FALSE], mean, na.rm=TRUE),
    sd = sapply(XYTrain[, numeric_vars, drop=FALSE], sd, na.rm=TRUE)
  )

  # Apply scaling
  XYTrain_scaled <- XYTrain
  XYTrain_scaled[, numeric_vars] <- scale(XYTrain[, numeric_vars, drop=FALSE])

  # Create model matrix
  x_train <- stats::model.matrix(~ . - 1, data = XYTrain_scaled[, expvars, drop = FALSE])
  n_features <- ncol(x_train)

  # Sort data by descending time for efficient loss calculation
  sort_order <- order(-XYTrain[[timevar]], -XYTrain[[eventvar]])
  x_train <- x_train[sort_order, ]
  time_sorted <- XYTrain[[timevar]][sort_order]
  event_sorted <- XYTrain[[eventvar]][sort_order]

  # Create cause-specific indicators for Fine-Gray style loss
  # For competing risks, we use a modified loss that accounts for censoring and competing events
  status_fg <- ifelse(event_sorted == 0, 0,  # censored
                     ifelse(event_sorted == failcode, 1, 2))  # event of interest or competing

  # ============================================================================
  # Neural Network Training with Fine-Gray Style Loss
  # ============================================================================
  # Initialize weights
  initial_w_vec <- unlist(initialize_weights(n_features, size))

  # Define objective function (negative log-likelihood for Fine-Gray + L2 penalty)
  objective_fn <- function(w_vec) {
    weights <- unpack_weights(w_vec, n_features, size)

    # Forward pass to get log-risk (subdistribution hazard)
    log_risk <- forward_pass(x_train, weights)$output

    # Fine-Gray style loss calculation for competing risks
    # This implements a simplified version of the Fine-Gray likelihood
    risk_scores <- exp(log_risk)

    # Calculate loss similar to Fine-Gray model
    # For each event time, compute the contribution
    loss <- 0
    n_events <- 0

    for (i in seq_along(time_sorted)) {
      if (status_fg[i] == 1) {  # Event of interest
        n_events <- n_events + 1

        # Risk set at this time (all subjects still at risk)
        risk_set <- risk_scores[i:length(risk_scores)]

        # Contribution to log-likelihood
        # log(risk_i / sum(risk_set))
        if (sum(risk_set) > 0) {
          loss <- loss + log_risk[i] - log(sum(risk_set))
        }
      }
    }

    # For censored and competing events, they contribute 0 to the likelihood
    # (they're handled implicitly through the risk set calculations)

    # L2 penalty
    l2_penalty <- decay * sum(w_vec^2)

    # Return negative log-likelihood + penalty
    -loss + l2_penalty
  }

  # Define gradient function
  gradient_fn <- function(w_vec) {
    weights <- unpack_weights(w_vec, n_features, size)

    # Forward pass
    pass <- forward_pass(x_train, weights)
    log_risk <- pass$output
    A1 <- pass$hidden_activations

    # Backward pass for Fine-Gray style loss
    risk_scores <- exp(log_risk)

    d_log_risk <- rep(0, length(log_risk))

    for (i in seq_along(time_sorted)) {
      if (status_fg[i] == 1) {  # Event of interest
        # Risk set at this time
        risk_set <- risk_scores[i:length(risk_scores)]
        risk_set_sum <- sum(risk_set)

        if (risk_set_sum > 0) {
          # Gradient contribution
          d_log_risk[i] <- d_log_risk[i] + 1  # From log_risk[i] term
          d_log_risk[i:length(log_risk)] <- d_log_risk[i:length(log_risk)] -
            risk_scores[i:length(log_risk)] / risk_set_sum  # From -log(sum) term
        }
      }
    }

    # Backpropagate through output layer
    dZ2 <- d_log_risk
    dW2 <- t(A1) %*% dZ2
    db2 <- sum(dZ2)
    dA1 <- dZ2 %*% t(weights$W2)

    # Backpropagate through sigmoid
    dZ1 <- dA1 * (A1 * (1 - A1))
    dW1 <- t(x_train) %*% dZ1
    db1 <- colSums(dZ1)

    # Add L2 gradient
    l2_grad <- 2 * decay * w_vec

    # Pack gradients
    c(as.vector(dW1), db1, as.vector(dW2), db2) + l2_grad
  }

  # Optimize
  if (verbose) cat("Optimizing neural network parameters...\n")
  opt_result <- optim(
    par = initial_w_vec,
    fn = objective_fn,
    gr = gradient_fn,
    method = "BFGS",
    control = list(maxit = maxit, trace = as.numeric(verbose))
  )

  final_weights <- unpack_weights(opt_result$par, n_features, size)

  # ============================================================================
  # Calculate Baseline Subdistribution Hazard (Fine-Gray style)
  # ============================================================================
  final_log_risk <- forward_pass(x_train, final_weights)$output
  final_risk_scores <- exp(final_log_risk)

  # Calculate baseline subdistribution hazard similar to Fine-Gray
  baseline_subhaz <- sapply(times, function(t_i) {
    # Find events at this time
    event_indices <- which(time_sorted == t_i & status_fg == 1)
    if (length(event_indices) > 0) {
      # Risk set: all subjects still at risk at time t_i
      risk_set_start <- min(event_indices)
      risk_set <- final_risk_scores[risk_set_start:length(final_risk_scores)]
      sum(1 / sum(risk_set))  # Simplified baseline hazard calculation
    } else {
      0
    }
  })

  cumulative_baseline_subhaz <- cumsum(baseline_subhaz)

  # ============================================================================
  # Return Results
  # ============================================================================
  if (verbose) cat("DeepSurv competing risks model fitting complete.\n")

  # Create a model object that mimics nnet structure for compatibility
  model <- list(
    weights = final_weights,
    baseline_subhaz = data.frame(time = times, subhaz = cumulative_baseline_subhaz),
    scaling_params = scaling_params,
    n = c(n_features, size, 1),  # Mimic nnet structure
    call = match.call()
  )
  class(model) <- "nnet"

  result <- list(
    model = model,
    times = times,
    varprof = varprof,
    expvars = expvars,
    factor_levels = factor_levels,
    failcode = failcode,
    model_type = "cr_deepsurv"
  )

  class(result) <- "ml4t2e_cr_deepsurv"
  return(result)
}

#' @title Predict_CRModel_DeepSurv
#'
#' @description Get predictions from a fitted DeepSurv competing risks model for new data.
#'
#' @param modelout the output from 'CRModel_DeepSurv'
#' @param newdata data frame with new observations for prediction
#' @param newtimes optional numeric vector of time points for prediction.
#'   If NULL (default), uses the baseline hazard times from training.
#'
#' @return a list containing:
#'   \item{CIFs}{predicted cumulative incidence function matrix
#'     (rows=times, cols=observations)}
#'   \item{Times}{the times at which CIFs are calculated}
#'
#' @importFrom stats model.matrix
#' @export
Predict_CRModel_DeepSurv <- function(modelout, newdata, newtimes = NULL) {

  # ============================================================================
  # Input Validation
  # ============================================================================
  if (missing(modelout)) {
    stop("argument \"modelout\" is missing")
  }
  if (missing(newdata)) {
    stop("argument \"newdata\" is missing")
  }
  if (!inherits(modelout, "ml4t2e_cr_deepsurv")) {
    stop("'modelout' must be output from CRModel_DeepSurv")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame")
  }

  # Check that required variables are present in newdata
  missing_vars <- setdiff(modelout$expvars, colnames(newdata))
  if (length(missing_vars) > 0) {
    stop("The following variables missing in newdata: ",
         paste(missing_vars, collapse = ", "))
  }

  # Validate newtimes if provided
  if (!is.null(newtimes)) {
    if (!is.numeric(newtimes)) {
      stop("'newtimes' must be numeric")
    }
    if (any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
  }

  # ============================================================================
  # Prepare newdata
  # ============================================================================
  # Ensure factor levels match training data
  newdata_prepared <- newdata[, modelout$expvars, drop = FALSE]

  for (var in modelout$expvars) {
    if (var %in% names(modelout$factor_levels)) {
      training_levels <- modelout$factor_levels[[var]]
      if (is.factor(newdata_prepared[[var]])) {
        newdata_prepared[[var]] <- factor(newdata_prepared[[var]], levels = training_levels)
      } else if (is.character(newdata_prepared[[var]])) {
        newdata_prepared[[var]] <- factor(newdata_prepared[[var]], levels = training_levels)
      }
    }
  }

  # Apply scaling to numeric variables
  numeric_vars <- names(modelout$model$scaling_params$mean)
  newdata_prepared[, numeric_vars] <- sweep(newdata_prepared[, numeric_vars, drop=FALSE],
                                           2, modelout$model$scaling_params$mean, "-")
  newdata_prepared[, numeric_vars] <- sweep(newdata_prepared[, numeric_vars, drop=FALSE],
                                           2, modelout$model$scaling_params$sd, "/")

  # Create model matrix
  x_new <- stats::model.matrix(~ . - 1, data = newdata_prepared)

  # ============================================================================
  # Make Predictions
  # ============================================================================
  # Get risk scores from neural network
  log_risk_scores <- forward_pass(x_new, modelout$model$weights)$output

  # Determine prediction times
  if (is.null(newtimes)) {
    pred_times <- c(0, modelout$model$baseline_subhaz$time)
    baseline_times <- pred_times
    baseline_cumsubhaz <- c(0, modelout$model$baseline_subhaz$subhaz)
  } else {
    pred_times <- sort(unique(newtimes))
    baseline_times <- c(0, modelout$model$baseline_subhaz$time)
    baseline_cumsubhaz <- c(0, modelout$model$baseline_subhaz$subhaz)
  }

  # Calculate CIF using Fine-Gray style: CIF(t|x) = 1 - exp(-H0(t) * exp(risk_score))
  # where H0(t) is the cumulative baseline subdistribution hazard
  cif_matrix <- 1 - exp(-outer(baseline_cumsubhaz, exp(as.vector(log_risk_scores))))

  # Ensure CIFs are properly bounded [0,1]
  cif_matrix <- pmin(pmax(cif_matrix, 0), 1)

  # ============================================================================
  # Apply Interpolation
  # ============================================================================
  if (is.null(newtimes)) {
    # Return predictions in native time grid: [observations, times]
    result_cifs <- t(cif_matrix)  # cif_matrix is [times, observations], transpose to [observations, times]
    result_times <- baseline_times
  } else {
    # Interpolate to new time points
    if (!is.numeric(newtimes) || any(newtimes < 0)) {
      stop("'newtimes' must be a numeric vector of non-negative values")
    }
    newtimes <- sort(unique(newtimes))

    # Use the standard CIF interpolation utility function
    pred_cifs <- cifMatInterpolaltor(
      probsMat = t(cif_matrix),  # cifMatInterpolaltor expects [observations, times]
      times = baseline_times,
      newtimes = newtimes
    )

    # cifMatInterpolaltor returns [newtimes, observations], transpose to [observations, newtimes]
    result_cifs <- t(pred_cifs)
    result_times <- newtimes
  }

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    CIFs = result_cifs,
    Times = result_times
  )

  return(result)
}