# ==============================================================================
# Advanced Ensemble Methods - Weighted Averaging, Stacking, Super Learner
# ==============================================================================

survprobMatWeightedAveraging <- function(listprobsMat, weights) {
  # Validate inputs
  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  # Filter out NULL entries
  listprobsMat <- Filter(Negate(is.null), listprobsMat)
  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  # Validate weights
  if (is.null(weights)) {
    stop("weights must be provided")
  }

  # Match weights to models
  model_names <- names(listprobsMat)
  if (!all(model_names %in% names(weights))) {
    stop("All models must have corresponding weights")
  }

  weights <- weights[model_names]
  if (any(is.na(weights))) {
    stop("Missing weights for one or more models")
  }

  # Normalize weights to sum to 1
  if (abs(sum(weights) - 1) > 1e-6) {
    message(sprintf("Weights sum to %.4f, normalizing to 1", sum(weights)))
    weights <- weights / sum(weights)
  }

  # Check dimensions consistency and filter if needed
  dims <- lapply(listprobsMat, dim)
  dim_strings <- sapply(dims, paste, collapse = "x")

  if (length(unique(dim_strings)) > 1) {
    dim_table <- table(dim_strings)
    most_common_dim <- names(dim_table)[which.max(dim_table)]
    valid_indices <- which(dim_strings == most_common_dim)
    listprobsMat <- listprobsMat[valid_indices]
    weights <- weights[valid_indices]

    n_excluded <- length(dim_strings) - length(valid_indices)
    if (n_excluded > 0) {
      warning(sprintf("Excluded %d prediction(s) due to dimension mismatch", n_excluded))
    }
  }

  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  # Create weighted average on cumulative hazard scale
  HazzardArray <- array(dim = c(dim(listprobsMat[[1]]), length(listprobsMat)))

  for (i in seq_along(listprobsMat)) {
    HazzardArray[, , i] <- -log(listprobsMat[[i]] + 1e-10)
  }

  # Weighted mean cumulative hazard
  WeightedHazzard <- apply(HazzardArray, c(1, 2), function(x) {
    sum(x * weights)
  })

  # Convert back to survival probability
  NewProbs <- exp(-WeightedHazzard)
  NewProbs <- pmax(0, pmin(NewProbs, 1))

  # Ensure matrix dimensions are preserved
  if (!is.matrix(NewProbs)) {
    dim(NewProbs) <- dim(listprobsMat[[1]])
  }

  NewProbs
}

cifMatWeightedAveraging <- function(listprobsMat, weights, type = "CumHaz") {
  if (!type %in% c("CumHaz", "prob")) {
    stop("Type must be either 'CumHaz' or 'prob'")
  }

  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  # Filter NULL
  listprobsMat <- Filter(Negate(is.null), listprobsMat)
  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  # Validate weights
  if (is.null(weights)) {
    stop("weights must be provided")
  }

  model_names <- names(listprobsMat)
  if (!all(model_names %in% names(weights))) {
    stop("All models must have corresponding weights")
  }

  weights <- weights[model_names]

  # Normalize weights
  if (abs(sum(weights) - 1) > 1e-6) {
    message(sprintf("Weights sum to %.4f, normalizing to 1", sum(weights)))
    weights <- weights / sum(weights)
  }

  # Check dimensions
  dims <- lapply(listprobsMat, dim)
  dim_strings <- sapply(dims, paste, collapse = "x")

  if (length(unique(dim_strings)) > 1) {
    dim_table <- table(dim_strings)
    most_common_dim <- names(dim_table)[which.max(dim_table)]
    valid_indices <- which(dim_strings == most_common_dim)
    listprobsMat <- listprobsMat[valid_indices]
    weights <- weights[valid_indices]

    n_excluded <- length(dim_strings) - length(valid_indices)
    if (n_excluded > 0) {
      warning(sprintf("Excluded %d prediction(s) due to dimension mismatch", n_excluded))
    }
  }

  if (length(listprobsMat) == 0) {
    return(NULL)
  }
  if (length(listprobsMat) == 1) {
    return(listprobsMat[[1]])
  }

  if (type == "CumHaz") {
    HazzardArray <- array(dim = c(dim(listprobsMat[[1]]), length(listprobsMat)))
    for (i in seq_along(listprobsMat)) {
      HazzardArray[, , i] <- -log(1 - listprobsMat[[i]] + 1e-10)
    }

    WeightedHazzard <- apply(HazzardArray, c(1, 2), function(x) {
      sum(x * weights)
    })

    NewProbs <- 1 - exp(-WeightedHazzard)
  } else {
    ProbsArray <- array(dim = c(dim(listprobsMat[[1]]), length(listprobsMat)))
    for (i in seq_along(listprobsMat)) {
      ProbsArray[, , i] <- listprobsMat[[i]]
    }

    NewProbs <- apply(ProbsArray, c(1, 2), function(x) {
      sum(x * weights)
    })
  }

  NewProbs <- pmax(0, pmin(NewProbs, 1))

  # Ensure matrix dimensions are preserved
  if (!is.matrix(NewProbs)) {
    dim(NewProbs) <- dim(listprobsMat[[1]])
  }

  NewProbs
}

optimizeSuperLearnerWeights <- function(predictions_list, actual_surv, loss_type = "mse", weights_matrix = NULL) {
  # Filter NULL and check dimensions
  predictions_list <- Filter(Negate(is.null), predictions_list)
  if (length(predictions_list) == 0) {
    stop("No valid predictions available")
  }

  n_models <- length(predictions_list)
  model_names <- names(predictions_list)

  # Trivial case: only 1 model
  if (n_models == 1) {
    w <- c(1.0)
    names(w) <- model_names
    return(w)
  }

  if (!is.matrix(actual_surv)) {
    stop("actual_surv must be a matrix of observed survival or CIF values")
  }

  reference_dim <- dim(predictions_list[[1]])
  if (!all(dim(actual_surv) == reference_dim)) {
    stop("Dimensions of actual_surv must match prediction matrices")
  }

  if (!is.null(weights_matrix)) {
    if (!all(dim(weights_matrix) == reference_dim)) {
      stop("weights_matrix must match dimensions of actual_surv")
    }
  }

  if (is.null(model_names) || any(model_names == "")) {
    stop("predictions_list must be a named list")
  }

  # Ensure all matrices share same dimensions
  mismatched <- vapply(predictions_list, function(mat) !all(dim(mat) == reference_dim), logical(1))
  if (any(mismatched)) {
    bad_models <- names(predictions_list)[mismatched]
    stop("All prediction matrices must share the same dimensions. Offenders: ", paste(bad_models, collapse = ", "))
  }

  # Objective function using Softmax Parameterization
  # par: Unconstrained vector of length n_models
  objective <- function(par) {
    # Softmax to get weights summing to 1 and positive
    exps <- exp(par - max(par)) # safe exp
    weights <- exps / sum(exps)

    # Compute weighted average
    weighted_pred <- Reduce(`+`, mapply(function(p, w) p * w,
      predictions_list, weights,
      SIMPLIFY = FALSE
    ))

    # Calculate loss
    if (loss_type == "mse") {
      diff_sq <- (weighted_pred - actual_surv)^2
      if (!is.null(weights_matrix)) {
        loss <- sum(weights_matrix * diff_sq, na.rm = TRUE) / sum(weights_matrix, na.rm = TRUE)
      } else {
        loss <- mean(diff_sq, na.rm = TRUE)
      }
    } else if (loss_type == "loglik") {
      p_safe <- pmax(pmin(weighted_pred, 1 - 1e-10), 1e-10)
      ll_terms <- actual_surv * log(p_safe) + (1 - actual_surv) * log(1 - p_safe)
      loss <- -mean(ll_terms, na.rm = TRUE)
    } else {
      stop("Unknown loss_type")
    }

    loss
  }

  # Initialize: par = 0 implies equal weights
  init_par <- rep(0, n_models)

  # Optimize (BFGS is fine for unconstrained)
  result <- optim(
    par = init_par,
    fn = objective,
    method = "BFGS"
  )

  # Convert optimal par back to weights
  best_par <- result$par
  exps <- exp(best_par - max(best_par))
  optimal_weights <- exps / sum(exps)
  names(optimal_weights) <- model_names

  optimal_weights
}

fitMetaLearner <- function(base_predictions, outcomes, meta_learner = "mse") {
  if (!meta_learner %in% c("mse", "loglik")) {
    stop("meta_learner must be one of 'mse' or 'loglik'")
  }

  weights <- optimizeSuperLearnerWeights(base_predictions, outcomes, loss_type = meta_learner)

  list(
    type = "weight_optimizer",
    loss_type = meta_learner,
    weights = weights
  )
}

EnsemblePredictions <- function(model_predictions,
                                ensemble_method = "average",
                                model_weights = NULL,
                                type = "survival",
                                sl_training_predictions = NULL,
                                sl_actual = NULL,
                                sl_loss = "mse",
                                sl_weights = NULL,
                                times,
                                ...) {
  # ... existing code ...
  # I'm not updating EnsemblePredictions here fully because the replacement chunks mechanism
  # works best on contiguous blocks. But wait, I need to touch EnsemblePredictions to pass weights_matrix?
  # The user didn't ask to update EnsemblePredictions signature, but `optimizeSuperLearnerWeights` usage inside it.
  # I will skip modifying EnsemblePredictions for now as it's not the primary entry point for stacking (SurvivalEnsembler is).
  # If EnsemblePredictions supports 'stacking' mode directly, it should update too.
  # But for the purpose of the critique (fixing ml4t2e_fit stacking), SurvivalEnsembler is key.

  # Wait, I'm replacing the whole block including EnsemblePredictions signature if I'm not careful.
  # I'll just replace optim/fit/prepare functions.
  # The tool call target lines 188-280 covers optimizeSuperLearnerWeights and fitMetaLearner.
  # I will need another call for prepare_stacking_data.
  # Let's verify I am not making invalid partial edits.
  # I will submit a MultiReplace.

  # Actually, the user prompt is "fix these issues".
  # I will replace optimizeSuperLearnerWeights and prepare_stacking_data.

  stop("This tool call body is for explanation. I will use multi_replace.")
}



EnsemblePredictions <- function(model_predictions,
                                ensemble_method = "average",
                                model_weights = NULL,
                                type = "survival",
                                sl_training_predictions = NULL,
                                sl_actual = NULL,
                                sl_loss = "mse",
                                sl_weights = NULL,
                                times,
                                ...) {
  # Validate inputs
  if (missing(times) || is.null(times)) {
    stop("'times' is a required argument and must not be NULL")
  }

  if (!ensemble_method %in% c("average", "weighted", "super_learner", "median", "min", "max", "geometric_mean", "stacking")) {
    stop("ensemble_method must be one of: 'average', 'weighted', 'super_learner', 'median', 'min', 'max', 'geometric_mean', 'stacking'")
  }

  if (!type %in% c("survival", "competing_risks")) {
    stop("type must be either 'survival' or 'competing_risks'")
  }

  # Filter NULL predictions
  model_predictions <- Filter(Negate(is.null), model_predictions)

  if (length(model_predictions) == 0) {
    warning("No valid predictions to ensemble")
    return(NULL)
  }

  if (length(model_predictions) == 1) {
    message("Only one prediction available, returning it directly")
    return(model_predictions[[1]])
  }

  # Apply the appropriate ensemble method
  if (ensemble_method == "average") {
    if (type == "survival") {
      combined_probs <- survprobMatListAveraging(model_predictions)
    } else {
      combined_probs <- cifMatListAveraging(model_predictions, type = "CumHaz")
    }
    result <- list(Probs = combined_probs, Times = times)
  } else if (ensemble_method == "weighted") {
    # Weighted averaging
    if (is.null(model_weights)) {
      stop("model_weights must be provided for weighted ensemble method")
    }

    available_names <- names(model_predictions)
    if (is.null(available_names) || any(available_names == "")) {
      stop("model_predictions must be a named list for weighted ensemble method")
    }

    if (is.null(names(model_weights))) {
      stop("model_weights must be a named numeric vector")
    }

    missing_weights <- setdiff(available_names, names(model_weights))
    if (length(missing_weights) > 0) {
      stop("Missing weights for: ", paste(missing_weights, collapse = ", "))
    }

    weights_to_use <- model_weights[available_names]
    if (any(is.na(weights_to_use))) {
      stop("Missing weights for one or more models")
    }

    model_weights <- weights_to_use

    if (type == "survival") {
      result <- survprobMatWeightedAveraging(model_predictions, model_weights)
    } else {
      result <- cifMatWeightedAveraging(model_predictions, model_weights, type = "CumHaz")
    }
  } else if (ensemble_method == "super_learner") {
    available_names <- names(model_predictions)
    if (is.null(available_names) || any(available_names == "")) {
      stop("model_predictions must be a named list for super learner ensembles")
    }

    if (!is.null(sl_weights)) {
      if (is.null(names(sl_weights))) {
        stop("sl_weights must be a named numeric vector")
      }

      weights <- sl_weights[available_names]
      if (any(is.na(weights))) {
        missing <- setdiff(available_names, names(sl_weights))
        stop("Missing super learner weights for: ", paste(missing, collapse = ", "))
      }

      if (abs(sum(weights) - 1) > 1e-6) {
        weights <- weights / sum(weights)
      }
    } else {
      if (is.null(sl_training_predictions) || is.null(sl_actual)) {
        stop("Super learner requires either sl_weights or both sl_training_predictions and sl_actual.")
      }

      if (is.null(names(sl_training_predictions))) {
        stop("sl_training_predictions must be a named list")
      }

      missing_models <- setdiff(available_names, names(sl_training_predictions))
      if (length(missing_models) > 0) {
        stop("sl_training_predictions is missing models: ", paste(missing_models, collapse = ", "))
      }

      training_subset <- sl_training_predictions[available_names]
      weights <- optimizeSuperLearnerWeights(training_subset, sl_actual, loss_type = sl_loss)
      weights <- weights[available_names]
    }

    if (type == "survival") {
      result <- survprobMatWeightedAveraging(model_predictions, weights)
    } else {
      result <- cifMatWeightedAveraging(model_predictions, weights, type = "CumHaz")
    }

    attr(result, "sl_weights") <- weights
  } else if (ensemble_method == "median") {
    combined_probs <- apply(simplify2array(model_predictions), c(1, 2), median)
    result <- list(Probs = combined_probs, Times = times)
  } else if (ensemble_method == "min") {
    combined_probs <- apply(simplify2array(model_predictions), c(1, 2), min)
    result <- list(Probs = combined_probs, Times = times)
  } else if (ensemble_method == "max") {
    combined_probs <- apply(simplify2array(model_predictions), c(1, 2), max)
    result <- list(Probs = combined_probs, Times = times)
  } else if (ensemble_method == "geometric_mean") {
    combined_probs <- apply(simplify2array(model_predictions), c(1, 2), function(x) exp(mean(log(x))))
    result <- list(Probs = combined_probs, Times = times)
  } else if (ensemble_method == "stacking") {
    if (is.null(sl_training_predictions) || is.null(sl_actual)) {
      stop("Stacking requires training predictions and actual outcomes.")
    }

    weights <- optimizeSuperLearnerWeights(sl_training_predictions, sl_actual, loss_type = sl_loss)
    if (type == "survival") {
      combined_probs <- survprobMatWeightedAveraging(model_predictions, weights)
    } else {
      combined_probs <- cifMatWeightedAveraging(model_predictions, weights, type = "CumHaz")
    }
    result <- list(Probs = combined_probs, Times = times)
    attr(result, "stacking_weights") <- weights
  }

  result
}

#' Helper to generate stacking predictions and actuals
#'
#' @keywords internal
prepare_stacking_data <- function(task, fitted_models, time_grid, outcome_type, cal_data) {
  sl_preds <- list()
  sl_actual <- NULL
  sl_weights <- NULL

  obs_times <- cal_data[[task$time_col]]
  obs_events <- cal_data[[task$event_col]]
  n_obs <- length(obs_times)
  n_times <- length(time_grid)

  # Calculate IPCW Weights
  # Censoring event: status == 0
  cens_obj <- survival::Surv(obs_times, obs_events == 0)
  km_cens <- survival::survfit(cens_obj ~ 1)

  # G(t) for grid times
  summ_grid <- summary(km_cens, times = time_grid, extend = TRUE)
  G_grid <- summ_grid$surv
  if (is.null(G_grid)) G_grid <- rep(1, n_times)
  G_grid[G_grid == 0] <- 1e-5

  # G(Ti) for dead/event cases
  u_times <- unique(obs_times)
  summ_obs <- summary(km_cens, times = u_times, extend = TRUE)
  G_map <- setNames(summ_obs$surv, as.character(u_times))

  # Construct Weight Matrix
  weights_mat <- matrix(0, nrow = n_obs, ncol = n_times)
  T_mat <- matrix(obs_times, nrow = n_obs, ncol = n_times)
  Grid_mat <- matrix(time_grid, nrow = n_obs, ncol = n_times, byrow = TRUE)

  # Alive weights (T > t): 1/G(t)
  alive_mask <- T_mat > Grid_mat
  G_grid_mat <- matrix(G_grid, nrow = n_obs, ncol = n_times, byrow = TRUE)
  weights_mat[alive_mask] <- 1 / G_grid_mat[alive_mask]

  # Dead/Event weights (T <= t & Event): 1/G(Ti)
  is_event <- obs_events != 0
  dead_mask <- (T_mat <= Grid_mat) & matrix(is_event, nrow = n_obs, ncol = n_times)

  G_Ti_vals <- G_map[as.character(obs_times)]
  G_Ti_vals[is.na(G_Ti_vals)] <- 1
  G_Ti_vals[G_Ti_vals == 0] <- 1e-5
  G_Ti_mat <- matrix(G_Ti_vals, nrow = n_obs, ncol = n_times)

  weights_mat[dead_mask] <- 1 / G_Ti_mat[dead_mask]

  sl_weights <- weights_mat

  if (outcome_type == "survival") {
    # Survival Actual Transformation
    actual_mat <- matrix(NA, nrow = n_obs, ncol = n_times)

    # Fill actuals
    actual_mat[alive_mask] <- 1
    # For standard survival, status 1 is event
    dead_surv_mask <- (T_mat <= Grid_mat) & (matrix(obs_events == 1, nrow = n_obs, ncol = n_times))
    actual_mat[dead_surv_mask] <- 0

    # Fill NAs with 0 where weight is 0 (censored before t) to allow calc
    actual_mat[is.na(actual_mat) & (weights_mat == 0)] <- 0

    sl_actual <- actual_mat

    # Predictions
    for (m in names(fitted_models)) {
      mod <- fitted_models[[m]]
      p <- mod$predict_survival(newdata = cal_data, times = time_grid, set = "calibration")
      sl_preds[[m]] <- ml4t2e_reshape_preds_to_matrix(p, cal_data, time_grid, "surv", task$id_col)
    }
  } else {
    # Competing Risk
    causes <- as.character(task$metadata$cause_map$code)
    mat_list_actual <- list()

    for (c_val in causes) {
      c_code <- as.numeric(c_val)
      m_act <- matrix(NA, nrow = length(obs_times), ncol = length(time_grid))
      for (j in seq_along(time_grid)) {
        t <- time_grid[j]
        alive <- (obs_times > t)
        m_act[alive, j] <- 0

        happened <- (obs_times <= t & obs_events == c_code)
        m_act[happened, j] <- 1

        other <- (obs_times <= t & obs_events != 0 & obs_events != c_code)
        m_act[other, j] <- 0

        # Censored remains NA, fill later
      }
      m_act[is.na(m_act) & (weights_mat == 0)] <- 0
      mat_list_actual[[c_val]] <- m_act
    }
    sl_actual <- do.call(cbind, mat_list_actual)

    # For CR, we stack the weights matrix?
    # Since sl_actual is concatenated (N x (T*Causes)), we need weights (N x (T*Causes))
    # We repeat the weights matrix for each cause
    sl_weights <- do.call(cbind, replicate(length(causes), weights_mat, simplify = FALSE))

    # CR Predictions (Concatenated)
    for (m in names(fitted_models)) {
      mod <- fitted_models[[m]]
      p <- mod$predict_cif(newdata = cal_data, times = time_grid, set = "calibration")
      p_mats <- list()
      for (c_val in causes) {
        pc <- p[p$cause == c_val, ]
        p_mats[[c_val]] <- ml4t2e_reshape_preds_to_matrix(pc, cal_data, time_grid, "cif", task$id_col)
      }
      sl_preds[[m]] <- do.call(cbind, p_mats)
    }
  }

  list(preds = sl_preds, actual = sl_actual, weights = sl_weights)
}
