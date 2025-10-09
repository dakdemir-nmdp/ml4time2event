# ==============================================================================
# Advanced Ensemble Methods - Weighted Averaging, Stacking, Super Learner
# ==============================================================================

#' @title Weighted Average for Survival Probability Matrices
#' @description Average survival probability matrices with user-specified weights
#' @param listprobsMat List of survival probability matrices (rows=times, cols=observations)
#' @param weights Named numeric vector of weights (must sum to 1)
#' @return Weighted average survival probability matrix
#' @noRd
survprobMatWeightedAveraging <- function(listprobsMat, weights) {
  # Validate inputs
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Filter out NULL entries
  listprobsMat <- Filter(Negate(is.null), listprobsMat)
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

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
  dim_strings <- sapply(dims, paste, collapse="x")

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

  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

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

#' @title Weighted Average for CIF Matrices
#' @description Average CIF matrices with user-specified weights
#' @param listprobsMat List of CIF matrices (rows=times, cols=observations)
#' @param weights Named numeric vector of weights (must sum to 1)
#' @param type Character, either "CumHaz" or "prob"
#' @return Weighted average CIF matrix
#' @noRd
cifMatWeightedAveraging <- function(listprobsMat, weights, type = "CumHaz") {
  if (!type %in% c("CumHaz", "prob")) {
    stop("Type must be either 'CumHaz' or 'prob'")
  }

  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

  # Filter NULL
  listprobsMat <- Filter(Negate(is.null), listprobsMat)
  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

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

  if (length(listprobsMat) == 0) return(NULL)
  if (length(listprobsMat) == 1) return(listprobsMat[[1]])

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

#' @title Optimize Super Learner Weights
#' @description Find optimal weights for ensemble using cross-validation loss
#' @param predictions_list List of prediction matrices from different models
#' @param actual_surv Actual survival matrix from training data
#' @param loss_type Type of loss function ("mse", "loglik")
#' @return Numeric vector of optimal weights
#' @importFrom stats optim
#' @noRd
optimizeSuperLearnerWeights <- function(predictions_list, actual_surv, loss_type = "mse") {
  n_models <- length(predictions_list)

  # Filter NULL and check dimensions
  predictions_list <- Filter(Negate(is.null), predictions_list)
  if (length(predictions_list) == 0) {
    stop("No valid predictions available")
  }

  if (!is.matrix(actual_surv)) {
    stop("actual_surv must be a matrix of observed survival or CIF values")
  }

  reference_dim <- dim(predictions_list[[1]])
  if (!all(dim(actual_surv) == reference_dim)) {
    stop("Dimensions of actual_surv must match prediction matrices")
  }

  model_names <- names(predictions_list)
  if (is.null(model_names) || any(model_names == "")) {
    stop("predictions_list must be a named list")
  }

  # Ensure all matrices share same dimensions
  mismatched <- vapply(predictions_list, function(mat) !all(dim(mat) == reference_dim), logical(1))
  if (any(mismatched)) {
    bad_models <- names(predictions_list)[mismatched]
    stop("All prediction matrices must share the same dimensions. Offenders: ", paste(bad_models, collapse = ", "))
  }

  # Objective function to minimize
  objective <- function(weights) {
    # Ensure weights sum to 1 and are non-negative
    if (any(weights < 0) || abs(sum(weights) - 1) > 1e-6) {
      return(.Machine$double.xmax)
    }

    # Compute weighted average
    weighted_pred <- Reduce(`+`, mapply(function(p, w) p * w,
                                         predictions_list, weights,
                                         SIMPLIFY = FALSE))

    # Calculate loss
    if (loss_type == "mse") {
      loss <- mean((weighted_pred - actual_surv)^2, na.rm = TRUE)
    } else if (loss_type == "loglik") {
      # Negative log-likelihood (for survival probabilities)
      loss <- -mean(log(pmax(weighted_pred, 1e-10)), na.rm = TRUE)
    } else {
      stop("Unknown loss_type")
    }

    loss
  }

  # Initialize with equal weights
  init_weights <- rep(1 / n_models, n_models)

  # Optimize with constraints
  result <- optim(
    par = init_weights,
    fn = objective,
    method = "L-BFGS-B",
    lower = rep(0, n_models),
    upper = rep(1, n_models)
  )

  # Normalize to sum to 1
  optimal_weights <- result$par / sum(result$par)
  names(optimal_weights) <- names(predictions_list)

  optimal_weights
}

#' @title Simple Meta-Learner for Stacking
#' @description Optimise ensemble weights for stacked survival predictions
#' @param base_predictions List of prediction matrices from individual base learners
#' @param outcomes Matrix of observed survival/CIF values aligned with the prediction matrices
#' @param meta_learner Loss type for optimisation ("mse" or "loglik")
#' @return List containing the optimised weights and loss specification
#' @noRd
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

#' @title Ensemble Predictions with Different Methods
#' @description Combine model predictions using specified ensemble method
#' @param model_predictions List of individual model predictions
#' @param ensemble_method Method to use: "average", "weighted", "super_learner"
#' @param model_weights Named numeric vector of weights (for weighted method)
#' @param type Type of predictions: "survival" or "competing_risks"
#' @param sl_training_predictions Optional list of training-set prediction matrices (for super learner optimisation)
#' @param sl_actual Optional matrix of observed survival/CIF values matching `sl_training_predictions`
#' @param sl_loss Loss function for super learner optimisation ("mse" or "loglik")
#' @param sl_weights Optional pre-computed super learner weights
#' @param times Optional vector of time points (required for some methods)
#' @param ... Additional arguments passed to specific methods
#' @return Combined ensemble predictions
#' @details When `ensemble_method = "super_learner"`, the returned matrix includes an
#'   attribute `sl_weights` containing the optimised weights used for the combination.
#' @export
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
    combined_probs <- survprobMatWeightedAveraging(model_predictions, weights)
    result <- list(Probs = combined_probs, Times = times)
    attr(result, "stacking_weights") <- weights
  }

  result
}
