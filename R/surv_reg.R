#' @title SurvModel_SurvReg
#'
#' @description Fit parametric survival models using survival::survreg with forward selection.
#' @param data data frame with explanatory and outcome variables
#' @param expvars character vector of names of explanatory variables in data
#' @param timevar character name of time variable in data
#' @param eventvar character name of event variable in data (needs to be 0/1)
#' @param dist distribution for the parametric model, one of "weibull", "exponential", "gaussian", "logistic","lognormal" or "loglogistic".
#'
#' @return a list of three items: survregOut: the final fitted survreg model object after forward selection,
#'  times: unique event times from the training data,
#'  varprof: profile of explanatory variables.
#'
#' @importFrom survival survreg Surv
#' @importFrom stats AIC as.formula predict quantile terms
#' @export
SurvModel_SurvReg <- function(data, expvars, timevar, eventvar, dist = "exponential") {
  # Assuming VariableProfile and survivalProbsInterpolator are loaded/available
  varprof <- VariableProfile(data, expvars) # Placeholder

  # Ensure event variable is numeric 0/1
  data[, eventvar] <- as.numeric(data[, eventvar] == 1)
  dataYX <- data[, c(timevar, eventvar, expvars), drop=FALSE]

  # Check if at least one variable is provided
  if (length(expvars) == 0) {
    stop("No explanatory variables provided (expvars). At least one variable is required for selection.")
  }

  # --- Forward Selection with AIC ---
  selected_vars <- c()
  candidate_vars <- expvars
  best_model <- NULL
  # Start with intercept-only model AIC
  null_formula <- stats::as.formula(paste("survival::Surv(", timevar, ",", eventvar, ") ~ 1"))
  null_model <- tryCatch(survival::survreg(null_formula, data = dataYX, dist = dist, x = FALSE, y = FALSE),
                         error = function(e) {warning("Failed to fit null model: ", e$message); NULL})
  if (is.null(null_model)) stop("Failed to fit intercept-only model.")
  best_aic <- stats::AIC(null_model)
  print(paste("Initial AIC (Intercept only):", round(best_aic, 2)))



  # Forward selection loop
  while (length(candidate_vars) > 0) {
    aic_values <- numeric(length(candidate_vars))
    names(aic_values) <- candidate_vars
    improved <- FALSE

    for (i in seq_along(candidate_vars)) {
      var <- candidate_vars[i]
      current_vars <- c(selected_vars, var)
      formula_str <- paste("survival::Surv(", timevar, ",", eventvar, ") ~", paste(current_vars, collapse = "+"))
      formula <- stats::as.formula(formula_str)
      model <- tryCatch(survival::survreg(formula, data = dataYX, dist = dist, x = FALSE, y = FALSE),
                        error = function(e) {warning("Failed to fit model with var", var, ": ", e$message); NULL})

      if (!is.null(model)) {
        aic_values[i] <- stats::AIC(model)
      } else {
        aic_values[i] <- Inf # Penalize models that fail to fit
      }
    }

    best_candidate_idx <- which.min(aic_values)
    best_candidate_aic <- aic_values[best_candidate_idx]
    best_candidate_var <- candidate_vars[best_candidate_idx]

    # Add variable if it improves AIC
    if (best_candidate_aic < best_aic) {
      print(paste("Adding", best_candidate_var, "AIC:", round(best_candidate_aic, 2), "(Improvement:", round(best_aic - best_candidate_aic, 2), ")"))
      selected_vars <- c(selected_vars, best_candidate_var)
      candidate_vars <- setdiff(candidate_vars, best_candidate_var)
      best_aic <- best_candidate_aic
      improved <- TRUE
    } else {
      # If no improvement, but no variable has been selected yet, pick the best single variable
      if (length(selected_vars) == 0) {
        print(paste("No improvement over intercept-only, but keeping best single variable:", best_candidate_var, "AIC:", round(best_candidate_aic, 2)))
        selected_vars <- best_candidate_var
      } else {
        print(paste("No improvement adding remaining variables. Best AIC:", round(best_aic, 2)))
      }
      break # Stop if no variable improves AIC
    }
  } # End while loop

  # Fit the final selected model
  if (length(selected_vars) > 0) {
      final_formula_str <- paste("survival::Surv(", timevar, ",", eventvar, ") ~", paste(selected_vars, collapse = "+"))
      final_model <- survival::survreg(stats::as.formula(final_formula_str), data = dataYX, dist = dist, x = TRUE, y = TRUE) # Keep x, y for prediction
  } else {
      # If no variables selected, use the intercept-only model
      warning("No variables selected by forward selection. Using intercept-only model.")
      final_model <- null_model
      # Need to refit with x=T, y=T if null_model didn't have them
      final_model <- survival::survreg(null_formula, data = dataYX, dist = dist, x = TRUE, y = TRUE)
      selected_vars <- character(0) # Ensure selected_vars is empty character vector
  }


  if (is.null(final_model)) {
    stop("No valid model could be fit even after selection. Please check the data and variables.")
  }

  # Get unique event times from training data
  times <- sort(unique(dataYX[dataYX[[eventvar]] == 1, timevar]))

  result <- list(survregOut = final_model,
              times = times,
              varprof = varprof)
  class(result) <- c("ml4t2e_surv_survreg", "list")
  return(result)
}


#' @title Predict_SurvModel_SurvReg
#'
#' @description Make predictions using a fitted parametric survival model.
#' @param modelout the output from 'SurvModel_SurvReg'
#' @param newdata the data for which the predictions are to be calculated
#' @param times optional vector of times for prediction. If NULL, uses times from modelout.
#'
#' @return a list containing the following objects:
#' Probs: predicted survival probability matrix (rows=times, cols=observations),
#' Times: the times at which the probabilities are calculated (including 0).
#'
#' @importFrom stats predict
#' @export
Predict_SurvModel_SurvReg <- function(modelout, newdata, new_times = NULL) {
  cat("Predict_SurvModel_SurvReg called with", nrow(newdata), "rows\n")
  # Select only the variables used in the final model
  selected_vars <- attr(terms(modelout$survregOut), "term.labels")

  # Check if selected_vars exist in newdata
  missing_pred_vars <- setdiff(selected_vars, colnames(newdata))
  if (length(missing_pred_vars) > 0) {
      stop("Predictor variables missing from newdata: ", paste(missing_pred_vars, collapse=", "))
  }
  # Subset newdata to include only necessary predictors (if any)
  if (length(selected_vars) > 0) {
      X <- newdata[, selected_vars, drop = FALSE]
  } else {
      # Handle intercept-only model case
      X <- newdata[1, character(0), drop = FALSE] # Create empty df for prediction if no vars
      # Need to replicate for all rows in newdata if prediction depends on it
      X <- X[rep(1, nrow(newdata)), , drop = FALSE]
  }

  # Get linear predictors using model matrix approach
  # Create model matrix for prediction using only the predictor terms
  if (length(selected_vars) > 0) {
    pred_formula <- as.formula(paste("~", paste(selected_vars, collapse = "+")))
    pred_terms <- terms(pred_formula)
    X_model <- model.matrix(pred_terms, data = X)
  } else {
    # Intercept-only model
    X_model <- matrix(1, nrow = nrow(X), ncol = 1)
    colnames(X_model) <- "(Intercept)"
  }
  
  coefficients <- coef(modelout$survregOut)
  
  # Ensure coefficient names match model matrix columns
  coef_names <- names(coefficients)
  mm_names <- colnames(X_model)
  
  # Match coefficients to model matrix columns
  matched_coefs <- rep(0, ncol(X_model))
  names(matched_coefs) <- mm_names
  
  for (name in coef_names) {
    if (name %in% mm_names) {
      coef_val <- coefficients[name]
      # Handle NA coefficients (unseen factor levels) by setting to 0
      if (!is.na(coef_val)) {
        matched_coefs[name] <- coef_val
      } else {
        warning("Coefficient for ", name, " is NA (unseen factor level), setting to 0")
        matched_coefs[name] <- 0
      }
    }
  }
  
  # Compute linear predictors
  linear_preds <- X_model %*% matched_coefs
  
    # Get model parameters
  dist <- modelout$survregOut$dist
  
  # Determine time points for prediction
  if (is.null(new_times)) {
    times <- modelout$times # Use times from training fit
  } else {
    times <- new_times
  }
  times <- sort(unique(c(0, times))) # Ensure 0 included
  
  # Compute survival probabilities based on distribution
  n_obs <- nrow(X)
  n_times <- length(times)
  surv_probs <- matrix(0.0, nrow = n_times, ncol = n_obs)  # Ensure numeric
  
  for (i in 1:n_obs) {
    lp <- linear_preds[i]
    if (dist == "exponential") {
      rate <- exp(-lp)
      surv_probs[, i] <- exp(-rate * times)
    } else {
      surv_probs[, i] <- NA
    }
  }

  # Result is matrix: rows=times, cols=observations
  return(list(Probs = surv_probs, Times = times))
}
