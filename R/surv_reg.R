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
#' @importFrom stats AIC as.formula predict quantile
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
       print(paste("No improvement adding remaining variables. Best AIC:", round(best_aic, 2)))
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

  return(list(survregOut = final_model,
              times = times,
              varprof = varprof))
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
Predict_SurvModel_SurvReg <- function(modelout, newdata, newtimes = NULL) {
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


  # Predict quantiles (survival times) for the new data
  p_seq <- seq(0.01, 0.99, length.out = 100)
  predicted_times_for_quantiles <- stats::predict(modelout$survregOut, newdata = X, type = "quantile", p = p_seq)

  # Determine time points for interpolation
  if (is.null(newtimes)) {
    times.interest <- modelout$times # Use times from training fit
  } else {
    times.interest <- sort(unique(c(0, newtimes))) # Use provided times, ensure 0 included
  }

  # Interpolate survival probabilities at desired times.interest
  # Assuming survivalProbsInterpolator is loaded/available
  predicttest <- apply(predicted_times_for_quantiles, 1, function(pred_times_row) {
      valid_idx <- !is.na(pred_times_row) & !is.infinite(pred_times_row)
      if (sum(valid_idx) < 2) return(rep(NA, length(times.interest)))

      sorted_idx <- order(pred_times_row[valid_idx])
      t_pred <- pred_times_row[valid_idx][sorted_idx]
      s_pred <- (1 - p_seq)[valid_idx][sorted_idx]

      t_interp <- c(0, t_pred)
      s_interp <- c(1, s_pred)

      survivalProbsInterpolator(times.interest, probs = s_interp, times = t_interp)
  })
  # Result is matrix: rows=times.interest, cols=observations

  return(list(Probs = predicttest, Times = times.interest))
}
