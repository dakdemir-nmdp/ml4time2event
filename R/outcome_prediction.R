#' @title PredictAllPossibleOutcomesSurvOrCifs
#'
#' @description Get predictions (Survival Probs or CIFs) for a given dataset from a list of fitted models.
#' Handles both survival ('SURV') and competing risks ('CR') models based on `modeltypes`.
#'
#' @param data Data frame for which to generate predictions.
#' @param modelslist A list where each element is a fitted model object (e.g., output from `RunSurvModels` or `RunCRModels`).
#' @param modeltypes A character vector of the same length as `modelslist`, indicating the type of each model ("SURV" or "CR").
#' @param times Numeric vector of time points at which to generate predictions.
#' @return A list of the same length as `modelslist`, where each element contains the prediction output
#'   (typically including interpolated probabilities/CIFs) from the corresponding `PredictSurvModels` or `PredictCRModels` function.
#'   Returns NA for unsupported model types.
#' @export
PredictAllPossibleOutcomesSurvOrCifs<-function(data, modelslist, modeltypes, times){

  if (length(modelslist) != length(modeltypes)) {
      stop("Length of 'modelslist' and 'modeltypes' must be equal.")
  }

  predictions_list <- lapply(1:length(modelslist), function(i){
    model_obj <- modelslist[[i]]
    model_type <- toupper(modeltypes[i]) # Ensure uppercase
    cat("Predicting using model", i, "(Type:", model_type, ")...\n")

    if (is.null(model_obj)) {
        warning("Model object at index", i, "is NULL. Skipping prediction.")
        return(NA)
    }

    if (model_type == "SURV"){
      # Assuming PredictSurvModels is loaded/available
      pred_out <- tryCatch(
          PredictSurvModels(models=model_obj, newdata=data, newtimes=times),
          error = function(e) {
              warning("Error predicting SURV model at index ", i, ": ", e$message)
              return(NA)
          }
      )
      return(pred_out)
    } else if (model_type == "CR"){
      # Assuming PredictCRModels is loaded/available
       pred_out <- tryCatch(
          PredictCRModels(models=model_obj, newdata=data, newtimes=times),
           error = function(e) {
              warning("Error predicting CR model at index ", i, ": ", e$message)
              return(NA)
          }
       )
       return(pred_out)
    } else {
      warning("Unsupported model type '", modeltypes[i], "' at index ", i, ". Skipping prediction.")
      return(NA)
    }
  })

  return(predictions_list)
}



#' @title CalculateExpectedTimeLost
#'
#' @description Calculate Restricted Mean Time Lost (RMTL) from predicted survival or CIF curves.
#' RMTL is the integral of the event probability (1-Survival or CIF) up to a time limit UL.
#'
#' @param PredictedCurves A list where each element is the prediction output from
#'   `PredictAllPossibleOutcomesSurvOrCifs` (or similar), containing an ensemble prediction matrix
#'   named `NewProbs` (rows=times, cols=observations).
#' @param modeltypes Character vector indicating the type ("SURV" or "CR") for each element in `PredictedCurves`.
#' @param times Numeric vector of time points corresponding to the rows of the prediction matrices.
#' @param UL Upper limit of integration for RMTL.
#' @param LL Lower limit of integration (default: 0).
#' @return A list of numeric vectors, where each vector contains the calculated RMTL for each observation
#'   corresponding to the respective input prediction curve.
#' @export
CalculateExpectedTimeLost<-function(PredictedCurves, modeltypes, times, UL, LL=0){
  if (length(PredictedCurves) != length(modeltypes)) {
      stop("Length of 'PredictedCurves' and 'modeltypes' must be equal.")
  }
  if (!is.numeric(UL) || length(UL) != 1 || UL <= LL) {
      stop("'UL' must be a single numeric value greater than 'LL'.")
  }
   if (!is.numeric(LL) || length(LL) != 1 || LL < 0) {
      stop("'LL' must be a single non-negative numeric value.")
  }

  ExpectedTimeLostList <- lapply(1:length(PredictedCurves), function(i) {
      pred_output <- PredictedCurves[[i]]
      model_type <- toupper(modeltypes[i])

      if (is.null(pred_output) || (is.atomic(pred_output) && length(pred_output) == 1 && is.na(pred_output)) || is.null(pred_output$NewProbs)) {
          warning("Invalid prediction output or missing 'NewProbs' at index ", i, ". Cannot calculate time lost.")
          return(NA) # Or perhaps NULL or vector of NAs?
      }

      ProbMatrix <- pred_output$NewProbs # Ensemble predictions (rows=times, cols=obs)

      # Calculate event probability: 1-Survival for SURV, CIF for CR
      if (model_type == "SURV"){
        EventProbMatrix <- 1 - ProbMatrix
      } else if (model_type == "CR"){
        EventProbMatrix <- ProbMatrix # Assuming NewProbs contains CIF for event of interest
      } else {
        warning("Unsupported model type '", modeltypes[i], "' at index ", i, ". Cannot calculate time lost.")
        return(NA)
      }

      # Integrate event probability for each observation using the Integrator helper
      # Assuming Integrator is loaded/available
      time_lost_vector <- tryCatch(
          apply(EventProbMatrix, 2, function(obs_event_probs){
              Integrator(times=times, scores=obs_event_probs, minmax=c(LL, UL), scale = FALSE) # Integrate, don't scale
          }),
          error = function(e) {
              warning("Error during integration for index ", i, ": ", e$message)
              return(rep(NA, ncol(EventProbMatrix))) # Return NA vector on error
          }
      )
      return(time_lost_vector)
  })

  return(ExpectedTimeLostList)
}
