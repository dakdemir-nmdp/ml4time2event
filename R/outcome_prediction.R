#'
#'
PredictAllPossibleOutcomesSurvOrCifs<-function(data, modelslist, modeltypes, times){

  if (length(modelslist) != length(modeltypes)) {
      stop("Length of 'modelslist' and 'modeltypes' must be equal.")
  }

  predictions_list <- lapply(1:length(modelslist), function(i){
    model_obj <- modelslist[[i]]
    model_type <- toupper(modeltypes[i]) # Ensure uppercase
    #cat("Predicting using model", i, "(Type:", model_type, ")...\n")

    if (is.null(model_obj)) {
        warning("Model object at index", i, "is NULL. Skipping prediction.")
        return(NA)
    }

    if (model_type == "SURV"){
      # Assuming PredictSurvModels is loaded/available
      pred_out <- tryCatch(
          PredictSurvModels(models=model_obj, newdata=data, new_times=times),
          error = function(e) {
              warning("Error predicting SURV model at index ", i, ": ", e$message)
              return(NA)
          }
      )
      return(pred_out)
    } else if (model_type == "CR"){
      # Assuming PredictCRModels is loaded/available
       pred_out <- tryCatch(
          PredictCRModels(models=model_obj, newdata=data, new_times=times),
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



