## Module: Dec Tools
## This file contains functions for auxiliary decision tools and scenario generation.
## Sections include:
##   - Conditional Categories Extraction (GetAllConditionalCategoriesData)
##   - HLA Scenarios Permuter (HLAVarsPermuter)
##   - Outcome Prediction Functions (PredictAllPossibleOutcomesSurvOrCifs, CalculateExpectedTimeLost)
##   - Permutation Importance and Design Generation (permutationImportance, makePermutedDesign)


#' @title GetAllCondidtionalCategoriesData
#'
#' @description  Get possible conditional categories  from data
#' all conditions are equality conditions expressed as data.frame
#'
#' @export
GetAllConditionalCategoriesData<-function(data, conditions=NULL, minmatch=8, mustmatchcols=NULL, freevars=NULL){

  if (!is.null(freevars)){
  data0<-data[,c(freevars)]
  }

  if (!is.null(conditions)){
    condMat<-matrix(0,ncol=ncol(conditions),nrow=nrow(data))
    colnames(condMat)<-colnames(conditions)
    for (vari in colnames(conditions)){
      condMat[,vari]<-as.numeric(data[,vari]==conditions[,vari])
    }
  }
  data<-data[which(rowSums(condMat)>=minmatch),]
  if (!is.null(mustmatchcols)){
    for (vari in mustmatchcols){
      data<-data[which(data[, vari]==conditions[,vari]),]
    }
  }

  if (!is.null(freevars)){
    datanew<-data0[,freevars]
    datanew<-datanew[!duplicated(datanew),]
    datanew<-expand.grid(datanew)
    datanew<-datanew[!duplicated(datanew),]
    nrowd<-nrow(datanew)
    dataall<-NULL
    for (repi in 1:nrowd){
      dataall<-rbind(dataall,data)
      print(dim(dataall))
    }
    dataall<-cbind(dataall[, !colnames(data)%in%freevars], datanew)
  }


  dataall[!duplicated(dataall),]
}



#' @title HLAVarsPermuter
#'
#' @description Write scenarios for HLA matching
#'
#'
#' @export
HLAVarsPermuter<-function(dataforpatient=NULL, HLAVarsdata=NULL){
  dataforpatient0<-dataforpatient[,!colnames(dataforpatient)%in%colnames(HLAVarsdata)]
  dataforpatientNew<-NULL
  if (!is.null(HLAVarsdata)){
    for (i in 1:nrow(HLAVarsdata)){
    dataforpatientNew<-rbind(dataforpatientNew,cbind(dataforpatient0, HLAVarsdata[i,]))
    }
  }
  dataforpatientNew
}




#' @title PredictAllPossibleOutcomesSurvOrCifs
#'
#' @description Given a list of models and data get predictions
#'
#'
#' @export
PredictAllPossibleOutcomesSurvOrCifs<-function(data, modelslist, modeltypes, times){

 lapply(1:length(data), function(vari){
    if (modeltypes[vari]=="SURV"){
      return(PredictSurvModels(models=modelslist[[vari]], newdata=data, newtimes=times))
    } else if (modeltypes[vari]=="CR"){
      return(PredictCRModels(models=modelslist[[vari]], newdata=data, newtimes=times))
    } else {
      NA
    }
  })

}



#' @title CalculateExpectedTimeLost
#'
#' @description Predicted time lost frm survival or cr curves
#'
#' @export
CalculateExpectedTimeLost<-function(PredictedCurves, modeltypes, times, UL,LL){
  for (vari in 1:length(PredictedCurves)){
    if (modeltypes[[vari]]=="SURV"){
      PredictedCurves[[vari]]<-1-PredictedCurves[[vari]]
    }
  }

  PredictedTimeLost<-lapply(PredictedCurves,
                            function(X){
    apply(X, 2,function(x){Integrator(times=times, scores=x, minmax=c(LL,UL))})
    }
    )
  return(PredictedTimeLost)
}









###################

#' @title computes permutation importance
#' @description computes the change in prediction error from permuting variables.
#'
#'
#' @param data a \code{data.frame} including both \code{y} and \code{vars}.
#' @param vars a character vector specifying columns of \code{data} to permute.
#' @param y a character vector giving the name of the target/outcome variable.
#' @param model an object with a predict method which returns a vector or matrix. presumably this object represents a model fit.
#' @param nperm positive integer giving the number of times to permute the indicated variables (default is 100).
#' @param predict.fun what function to generate predictions using \code{model}. default is the predict method for \code{model}. the function must take two arguments, \code{object} and \code{newdata} and should return a vector or matrix.
#' @param loss.fun what loss function to use to measure prediction errors. default is mean squared-error for ordered predictions and mean misclassification error for unordered prediction errors. this function must take two arguments, \dQuote{x} and \dQuote{y}, which operate on the output of \code{predict.fun} and \code{data[, y]}.
#' @param contrast.fun what function to use to contrast the permuted and unpermuted predictions. default is the difference. this function takes two arguments \dQuote{x} and \dQuote{y}, which are the output of the \code{loss.fun}.
#'
#' @return a numeric vector or matrix, depending on \code{contrast.fun} and \code{loss.fun}, giving the change in prediction error from \code{nperm} permutations of \code{vars}.
#'
#' @examples
#' X = replicate(3, rnorm(100))
#' y = X %*% runif(3)
#' data = data.frame(X, y)
#' fit = lm(y ~ -1 + X1 + X2 + X3, data)
#'
#' permutationImportance(data, "X1", "y", fit)
#' @export
permutationImportance = function(data, vars, y, model,
                                 nperm = 100L,
                                 predict.fun = function(object, newdata) predict(object, newdata = newdata),
                                 loss.fun = function(x, y) sum((x- y)^2),
                                 contrast.fun = function(x, y) x - y) {


  design = makePermutedDesign(data, vars, nperm)
  unpermuted = loss.fun(predict.fun(model, data), data[, y])
  permuted.predictions = predict.fun(model, design)
  permuted = loss.fun(permuted.predictions, design[, y])
  if (length(permuted) == nperm * nrow(data)) {
    permuted <- array(permuted, c(nrow(data), nperm))
  }
  contrast.fun(permuted, unpermuted)
}

#' @title creates a \code{data.frame} with some columns permuted
#' @description takes an input data.frame, permutes some variables, and stacks the resulting \code{data.frame}s.
#'
#'
#' @param data a \code{data.frame} a subset of which must be \code{vars}.
#' @param vars a character vector indicating columns in \code{data} to permute.
#' @param nperm an integer specifying the number of times to permute the columns indicated by \code{vars}.
#'
#' @return a \code{data.frame} with number of rows equal to \code{nrow(data) * nperm}
#'
#' @examples
#' data = data.frame(x = 1:3, y = letters[1:3])
#' makePermutedDesign(data, "x", 3)
#' @export
makePermutedDesign = function(data, vars, nperm) {
  design = data[rep(1:nrow(data), times = nperm), ]
  idx = lapply(1:nperm, function(x) sample(1:nrow(data)) + (x - 1) * nrow(data))
  idx = unlist(idx)
  design[, vars] = design[idx, vars]
  design
}
