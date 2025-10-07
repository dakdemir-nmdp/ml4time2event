#' ZeroOneScalerData
#'
#' Scale numeric variables in a data frame to the range [0, 1].
#' Records the min and max values used for scaling.
#'
#' @param data Data frame with numeric columns to scale.
#' @return A list containing:
#'   - data: The data frame with numeric columns scaled.
#'   - minxvec: A named numeric vector of the minimum values used for scaling each column.
#'   - maxxvec: A named numeric vector of the maximum values used for scaling each column.
#' @export
ZeroOneScalerData<-function(data){
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  datanew <- data
  minxvec <- rep(NA_real_, ncol(data))
  maxxvec <- rep(NA_real_, ncol(data))
  names(minxvec) <- names(maxxvec) <- colnames(data)

  for (i in 1:ncol(datanew)){
    col_data <- datanew[[i]]
    if (is.numeric(col_data)){
      minx <- min(col_data, na.rm = TRUE)
      maxx <- max(col_data, na.rm = TRUE)
      range_val <- maxx - minx

      # Avoid division by zero if range is 0 (column is constant)
      if (range_val > .Machine$double.eps) { # Use machine epsilon for comparison
          datanew[[i]] <- (col_data - minx) / range_val
      } else {
          # If range is zero, set scaled value to 0, 0.5, or 1?
          # Setting to 0.5 might be reasonable, or just leave as is (0 if minx=maxx=0).
          # Let's set to 0 if range is effectively zero.
          datanew[[i]] <- rep(0, length(col_data))
          warning("Column '", colnames(datanew)[i], "' has zero range. Scaled values set to 0.")
      }
      minxvec[i] <- minx
      maxxvec[i] <- maxx
    }
  }
  return(list(data=datanew, minxvec=minxvec, maxxvec=maxxvec ))
}



#' ZeroOneScalerApplierData
#'
#' Apply a pre-calculated 0-1 scaling to numeric variables in a new data frame.
#' Uses minimum and maximum values provided (e.g., from a training set).
#'
#' @param data Data frame with numeric columns to scale.
#' @param mins A named numeric vector of minimum values (names must match columns in 'data').
#' @param maxs A named numeric vector of maximum values (names must match columns in 'data').
#' @return A data frame with numeric columns scaled according to the provided mins and maxs.
#' @export
ZeroOneScalerApplierData<-function(data, mins, maxs){
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (is.null(names(mins)) || is.null(names(maxs))) stop("'mins' and 'maxs' must be named vectors.")

  datanew <- data
  numeric_cols <- names(which(sapply(datanew, is.numeric)))

  for (i in numeric_cols){
    if (i %in% names(mins) && i %in% names(maxs)) {
        minx <- mins[i]
        maxx <- maxs[i]
        range_val <- maxx - minx

        if (!is.na(minx) && !is.na(maxx)) { # Ensure min/max are valid
            col_data <- datanew[[i]]
            if (range_val > .Machine$double.eps) {
                datanew[[i]] <- (col_data - minx) / range_val
            } else {
                # Handle zero range case consistent with ZeroOneScalerData
                datanew[[i]] <- rep(0, length(col_data))
            }
            # Optional: Add capping for values outside the original training range?
            # datanew[[i]][datanew[[i]] < 0] <- 0
            # datanew[[i]][datanew[[i]] > 1] <- 1
        } else {
            warning("Missing min/max value for numeric column '", i, "'. Skipping scaling.")
        }
    } else {
        warning("Min/max value not provided for numeric column '", i, "'. Skipping scaling.")
    }
  }
  return(datanew)
}

#' UndoZeroOneScalerApplierData
#'
#' Reverse a 0-1 scaling transformation on numeric variables using the original min/max values.
#'
#' @param data Data frame with scaled numeric columns (values expected between 0 and 1).
#' @param mins A named numeric vector of the original minimum values used for scaling.
#' @param maxs A named numeric vector of the original maximum values used for scaling.
#' @return A data frame with numeric columns un-scaled back to their original range.
#' @export
UndoZeroOneScalerApplierData<-function(data, mins, maxs){
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (is.null(names(mins)) || is.null(names(maxs))) stop("'mins' and 'maxs' must be named vectors.")

  datanew <- data
  numeric_cols <- names(which(sapply(datanew, is.numeric)))

  for (i in numeric_cols){
     if (i %in% names(mins) && i %in% names(maxs)) {
        minx <- mins[i]
        maxx <- maxs[i]
        range_val <- maxx - minx

         if (!is.na(minx) && !is.na(maxx)) { # Ensure min/max are valid
            col_data <- datanew[[i]]
            if (range_val > .Machine$double.eps) {
                datanew[[i]] <- (col_data * range_val) + minx
            } else {
                # If range was zero, original value was minx (or maxx)
                datanew[[i]] <- rep(minx, length(col_data))
            }
         } else {
             warning("Missing min/max value for numeric column '", i, "'. Skipping un-scaling.")
         }
     } else {
         warning("Min/max value not provided for numeric column '", i, "'. Skipping un-scaling.")
     }
  }
  return(datanew)
}





#' NumVarstCatsData
#'
#' Convert numeric variables to categorical (factor) variables based on quantiles or specified cuts.
#'
#' @param data Data frame containing numeric variables to categorize.
#' @param numgroups Integer. If specified, numeric variables (with enough unique values) are cut into this many quantile groups.
#' @param cuts Numeric vector. If specified, numeric variables are cut using these specific values as breaks. `numgroups` is ignored if `cuts` is provided.
#' @param min_unique_vals Integer. Minimum number of unique non-NA values a numeric variable must have to be categorized (default: 5).
#' @return Data frame with specified numeric variables converted to factors.
#' @importFrom stats quantile
#' @export
NumVarstCatsData<-function(data, numgroups=NULL, cuts=NULL, min_unique_vals = 5){
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (is.null(numgroups) && is.null(cuts)) stop("Must provide either 'numgroups' or 'cuts'.")
  if (!is.null(numgroups) && !is.null(cuts)) warning("'cuts' provided, 'numgroups' will be ignored.")

  datanew <- data
  numeric_cols <- names(which(sapply(datanew, is.numeric)))

  # --- Internal Helper: quantcat (modified from Hmisc::cut2) ---
  # Creates quantile groups, handling duplicated quantiles by making them distinct points.
  quantcat <- function (x, q = 4) {
      q_probs <- seq(0, 1, length.out = q + 1)
      quant <- stats::quantile(x, q_probs, na.rm = TRUE, type = 7) # Use default quantile type
      dups <- duplicated(quant)

      if (any(dups)) {
          # Handle duplicated quantiles - make them distinct points
          unique_quants <- unique(quant)
          new_breaks <- sort(unique_quants)
          # Ensure min and max are included if they were duplicated away
          if (!min(x, na.rm=TRUE) %in% new_breaks) new_breaks <- c(min(x, na.rm=TRUE), new_breaks)
          if (!max(x, na.rm=TRUE) %in% new_breaks) new_breaks <- c(new_breaks, max(x, na.rm=TRUE))
          new_breaks <- unique(sort(new_breaks))

          if (length(new_breaks) < 2) { # Cannot cut if only one unique value
              warning("Cannot create quantile groups; too few unique values after handling duplicates.")
              return(factor(rep(paste(range(x, na.rm=TRUE), collapse="-"), length(x)))) # Return single level factor
          }
          retval <- cut(x, breaks = new_breaks, include.lowest = TRUE, right = FALSE, dig.lab = 4) # Use right=FALSE for [a,b) intervals
      } else {
          # No duplicates, standard cut
          if (length(quant) < 2) {
               warning("Cannot create quantile groups; too few unique quantile values.")
               return(factor(rep(paste(range(x, na.rm=TRUE), collapse="-"), length(x))))
          }
          retval <- cut(x, breaks = quant, include.lowest = TRUE, right = FALSE, dig.lab = 4)
      }
      # Clean up level names (optional)
      levels(retval) <- gsub(",", ", ", levels(retval)) # Add space after comma
      return(retval)
  }

  # --- Internal Helper: cutcat ---
  # Cuts based on provided break points.
  cutcat <- function (x, cuts) {
      breaks <- sort(unique(c(-Inf, cuts, Inf))) # Add Inf/-Inf for full range
      if (length(breaks) < 2) {
          warning("Cannot cut with provided breaks; too few unique break points.")
          return(factor(rep(paste(range(x, na.rm=TRUE), collapse="-"), length(x))))
      }
      retval <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE, dig.lab = 4)
      levels(retval) <- gsub(",", ", ", levels(retval))
      return(retval)
  }

  # --- Apply Categorization ---
  for (i in numeric_cols){
    col_data <- datanew[[i]]
    if (length(unique(col_data[!is.na(col_data)])) >= min_unique_vals) {
        cat("Categorizing variable '", i, "'...\n")
        if (!is.null(cuts)) {
            datanew[[i]] <- cutcat(col_data, cuts = cuts)
        } else { # Use numgroups
            datanew[[i]] <- quantcat(col_data, q = numgroups)
        }
    } else {
        cat("Skipping categorization for variable '", i, "' (less than ", min_unique_vals, " unique values).\n")
    }
  }
  return(datanew)
}
