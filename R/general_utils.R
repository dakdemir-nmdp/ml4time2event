#' @title VariableProfile
#' @description Creates a profile of variables in a dataset, summarizing factors and numeric ranges.
#' @param data data frame
#' @param expvars character vector of variable names to profile
#' @return a list named by expvars, containing tables for factors/characters or min/max for numerics.
#' @export
VariableProfile<-function(data, expvars){
  varprofile<-vector(mode="list", length=length(expvars))
  names(varprofile)<-expvars
  for (vari in expvars){
    if (vari %in% colnames(data)) { # Check if variable exists in data
        col_data <- data[[vari]] # Use [[ ]] for safer column access
        if (is.factor(col_data)){
            varprofile[[vari]]<-table(col_data, useNA = "ifany") # Include NAs in table
        } else if (is.numeric(col_data)){
            varprofile[[vari]]<-c(min=min(col_data, na.rm = TRUE), max=max(col_data, na.rm = TRUE)) # Add names
        } else if (is.character(col_data)){
            varprofile[[vari]]<-table(col_data, useNA = "ifany") # Include NAs
        } else {
            varprofile[[vari]] <- paste("Unsupported type:", class(col_data))
        }
    } else {
        varprofile[[vari]] <- "Variable not found in data"
    }
  }
  varprofile
}

#' @title listrules
#' @description Extract rules from a partykit tree object.
#' (Adapted from partykit:::.list.rules.party)
#' @param x A party object representing a tree.
#' @param i Node ID(s) to extract rules for (default: terminal nodes).
#' @param ... Additional arguments (not used).
#' @return A character vector or list of character vectors representing the rules.
#' @importFrom partykit nodeids data_party id_node kids_node split_node varid_split index_split breaks_split right_split node_party
#' @noRd
listrules<-function (x, i = NULL, ...)
{
  # Get terminal node IDs if not specified
  if (is.null(i))
    i <- partykit::nodeids(x, terminal = TRUE)

  # If multiple nodes, apply recursively
  if (length(i) > 1) {
    ret <- sapply(i, listrules, x = x)
    # Use node IDs as names if available and numeric, otherwise default names
    names(ret) <- if (is.numeric(i) && !is.null(names(x))) names(x)[i] else i
    return(ret)
  }

  # Ensure single numeric node ID
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)

  # Get data associated with the node
  dat <- partykit::data_party(x, i)

  # Handle fitted values if present (though typically not needed for rule extraction)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    # fit <- dat[, findx:ncol(dat), drop = FALSE] # Not used for rules
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data # Fallback to original data if needed
  } else {
    # fit <- NULL
    dat <- x$data # Use original data if no fitted values stored
  }

  # Recursive function to build the rule string
  rule <- c()
  recFun <- function(node) {
    # Base case: reached the target node i
    if (partykit::id_node(node) == i)
      return(NULL)

    # Find which child node leads towards the target node i
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    # Ensure kid is numeric before comparison
    if (!is.numeric(kid)) {
        warning(paste("Non-numeric kid node ID found:", kid, "at node", partykit::id_node(node)))
        return(NULL) # Stop recursion for this branch
    }
    # Check if i is reachable from this node's children
    if (!any(kid <= i)) {
         warning(paste("Target node", i, "not reachable from node", partykit::id_node(node)))
         return(NULL) # Stop recursion
    }
    whichkid <- max(which(kid <= i)) # Find the child subtree containing i

    # Get split information
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar] # Variable name used for splitting

    # Handle missing variable names gracefully
    if (is.na(svar) || is.null(svar)) {
        warning(paste("Missing variable name for varid", ivar, "at node", partykit::id_node(node)))
        svar <- paste0("var", ivar) # Use a placeholder name
    }


    index <- partykit::index_split(split) # Index for categorical splits

    # Construct rule based on variable type
    if (is.factor(dat[[svar]])) { # Use [[ ]] for safer access
      # Factor split
      if (is.null(index)) # Binary split? Infer index based on breaks
        index <- ((1:nlevels(dat[[svar]])) > partykit::breaks_split(split)) + 1
      slevels <- levels(dat[[svar]])[index == whichkid]
      # Escape special characters in levels for the rule string
      slevels_escaped <- gsub("\"", "\\\"", slevels)
      srule <- paste0(svar, " %in% c(\"", paste(slevels_escaped, collapse = "\", \""), "\")")
    } else {
      # Numeric or other type split
      if (is.null(index)) # If index is null, assume binary split
        index <- 1:length(kid) # Should correspond to kids
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), Inf))
      # Ensure breaks matrix has correct dimensions if only one break
       if (length(partykit::breaks_split(split)) == 1 && nrow(breaks) == 1) {
           breaks <- rbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), Inf))
       }
      # Check if index is valid for breaks
      if (whichkid > nrow(breaks)) {
          warning(paste("Invalid 'whichkid' index", whichkid, "for breaks at node", partykit::id_node(node)))
          srule <- paste0(svar, " ERROR_invalid_split_index")
      } else {
          sbreak <- breaks[whichkid, ]
          right <- partykit::right_split(split) # Which side is TRUE (right or left)?
          srule <- c()
          if (is.finite(sbreak[1]))
            srule <- c(srule, paste(svar, ifelse(right, ">", ">="), format(sbreak[1], digits = 3))) # Format numbers
          if (is.finite(sbreak[2]))
            srule <- c(srule, paste(svar, ifelse(right, "<=", "<"), format(sbreak[2], digits = 3))) # Format numbers
          srule <- paste(srule, collapse = " & ")
      }
    }

    # Prepend the current rule segment
    rule <<- c(rule, srule)

    # Recurse down the correct child node
    # Check if the child node exists before recursing
    child_node <- node[[whichkid]]
    if (!is.null(child_node)) {
        return(recFun(child_node))
    } else {
        warning(paste("Child node", whichkid, "is NULL at node", partykit::id_node(node)))
        return(NULL) # Stop recursion
    }
  }

  # Start recursion from the root node
  recFun(partykit::node_party(x))

  # Combine rule segments
  paste(rev(rule), collapse = " & ") # Reverse rules to get path from root
}
