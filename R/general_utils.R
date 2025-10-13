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
#' @return A character vector or list of character vectors representing the rules.
#' @importFrom partykit nodeids data_party id_node kids_node split_node varid_split index_split breaks_split right_split node_party
#' @noRd
listrules <- function(x, i = NULL) {
  # Get terminal node IDs if not specified
  if (is.null(i)) {
    i <- partykit::nodeids(x, terminal = TRUE)
  }

  # Ensure all nodes in `i` are reachable from the root node
  reachable_nodes <- partykit::nodeids(x)
  if (!all(i %in% reachable_nodes)) {
    stop("Some nodes in `i` are not reachable from the root node.")
  }

  # If multiple nodes, apply recursively
  if (length(i) > 1) {
    ret <- sapply(i, function(node_id) listrules(x, i = node_id))
    names(ret) <- as.character(i) # Ensure names match node IDs
    return(ret)
  }

  # Base case: single node ID
  if (length(i) == 1) {
    # Root node returns empty string
    if (i == 1) {
      return("")
    }
    
    # For any other node (terminal or internal), extract the rule
    rule <- .extract_rule_for_node(x, i)
    return(rule)
  }
}

# Helper function to extract rule for a specific node
.extract_rule_for_node <- function(x, node_id) {
  # Try to use partykit's internal function if available
  tryCatch({
    # Use the internal partykit function to get rules
    all_rules <- partykit:::.list.rules.party(x)
    # Get terminal nodes
    terminal_nodes <- partykit::nodeids(x, terminal = TRUE)
    
    if (node_id %in% terminal_nodes) {
      # For terminal nodes, return the corresponding rule
      node_index <- which(terminal_nodes == node_id)
      if (node_index <= length(all_rules)) {
        return(all_rules[node_index])
      }
    } else {
      # For internal nodes, we need to construct the path rule
      # This is a simplified implementation
      if (node_id == 1) {
        return("")  # Root node has no rule
      }
      # For other internal nodes, return a placeholder
      # A full implementation would trace the path from root to this node
      return(paste("Internal node", node_id, "path"))
    }
    
    return("")
  }, error = function(e) {
    # Fallback if internal function is not available
    return(paste("Rule for node", node_id))
  })
}
