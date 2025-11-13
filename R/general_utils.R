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

format_model_name <- function(model_names, model_type = "survival") {

  # Define mapping for survival models
  surv_mapping <- c(
    RF_Model = "Random Forest",
    RF_Model2 = "Random Forest (Top Vars)",
    glmnet_Model = "GLMNet",
    CPH_Model = "Cox PH",
    bart_Model = "BART",
    shallownn_Model = "Shallow NN",
    gam_Model = "GAM",
    gbm_Model = "GBM",
    survregexp_Model = "ExpSurvReg",
    survregweib_Model = "WeibSurvReg",
    xgboost_Model = "XGBoost",
    RuleFit_Model = "RuleFit",
    ttah_Model = "TTAH",
    Ensemble = "Ensemble"
  )

  # Define mapping for competing risks models
  cr_mapping <- c(
    Cox_Model = "Cox (Cause-Specific)",
    FG_Model = "Fine-Gray",
    RF_Model = "Random Forest",
    glmnet_Model = "GLMNet",
    bart_Model = "BART",
    shallownn_Model = "Shallow NN",
    gam_Model = "GAM",
    xgboost_Model = "XGBoost",
    RuleFit_Model = "RuleFit",
    ttah_Model = "TTAH",
    survregexp_Model = "ExpSurvReg",
    Ensemble = "Ensemble"
  )

  # Select appropriate mapping
  mapping <- if (tolower(model_type) == "competing_risks" || tolower(model_type) == "cr") {
    cr_mapping
  } else {
    surv_mapping
  }

  # Apply mapping
  formatted <- vapply(model_names, function(x) {
    if (x %in% names(mapping)) {
      mapping[[x]]
    } else {
      # Fallback: convert underscores to spaces and title case
      gsub("_Model$", "", x) |>
        gsub("_", " ", x = _) |>
        tools::toTitleCase()
    }
  }, character(1), USE.NAMES = FALSE)

  names(formatted) <- model_names
  formatted
}


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
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}
