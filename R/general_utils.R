VariableProfile <- function(data, expvars) {
  varprofile <- vector(mode = "list", length = length(expvars))
  names(varprofile) <- expvars
  for (vari in expvars) {
    if (vari %in% colnames(data)) { # Check if variable exists in data
      col_data <- data[[vari]] # Use [[ ]] for safer column access
      if (is.factor(col_data)) {
        varprofile[[vari]] <- table(col_data, useNA = "ifany") # Include NAs in table
      } else if (is.numeric(col_data)) {
        varprofile[[vari]] <- c(min = min(col_data, na.rm = TRUE), max = max(col_data, na.rm = TRUE)) # Add names
      } else if (is.character(col_data)) {
        varprofile[[vari]] <- table(col_data, useNA = "ifany") # Include NAs
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
  tryCatch(
    {
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
          return("") # Root node has no rule
        }
        # For other internal nodes, return a placeholder
        # A full implementation would trace the path from root to this node
        return(paste("Internal node", node_id, "path"))
      }

      return("")
    },
    error = function(e) {
      # Fallback if internal function is not available
      return(paste("Rule for node", node_id))
    }
  )
}
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

#' Reshape Predictions to Matrix
#'
#' Optimized helper to convert tidy predictions to matrix for score computation.
#'
#' @param preds Tidy prediction data frame.
#' @param data Data frame (to align rows).
#' @param time_grid Numeric vector of times.
#' @param val_col Character, name of value column ("surv" or "cif").
#' @param id_col Optional character, name of ID column.
#'
#' @return Matrix (n x m).
#' @keywords internal
ml4t2e_reshape_preds_to_matrix <- function(preds, data, time_grid, val_col, id_col = NULL) {
  # Input validation
  if (!all(c("time", val_col) %in% names(preds))) {
    stop(sprintf("preds must contain 'time' and '%s' columns", val_col))
  }

  # Handle IDs
  if (is.null(id_col)) {
    data_ids <- seq_len(nrow(data))
    if (!"id" %in% names(preds)) {
      stop("preds must have an 'id' column")
    }
  } else {
    data_ids <- data[[id_col]]
  }

  # Standardize times to character for consistent column naming
  # Use high precision formatting to avoid collision
  fmt_time <- function(t) sprintf("%.10g", t)

  preds_clean <- preds
  # Initialize to avoid lint NOTE
  time_str <- NULL
  preds_clean$time_str <- fmt_time(as.numeric(preds$time))
  time_grid_str <- fmt_time(as.numeric(time_grid))

  # Filter to only times in the grid (ignores extra predictions)
  preds_clean <- preds_clean[preds_clean$time_str %in% time_grid_str, ]

  wide_df <- preds_clean %>%
    dplyr::select(id, time_str, dplyr::all_of(val_col)) %>%
    tidyr::pivot_wider(
      names_from = "time_str",
      values_from = dplyr::all_of(val_col),
      values_fn = mean # Handle duplicates if any (shouldn't be)
    )

  # Create a template dataframe to ensure all IDs and Times exist
  col_order <- time_grid_str

  # Ensure all target columns exist
  missing_cols <- setdiff(col_order, names(wide_df))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) wide_df[[mc]] <- NA_real_
  }

  # Join with original IDs to ensure correct row order and missing IDs are handled
  df_map <- data.frame(id = data_ids, stringsAsFactors = FALSE)

  # Harmonize ID types
  if (!identical(class(df_map$id), class(wide_df$id))) {
    wide_df$id <- as(wide_df$id, class(df_map$id))
  }

  df_final <- dplyr::left_join(df_map, wide_df, by = "id")

  # Extract matrix
  mat <- as.matrix(df_final[, col_order])

  # Ensure numeric (matrix might become character if NAs handled poorly?)
  mode(mat) <- "numeric"

  if (ncol(mat) != length(time_grid)) {
    stop("Reshaped matrix has wrong number of columns")
  }

  mat
}
