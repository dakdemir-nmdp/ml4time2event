#' Data Splitting Utilities for Training
#'
#' Functions to handle data splitting for stacking, conformal calibration,
#' and their combinations. Provides validation and clear error messages.
#'
#' @keywords internal
#' @noRd

#' Split data for ensemble stacking and/or conformal calibration
#'
#' Handles three scenarios:
#' 1. Stacking only: Splits data into train + stack (meta-learner)
#' 2. Conformal only: Splits data into train + conformal calibration
#' 3. Both: Splits data into train + stack + conformal (holdout split 50/50)
#'
#' @param task Task object with $data field
#' @param is_stacking Logical; whether ensemble stacking is requested
#' @param is_conformal Logical; whether conformal calibration is requested
#' @param conformal_ratio Numeric (0,1); fraction for calibration
#' @return List with train_data, stack_data, conf_data, and metadata
#' @keywords internal
.split_data_for_training <- function(task, is_stacking, is_conformal, conformal_ratio = 0.2) {
    n_total <- nrow(task$data)

    # Validate inputs
    if (!is.logical(is_stacking) || !is.logical(is_conformal)) {
        rlang::abort(
            "Invalid arguments: is_stacking and is_conformal must be logical.",
            class = "invalid_argument_error"
        )
    }

    if (is_conformal && (conformal_ratio <= 0 || conformal_ratio >= 1)) {
        rlang::abort(
            glue::glue(
                "Invalid conformal_ratio: {conformal_ratio}. Must be between 0 and 1 (exclusive)."
            ),
            class = "invalid_argument_error"
        )
    }

    # No splitting needed
    if (!is_stacking && !is_conformal) {
        return(list(
            train_data = task$data,
            stack_data = NULL,
            conf_data = NULL,
            split_type = "none",
            split_info = glue::glue("Using all {n_total} observations for training.")
        ))
    }

    # Both stacking and conformal - CRITICAL: Use 3-way split to prevent data leakage
    if (is_stacking && is_conformal) {
        # IMPORTANT: We need THREE separate, non-overlapping sets:
        # 1. Training set (60%) - for fitting base models
        # 2. Stacking set (20%) - for training meta-learner
        # 3. Conformal set (20%) - for calibrating prediction intervals
        #
        # Using the same data for stacking and conformal creates leakage because:
        # - Meta-learner learns weights from holdout
        # - Conformal bands are computed on same holdout
        # - Result: Overly optimistic confidence intervals

        # Use a fixed 60-20-20 split when both are requested
        # (ignoring conformal_ratio parameter to ensure proper 3-way split)
        train_ratio <- 0.60
        stack_ratio <- 0.20

        num_train <- floor(n_total * train_ratio)
        num_stack <- floor(n_total * stack_ratio)
        num_conf <- n_total - num_train - num_stack # Remainder for conformal

        # Validate we have enough data for 3-way split
        if (num_train < 5 || num_stack < 2 || num_conf < 2) {
            rlang::abort(
                glue::glue(
                    "Insufficient data for 3-way split (training + stacking + conformal).\\n",
                    "Total observations: {n_total}\\n",
                    "Required split: Train={num_train} (min 5), Stack={num_stack} (min 2), Conformal={num_conf} (min 2)\\n",
                    "Minimum total required: ~10 observations\\n",
                    "Suggestion: Use a larger dataset, or disable either stacking or conformal calibration."
                ),
                class = "insufficient_data_error"
            )
        }

        # Create random permutation for 3-way split
        all_indices <- sample(seq_len(n_total))

        # Assign to three SEPARATE, NON-OVERLAPPING sets
        train_indices <- all_indices[1:num_train]
        stack_indices <- all_indices[(num_train + 1):(num_train + num_stack)]
        conf_indices <- all_indices[(num_train + num_stack + 1):n_total]

        # Sanity check: no overlap
        stopifnot(
            length(intersect(train_indices, stack_indices)) == 0,
            length(intersect(train_indices, conf_indices)) == 0,
            length(intersect(stack_indices, conf_indices)) == 0,
            length(train_indices) + length(stack_indices) + length(conf_indices) == n_total
        )

        return(list(
            train_data = task$data[train_indices, ],
            stack_data = task$data[stack_indices, ],
            conf_data = task$data[conf_indices, ],
            split_type = "threeway",
            split_info = glue::glue(
                "3-way split to prevent data leakage: ",
                "Train={length(train_indices)} ({round(100*length(train_indices)/n_total)}%), ",
                "Stack={length(stack_indices)} ({round(100*length(stack_indices)/n_total)}%), ",
                "Conformal={length(conf_indices)} ({round(100*length(conf_indices)/n_total)}%)"
            ),
            # Include split details for transparency
            split_details = list(
                n_total = n_total,
                n_train = length(train_indices),
                n_stack = length(stack_indices),
                n_conformal = length(conf_indices),
                train_indices = train_indices,
                stack_indices = stack_indices,
                conformal_indices = conf_indices
            )
        ))
    }

    # Stacking only
    if (is_stacking) {
        ratio <- 0.2 # Fixed 20% for stacking
        num_stack <- floor(n_total * ratio)

        if (num_stack < 2) {
            rlang::abort(
                glue::glue(
                    "Insufficient data for stacking.\n",
                    "Total observations: {n_total}\n",
                    "Stacking split would have: {num_stack} observations\n",
                    "Minimum required: 2\n",
                    "Suggestion: Use a larger dataset or disable stacking."
                ),
                class = "insufficient_data_error"
            )
        }

        stack_indices <- sample(seq_len(n_total), size = num_stack)

        return(list(
            train_data = task$data[-stack_indices, ],
            stack_data = task$data[stack_indices, ],
            conf_data = NULL,
            split_type = "stacking_only",
            split_info = glue::glue(
                "Using 20% ({length(stack_indices)} obs) for stacking meta-learner, ",
                "{nrow(task$data[-stack_indices,])} for training."
            )
        ))
    }

    # Conformal only
    if (is_conformal) {
        ratio <- conformal_ratio
        num_conf <- floor(n_total * ratio)

        if (num_conf < 2) {
            rlang::abort(
                glue::glue(
                    "Insufficient data for conformal calibration.\n",
                    "Total observations: {n_total}\n",
                    "Calibration split would have: {num_conf} observations\n",
                    "Requested ratio: {ratio*100}%\n",
                    "Minimum required: 2\n",
                    "Suggestion: Reduce conformal_calibration ratio or use more data."
                ),
                class = "insufficient_data_error"
            )
        }

        conf_indices <- sample(seq_len(n_total), size = num_conf)

        return(list(
            train_data = task$data[-conf_indices, ],
            stack_data = NULL,
            conf_data = task$data[conf_indices, ],
            split_type = "conformal_only",
            split_info = glue::glue(
                "Using {ratio*100}% ({length(conf_indices)} obs) for conformal calibration, ",
                "{nrow(task$data[-conf_indices,])} for training."
            )
        ))
    }
}

#' Validate data split feasibility
#'
#' Checks if requested split ratios are feasible given data size.
#' Returns informative error if not.
#'
#' @param n_total Integer; total number of observations
#' @param stacking_requested Logical
#' @param conformal_requested Logical
#' @param conformal_ratio Numeric (0,1)
#' @return Invisible TRUE if valid, otherwise aborts with error
#' @keywords internal
.validate_split_feasibility <- function(n_total, stacking_requested, conformal_requested, conformal_ratio = 0.2) {
    min_required <- 10 # Absolute minimum for any analysis

    if (n_total < min_required) {
        rlang::abort(
            glue::glue(
                "Dataset too small for time-to-event analysis.\n",
                "Observations: {n_total}\n",
                "Minimum required: {min_required}\n",
                "Suggestion: Collect more data before model fitting."
            ),
            class = "insufficient_data_error"
        )
    }

    # Check if splits are feasible
    if (stacking_requested && conformal_requested) {
        holdout_size <- floor(n_total * conformal_ratio)
        train_size <- n_total - holdout_size

        if (train_size < 5) {
            rlang::abort(
                glue::glue(
                    "Requested splits leave too few training observations.\n",
                    "Total: {n_total}, Holdout: {holdout_size}, Remaining for train: {train_size}\n",
                    "Minimum training size: 5\n",
                    "Suggestion: Reduce conformal_calibration ratio or disable one of stacking/conformal."
                ),
                class = "insufficient_data_error"
            )
        }
    } else if (conformal_requested) {
        cal_size <- floor(n_total * conformal_ratio)
        train_size <- n_total - cal_size

        if (train_size < 5 || cal_size < 2) {
            rlang::abort(
                glue::glue(
                    "Conformal split infeasible.\n",
                    "Total: {n_total}, Ratio: {conformal_ratio}\n",
                    "Would create: Train={train_size}, Calibration={cal_size}\n",
                    "Minimum: Train=5, Calibration=2\n",
                    "Suggestion: Reduce conformal_calibration ratio."
                ),
                class = "insufficient_data_error"
            )
        }
    } else if (stacking_requested) {
        stack_size <- floor(n_total * 0.2)
        train_size <- n_total - stack_size

        if (train_size < 5 || stack_size < 2) {
            rlang::abort(
                glue::glue(
                    "Stacking split infeasible.\n",
                    "Total: {n_total}\n",
                    "Would create: Train={train_size}, Stack={stack_size}\n",
                    "Minimum: Train=5, Stack=2\n",
                    "Suggestion: Use a larger dataset or disable stacking."
                ),
                class = "insufficient_data_error"
            )
        }
    }

    invisible(TRUE)
}

#' Check if data splitting is stratifiable
#'
#' For survival/competing risks, checks if event rates are sufficient
#' to support stratified splitting.
#'
#' @param data Data frame
#' @param event_col Character; name of event column
#' @param min_events Integer; minimum events required
#' @return Logical; TRUE if stratification is recommended
#' @keywords internal
.check_stratification_feasibility <- function(data, event_col, min_events = 10) {
    if (is.null(event_col) || !event_col %in% names(data)) {
        return(FALSE)
    }

    events <- data[[event_col]]
    n_events <- sum(events > 0, na.rm = TRUE)

    # Need at least min_events to stratify
    n_events >= min_events
}
