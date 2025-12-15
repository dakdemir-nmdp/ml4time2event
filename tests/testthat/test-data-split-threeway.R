# Test: Data Splitting for Stacking and Conformal
#
# Tests that data splitting correctly handles:
# 1. No split (neither stacking nor conformal)
# 2. Stacking only
#  3. Conformal only
# 4. **CRITICAL**: 3-way split for both stacking + conformal (prevents data leakage)

test_that("No splitting when neither stacking nor conformal requested", {
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(100),
            event = rbinom(100, 1, 0.7),
            x1 = rnorm(100)
        ),
        time = "time",
        event = "event"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = FALSE,
        is_conformal = FALSE,
        conformal_ratio = 0.2
    )

    expect_equal(nrow(result$train_data), 100)
    expect_null(result$stack_data)
    expect_null(result$conf_data)
    expect_equal(result$split_type, "none")
})

test_that("Stacking-only split uses 20% for meta-learner", {
    set.seed(123)
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(100),
            event = rbinom(100, 1, 0.7),
            x1 = rnorm(100)
        ),
        time = "time",
        event = "event"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = TRUE,
        is_conformal = FALSE,
        conformal_ratio = 0.2
    )

    expect_equal(nrow(result$stack_data), 20) # 20% of 100
    expect_equal(nrow(result$train_data), 80) # Remaining 80%
    expect_null(result$conf_data)
    expect_equal(result$split_type, "stacking_only")
})

test_that("Conformal-only split uses specified ratio", {
    set.seed(124)
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(100),
            event = rbinom(100, 1, 0.7),
            x1 = rnorm(100)
        ),
        time = "time",
        event = "event"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = FALSE,
        is_conformal = TRUE,
        conformal_ratio = 0.15
    )

    expect_equal(nrow(result$conf_data), 15) # 15% of 100
    expect_equal(nrow(result$train_data), 85) # Remaining 85%
    expect_null(result$stack_data)
    expect_equal(result$split_type, "conformal_only")
})

test_that("3-way split prevents data leakage (CRITICAL TEST)", {
    set.seed(125)
    n <- 100
    task <- ml4t2e_task_surv(
        data = data.frame(
            id = 1:n,
            time = rexp(n),
            event = rbinom(n, 1, 0.7),
            x1 = rnorm(n)
        ),
        time = "time",
        event = "event",
        id = "id"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = TRUE,
        is_conformal = TRUE,
        conformal_ratio = 0.2 # This parameter is ignored; fixed 60-20-20 used
    )

    # Check split type
    expect_equal(result$split_type, "threeway")

    # Check split details are included
    expect_true(!is.null(result$split_details))
    expect_equal(result$split_details$n_total, n)

    # Check approximate 60-20-20 split
    n_train <- nrow(result$train_data)
    n_stack <- nrow(result$stack_data)
    n_conf <- nrow(result$conf_data)

    expect_equal(n_train, 60) # floor(100 * 0.6)
    expect_equal(n_stack, 20) # floor(100 * 0.2)
    expect_equal(n_conf, 20) # Remainder

    # CRITICAL: Check no overlap between sets
    train_ids <- result$train_data$id
    stack_ids <- result$stack_data$id
    conf_ids <- result$conf_data$id

    # No overlap between train and stack
    expect_equal(length(intersect(train_ids, stack_ids)), 0)

    # No overlap between train and conformal
    expect_equal(length(intersect(train_ids, conf_ids)), 0)

    # No overlap between stack and conformal (THIS IS THE KEY!)
    expect_equal(length(intersect(stack_ids, conf_ids)), 0)

    # All observations accounted for
    all_split_ids <- c(train_ids, stack_ids, conf_ids)
    expect_equal(length(all_split_ids), n)
    expect_equal(length(unique(all_split_ids)), n)
    expect_setequal(all_split_ids, 1:n)
})

test_that("3-way split uses stored indices correctly", {
    set.seed(126)
    n <- 50
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(n),
            event = rbinom(n, 1, 0.7),
            x1 = rnorm(n)
        ),
        time = "time",
        event = "event"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = TRUE,
        is_conformal = TRUE,
        conformal_ratio = 0.2
    )

    # Check that split_details contains indices
    expect_true(!is.null(result$split_details$train_indices))
    expect_true(!is.null(result$split_details$stack_indices))
    expect_true(!is.null(result$split_details$conformal_indices))

    # Verify indices match actual data
    expect_equal(
        task$data[result$split_details$train_indices, ],
        result$train_data
    )

    expect_equal(
        task$data[result$split_details$stack_indices, ],
        result$stack_data
    )

    expect_equal(
        task$data[result$split_details$conformal_indices, ],
        result$conf_data
    )
})

test_that("3-way split fails gracefully with insufficient data", {
    # Too small for 3-way split (needs ~10+ observations)
    small_task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(8),
            event = rbinom(8, 1, 0.7),
            x1 = rnorm(8)
        ),
        time = "time",
        event = "event"
    )

    expect_error(
        .split_data_for_training(
            task = small_task,
            is_stacking = TRUE,
            is_conformal = TRUE,
            conformal_ratio = 0.2
        ),
        "Insufficient data for 3-way split"
    )
})

test_that("Data split validation catches insufficient data", {
    # Test .validate_split_feasibility for various scenarios

    # Too small overall
    expect_error(
        .validate_split_feasibility(
            n_total = 5,
            stacking_requested = TRUE,
            conformal_requested = TRUE,
            conformal_ratio = 0.2
        ),
        "Dataset too small"
    )

    # Enough for conformal only
    expect_error(
        .validate_split_feasibility(
            n_total = 15,
            stacking_requested = FALSE,
            conformal_requested = TRUE,
            conformal_ratio = 0.2
        ),
        NA # Should not error
    )

    # Enough for stacking only
    expect_error(
        .validate_split_feasibility(
            n_total = 15,
            stacking_requested = TRUE,
            conformal_requested = FALSE,
            conformal_ratio = 0.2
        ),
        NA # Should not error
    )
})

test_that("3-way split message clearly indicates prevention of data leakage", {
    set.seed(127)
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = rexp(100),
            event = rbinom(100, 1, 0.7),
            x1 = rnorm(100)
        ),
        time = "time",
        event = "event"
    )

    result <- .split_data_for_training(
        task = task,
        is_stacking = TRUE,
        is_conformal = TRUE,
        conformal_ratio = 0.2
    )

    # Check that message mentions data leakage prevention
    expect_match(result$split_info, "3-way split")
    expect_match(result$split_info, "prevent.*leakage|leakage", ignore.case = TRUE)
    expect_match(result$split_info, "Train=")
    expect_match(result$split_info, "Stack=")
    expect_match(result$split_info, "Conformal=")
})

test_that("3-way split is reproducible with same seed", {
    n <- 100
    task <- ml4t2e_task_surv(
        data = data.frame(
            id = 1:n,
            time = rexp(n),
            event = rbinom(n, 1, 0.7),
            x1 = rnorm(n)
        ),
        time = "time",
        event = "event",
        id = "id"
    )

    set.seed(999)
    result1 <- .split_data_for_training(task, TRUE, TRUE, 0.2)

    set.seed(999)
    result2 <- .split_data_for_training(task, TRUE, TRUE, 0.2)

    # Should get identical splits with same seed
    expect_equal(result1$train_data$id, result2$train_data$id)
    expect_equal(result1$stack_data$id, result2$stack_data$id)
    expect_equal(result1$conf_data$id, result2$conf_data$id)
})
