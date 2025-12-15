# Test: Input Validation
#
# Tests that validation functions properly catch invalid inputs
# and provide clear error messages.

test_that("ml4t2e_task_surv() validates negative times", {
    # Create data with negative times
    bad_data <- data.frame(
        time = c(-10, -5, 100, 200),
        event = c(1, 1, 0, 1),
        x1 = c(1, 2, 3, 4)
    )

    expect_error(
        ml4t2e_task_surv(bad_data, time = "time", event = "event"),
        "negative"
    )
})

test_that("ml4t2e_task_surv() validates NA in time column", {
    bad_data <- data.frame(
        time = c(NA, 5, 100, 200),
        event = c(1, 1, 0, 1),
        x1 = c(1, 2, 3, 4)
    )

    expect_error(
        ml4t2e_task_surv(bad_data, time = "time", event = "event"),
        "NA value"
    )
})

test_that("ml4t2e_task_surv() validates empty data", {
    empty_data <- data.frame(
        time = numeric(0),
        event = integer(0),
        x1 = numeric(0)
    )

    expect_error(
        ml4t2e_task_surv(empty_data, time = "time", event = "event"),
        "empty"
    )
})

test_that("ml4t2e_task_surv() aborts if all observations are censored", {
    censored_data <- data.frame(
        time = c(10, 20, 30),
        event = c(0, 0, 0), # All censored
        x1 = c(1, 2, 3)
    )

    expect_error(
        ml4t2e_task_surv(censored_data, time = "time", event = "event"),
        "censored"
    )
})

test_that("ml4t2e_task_surv() warns with very few events", {
    few_events_data <- data.frame(
        time = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
        event = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), # Only 1 event
        x1 = 1:10
    )

    expect_warning(
        ml4t2e_task_surv(few_events_data, time = "time", event = "event"),
        "Very few events"
    )
})

test_that("ml4t2e_task_surv() validates duplicate IDs if provided", {
    dup_id_data <- data.frame(
        id = c(1, 1, 2, 3), # Duplicate ID
        time = c(10, 20, 30, 40),
        event = c(1, 1, 0, 1),
        x1 = c(1, 2, 3, 4)
    )

    expect_error(
        ml4t2e_task_surv(dup_id_data, time = "time", event = "event", id = "id"),
        "Duplicate IDs"
    )
})

test_that("ml4t2e_task_cr() validates negative times", {
    bad_cr_data <- data.frame(
        time = c(-10, 5, 100),
        status = c(1, 1, 0),
        cause = c("A", "B", NA),
        x1 = c(1, 2, 3)
    )

    expect_error(
        ml4t2e_task_cr(bad_cr_data, time = "time", status = "status", cause = "cause"),
        "negative"
    )
})

test_that("ml4t2e_task_cr() validates NA in time column", {
    bad_cr_data <- data.frame(
        time = c(NA, 5, 100),
        status = c(1, 1, 0),
        cause = c("A", "B", NA),
        x1 = c(1, 2, 3)
    )

    expect_error(
        ml4t2e_task_cr(bad_cr_data, time = "time", status = "status", cause = "cause"),
        "NA value"
    )
})

test_that("ml4t2e_fit() validates task structure", {
    # Test with non-task object
    expect_error(
        ml4t2e_fit(task = list(data = data.frame()), models = "cox"),
        "Invalid task object"
    )

    # Test with NULL task
    expect_error(
        ml4t2e_fit(task = NULL, models = "cox"),
        "Task is NULL"
    )
})

test_that("ml4t2e_fit() validates model names against registry", {
    # Create valid task
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = c(10, 20, 30, 40, 50),
            event = c(1, 1, 0, 1, 0),
            x1 = c(1, 2, 3, 4, 5)
        ),
        time = "time",
        event = "event"
    )

    # Test with invalid model name
    expect_error(
        ml4t2e_fit(task, models = "nonexistent_model", keep_data = TRUE),
        "Unknown model"
    )

    # Test with mix of valid and invalid models
    expect_error(
        ml4t2e_fit(task, models = c("cox", "fake_model"), keep_data = TRUE),
        "Unknown model"
    )
})

test_that("ml4t2e_fit() validates ensemble strategy", {
    task <- ml4t2e_task_surv(
        data = data.frame(
            time = c(10, 20, 30, 40, 50),
            event = c(1, 1, 0, 1, 0),
            x1 = c(1, 2, 3, 4, 5)
        ),
        time = "time",
        event = "event"
    )

    # Valid ensemble strategies should work
    expect_error(
        ml4t2e_fit(task, models = "cox", ensemble = "simple", keep_data = TRUE),
        NA
    )

    expect_error(
        ml4t2e_fit(task, models = "cox", ensemble = FALSE, keep_data = TRUE),
        NA
    )
})

test_that(".validate_conformal_alpha() validates alpha range", {
    # Alpha too low
    expect_error(
        .validate_conformal_alpha(alpha = 0, allow_null = FALSE),
        "between 0 and 1"
    )

    # Alpha too high
    expect_error(
        .validate_conformal_alpha(alpha = 1, allow_null = FALSE),
        "between 0 and 1"
    )

    # Alpha negative
    expect_error(
        .validate_conformal_alpha(alpha = -0.1, allow_null = FALSE),
        "between 0 and 1"
    )

    # Valid alpha
    expect_error(
        .validate_conformal_alpha(alpha = 0.1, allow_null = FALSE),
        NA
    )
})

test_that(".validate_times() validates time vector", {
    # Non-numeric times
    expect_error(
        .validate_times(times = c("a", "b", "c"), allow_null = FALSE),
        "must be numeric"
    )

    # NA in times
    expect_error(
        .validate_times(times = c(1, NA, 3), allow_null = FALSE),
        "NA values"
    )

    # Negative times
    expect_error(
        .validate_times(times = c(-1, 0, 1), allow_null = FALSE, min_val = 0),
        "below minimum"
    )

    # Valid times
    expect_error(
        .validate_times(times = c(0, 1, 2, 3), allow_null = FALSE, min_val = 0),
        NA
    )
})

test_that(".validate_data_quality() catches all data issues", {
    # Missing required column
    expect_error(
        .validate_data_quality(
            data = data.frame(x = 1:5),
            required_cols = c("time", "event")
        ),
        "Required column.*missing"
    )

    # Non-data.frame input
    expect_error(
        .validate_data_quality(
            data = list(time = 1:5),
            required_cols = c("time")
        ),
        "must be a data.frame"
    )

    # Empty data.frame
    expect_error(
        .validate_data_quality(
            data = data.frame(time = numeric(0), event = integer(0)),
            required_cols = c("time", "event")
        ),
        "empty"
    )
})

test_that("validation error messages are informative", {
    # Check that error messages include context
    bad_data <- data.frame(
        time = c(-1, 5, 10),
        event = c(1, 1, 0),
        x1 = 1:3
    )

    err <- tryCatch(
        ml4t2e_task_surv(bad_data, time = "time", event = "event"),
        error = function(e) e
    )

    # Error message should mention the context
    expect_match(err$message, "ml4t2e_task_surv")
    expect_match(err$message, "negative")
})

test_that("validation passes for valid inputs", {
    # Valid survival task
    valid_data <- data.frame(
        time = c(10, 20, 30, 40, 50),
        event = c(1, 1, 0, 1, 0),
        x1 = c(1, 2, 3, 4, 5),
        x2 = c(2.1, 3.2, 4.3, 5.4, 6.5)
    )

    task <- expect_error(
        ml4t2e_task_surv(valid_data, time = "time", event = "event"),
        NA
    )

    expect_s3_class(task, "t2e_task_surv")
    expect_s3_class(task, "t2e_task")

    # Valid fit
    fit <- expect_error(
        ml4t2e_fit(task, models = "cox", keep_data = TRUE),
        NA
    )

    expect_s3_class(fit, "t2e_fit")
})

test_that("validation allows NULL where appropriate", {
    # .validate_conformal_alpha allows NULL when allow_null = TRUE
    expect_error(
        .validate_conformal_alpha(alpha = NULL, allow_null = TRUE),
        NA
    )

    # .validate_times allows NULL when allow_null = TRUE
    expect_error(
        .validate_times(times = NULL, allow_null = TRUE),
        NA
    )
})
