library(testthat)
library(ml4time2event)
library(survival)

test_that("Survival: Stacking Ensemble Integration", {
    # Survival Data
    lung_df <- survival::lung
    lung_df$status <- lung_df$status - 1
    lung_df$id <- seq_len(nrow(lung_df))
    lung_df <- na.omit(lung_df)

    # Use a large enough dataset to allow splitting
    task <- ml4t2e_task_surv(
        data = lung_df,
        time = "time",
        event = "status",
        features = c("age", "sex", "ph.ecog"), # minimal features
        id = "id"
    )

    # Fit with stacking
    # Force controls to speed up RF
    fit_stack <- ml4t2e_fit(keep_data = TRUE, 
        task = task,
        models = c("cox", "random_forest"),
        ensemble = "stack",
        controls = list(
            random_forest = list(ntree = 10),
            times = seq(0, 500, length.out = 10)
        ),
        seed = 42
    )

    expect_equal(fit_stack$ensemble_strategy, "stack")
    expect_s3_class(fit_stack$ensemble, "SurvivalEnsembler")

    # Check weights
    weights <- fit_stack$ensemble$weights
    expect_false(is.null(weights))
    expect_true(is.numeric(weights))
    expect_equal(sum(weights), 1, tolerance = 1e-4)
    expect_named(weights, c("cox", "random_forest"))

    # Predict
    preds <- predict(fit_stack, newdata = lung_df[1:5, ], times = 100)
    expect_s3_class(preds, "t2e_pred")

    ens_preds <- preds[preds$model == "ensemble", ]
    expect_gt(nrow(ens_preds), 0)
})

test_that("Competing Risk: Stacking Ensemble Integration", {
    # CR Data
    bmt <- get_bmt_competing_risks_data()
    bmt$id <- seq_len(nrow(bmt))

    # Prepare columns: status must be binary, cause must be integer code
    bmt$cause <- ifelse(bmt$status == 0, NA_integer_, bmt$status)
    bmt$binary_status <- ifelse(bmt$status == 0, 0L, 1L)

    task <- ml4t2e_task_cr(
        data = bmt,
        time = "ftime",
        status = "binary_status",
        cause = "cause",
        features = c("age", "sex"),
        id = "id"
    )

    # Fit with stacking
    fit_stack <- ml4t2e_fit(keep_data = TRUE, 
        task = task,
        models = c("cox", "cr_fine_gray"), # cr_random_forest might be slow?
        ensemble = "stack",
        controls = list(
            times = seq(0, 50, length.out = 10)
        ),
        seed = 99
    )

    expect_equal(fit_stack$ensemble_strategy, "stack")
    expect_s3_class(fit_stack$ensemble, "CompetingRiskEnsembler")

    weights <- fit_stack$ensemble$weights
    expect_false(is.null(weights))
    expect_equal(sum(weights), 1, tolerance = 1e-4)

    # Predict
    preds <- predict(fit_stack, newdata = bmt[1:5, ], times = 10)
    expect_s3_class(preds, "t2e_pred")

    ens_preds <- preds[preds$model == "ensemble", ]
    expect_gt(nrow(ens_preds), 0)
})
