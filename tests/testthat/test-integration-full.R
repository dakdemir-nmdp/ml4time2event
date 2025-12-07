test_that("Full pipeline integration test with save/load cycle", {
    skip_on_cran()

    # Use lung data
    lung <- survival::lung
    lung$status <- ifelse(lung$status == 2, 1, 0)
    lung <- na.omit(lung)

    # create pipeline
    pipe <- ml4t2e_pipeline(
        outcome = list(type = "survival", time = "time", event = "status"),
        models = c("cox"),
        ensemble = "none",
        recipe = function(x) {
            recipes::recipe(status ~ ., data = x) %>%
                recipes::step_impute_median(recipes::all_numeric_predictors())
        },
        resampling = rsample::vfold_cv(lung, v = 2)
    )

    # fit
    expect_no_error(pipe$fit(lung))

    # predict
    preds_orig <- pipe$predict(head(lung), times = c(100, 200))

    # save and load
    tmp <- tempfile(fileext = ".rds")
    ml4t2e_save(pipe, tmp)

    pipe_loaded <- ml4t2e_load(tmp)
    expect_true(inherits(pipe_loaded, "ml4t2e_pipeline"))

    # predict with loaded
    preds_loaded <- pipe_loaded$predict(head(lung), times = c(100, 200))

    expect_equal(preds_orig, preds_loaded)

    unlink(tmp)
})

test_that("ml4t2e_explain works with pipelines", {
    skip_on_cran()

    lung <- survival::lung
    lung$status <- ifelse(lung$status == 2, 1, 0)
    lung <- na.omit(lung)

    pipe <- ml4t2e_pipeline(
        outcome = list(type = "survival", time = "time", event = "status", features = c("age", "sex")),
        models = c("cox"),
        ensemble = "none"
    )
    pipe$fit(lung)

    # Explain
    expl <- ml4t2e_explain(pipe, newdata = head(lung, 10), times = 100)
    expect_s3_class(expl, "t2e_explain")
    expect_true("importance" %in% names(expl))
})

test_that("ml4t2e_calculate_shap works with pipelines", {
    skip_on_cran()
    skip_if_not_installed("fastshap")

    lung <- survival::lung
    lung$status <- ifelse(lung$status == 2, 1, 0)
    lung <- na.omit(lung)
    # reduce size for speed
    lung <- head(lung, 50)

    pipe <- ml4t2e_pipeline(
        outcome = list(type = "survival", time = "time", event = "status", features = c("age", "sex")),
        models = c("cox"),
        ensemble = "none"
    )
    pipe$fit(lung)

    # SHAP
    shap_res <- ml4t2e_calculate_shap(pipe, data = head(lung, 5), time_horizon = 100, nsim = 5)
    expect_s3_class(shap_res, "ml4t2e_shap")
})

test_that("ml4t2e_leaderboard works", {
    skip_on_cran()

    lung <- survival::lung
    lung$status <- ifelse(lung$status == 2, 1, 0)
    lung <- na.omit(lung)
    # reduce size
    lung <- head(lung, 40)

    pipe <- ml4t2e_pipeline(
        outcome = list(type = "survival", time = "time", event = "status"),
        models = c("cox"), # could add more but for speed just one
        ensemble = "none",
        resampling = rsample::vfold_cv(lung, v = 2)
    )
    pipe$fit(lung)

    lb <- ml4t2e_leaderboard(pipe)
    expect_s3_class(lb, "t2e_leaderboard")
    expect_true("mean" %in% names(lb))
    expect_true(nrow(lb) > 0)

    p <- autoplot(lb)
    expect_s3_class(p, "ggplot")
})
