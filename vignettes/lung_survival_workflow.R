## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----load-libraries-----------------------------------------------------------
library(ml4time2event)
library(dplyr)
library(ggplot2)
library(recipes)
library(rsample)

lung <- get_lung_survival_data()

# Convert status from 1/2 to 0/1 format (1 = event, 0 = censored)
# The lung data uses 1 = censored, 2 = death, but we need 0 = censored, 1 = death
if (max(lung$status, na.rm = TRUE) > 1) {
  lung$status <- ifelse(lung$status == 2, 1L, 0L)
}

glimpse(lung)

## ----load-data----------------------------------------------------------------
table(lung$status, useNA = "ifany")



## ----split-data---------------------------------------------------------------
set.seed(2025)
split <- initial_split(lung, prop = 0.75, strata = "status")
lung_train <- training(split)
lung_test  <- testing(split)

feature_cols <- c("age", "sex", "ph.ecog", "ph.karno", "pat.karno",
                  "meal.cal", "wt.loss")


## ----task-and-fit-------------------------------------------------------------
surv_task <- ml4t2e_task_surv(
  data = lung_train,
  time = "time",
  event = "status",
  features = feature_cols,
  time_units = "days"
)

surv_fit <- ml4t2e_fit(
  task = surv_task,
  models = c("cox", "random_forest", "gbm"),
  ensemble = "auto",
  controls = list(times = seq(0, 1000, length.out = 60))
)

ml4t2e_evaluate(
  surv_fit,
  metrics = c("c_index", "ibs"),
  include = c("ensemble", "cox", "random_forest", "gbm")
)


## ----validation---------------------------------------------------------------
surv_task_val <- ml4t2e_task_surv(
  data = lung_test,
  time = "time",
  event = "status",
  features = feature_cols,
  time_units = "days"
)

surv_preds <- predict(
  surv_fit,
  newdata = lung_test,
  times = seq(0, 1000, length.out = 60),
  type = "survival",
  include = c("ensemble", "cox", "random_forest", "gbm")
)

ml4t2e_evaluate(
  surv_preds,
  task = surv_task_val,
  metrics = c("c_index", "ibs")
)


## ----visualise---------------------------------------------------------------
autoplot(surv_fit, include = c("ensemble", "cox"))


## ----importance--------------------------------------------------------------
surv_importance <- ml4t2e_explain(
  surv_fit,
  newdata = lung_test,
  top_n = 6
)
autoplot(surv_importance)


## ----pipeline-----------------------------------------------------------------
surv_pipeline <- ml4t2e_pipeline(
  outcome = list(type = "survival", time = "time", event = "status"),
  models = c("cox", "random_forest"),
  ensemble = "auto",
  recipe = recipe(status ~ ., data = lung_train) |>
    step_impute_median(all_numeric_predictors()) |>
    step_dummy(all_nominal_predictors()),
  resampling = vfold_cv(lung_train, v = 3)
)

surv_pipeline$fit(lung_train)
surv_pipeline$summary()

surv_pipeline$evaluate(
  newdata = lung_test,
  metrics = c("c_index", "ibs"),
  include = "ensemble"
)


## ----persistence--------------------------------------------------------------
tmp_fit <- tempfile(fileext = ".rds")
ml4t2e_save(surv_fit, tmp_fit)
ml4t2e_load(tmp_fit)

tmp_pipe <- tempfile(fileext = ".rds")
saveRDS(surv_pipeline, tmp_pipe)
readRDS(tmp_pipe)$summary()
