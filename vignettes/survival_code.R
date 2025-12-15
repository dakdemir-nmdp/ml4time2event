knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(ml4time2event)
library(dplyr)
library(ggplot2)
library(survival)

lung <- get_lung_survival_data()

# Convert status from 1/2 to 0/1 format (1 = event, 0 = censored)
if (max(lung$status, na.rm = TRUE) > 1) {
  lung$status <- ifelse(lung$status == 2, 1L, 0L)
}

set.seed(2025)
train_rows <- sample.int(nrow(lung), size = floor(0.7 * nrow(lung)))
lung_train <- lung[train_rows, ]
lung_test <- lung[-train_rows, ]

feature_cols <- c(
  "age", "sex", "ph.ecog", "ph.karno", "pat.karno",
  "meal.cal", "wt.loss"
)

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

autoplot(surv_fit, include = c("ensemble", "cox"))

surv_importance <- ml4t2e_explain(
  surv_fit,
  newdata = lung_test,
  top_n = 5
)
autoplot(surv_importance)

tmp_fit <- tempfile(fileext = ".rds")
ml4t2e_save(surv_fit, tmp_fit)
ml4t2e_load(tmp_fit)
