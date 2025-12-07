knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(ml4time2event)
library(dplyr)

bmt <- get_bmt_competing_risks_data() |>
  mutate(
    cause = if_else(status == 0L, NA_integer_, status),
    status = if_else(status == 0L, 0L, 1L)
  )

set.seed(2025)
train_rows <- sample.int(nrow(bmt), size = floor(0.7 * nrow(bmt)))
bmt_train <- bmt[train_rows, ]
bmt_test  <- bmt[-train_rows, ]

feature_cols <- c("sex", "d", "phase", "age", "source")

cr_task <- ml4t2e_task_cr(
  data = bmt_train,
  time = "ftime",
  status = "status",
  cause = "cause",
  features = feature_cols,
  time_units = "years"
)

cr_fit <- ml4t2e_fit(
  task = cr_task,
  models = c("cox", "fine_gray", "cr_random_forest"),
  ensemble = "auto",
  controls = list(times = seq(0, 120, length.out = 50))
)

ml4t2e_evaluate(
  cr_fit,
  metrics = c("c_index", "ibs"),
  include = c("ensemble", "cox", "fine_gray", "cr_random_forest")
)

cr_task_val <- ml4t2e_task_cr(
  data = bmt_test,
  time = "ftime",
  status = "status",
  cause = "cause",
  features = feature_cols,
  time_units = "years"
)

cr_preds <- predict(
  cr_fit,
  newdata = bmt_test,
  times = seq(0, 120, length.out = 50),
  type = "cif",
  include = c("ensemble", "cox", "fine_gray", "cr_random_forest")
)

ml4t2e_evaluate(
  cr_preds,
  task = cr_task_val,
  metrics = c("c_index", "ibs")
)

autoplot(cr_fit, include = c("ensemble", "cox"))

tmp_fit <- tempfile(fileext = ".rds")
ml4t2e_save(cr_fit, tmp_fit)
ml4t2e_load(tmp_fit)
