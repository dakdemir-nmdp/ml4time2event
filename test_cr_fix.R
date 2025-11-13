#!/usr/bin/env Rscript

# Test script to verify the competing risks fix
devtools::load_all()
library(dplyr)

# Load and prepare data
bmt_raw <- get_bmt_competing_risks_data()
bmt_df <- bmt_raw %>%
  mutate(
    cause = if_else(status == 0L, NA_integer_, status),
    status = if_else(status == 0L, 0L, 1L)
  )

# Create train split
set.seed(123)
train_rows <- sample.int(nrow(bmt_df), size = floor(0.7 * nrow(bmt_df)))
bmt_train <- bmt_df[train_rows, ]

# Create task
cr_task <- ml4t2e_task_cr(
  data = bmt_train,
  time = "ftime",
  status = "status",
  cause = "cause",
  features = c("sex", "d", "phase", "age", "source"),
  time_units = "years"
)

# Fit models (this was failing before the fix)
cr_fit <- ml4t2e_fit(
  task = cr_task,
  models = c("cox", "fine_gray", "cr_random_forest"),
  ensemble = "auto",
  controls = list(times = seq(0, 120, length.out = 60))
)

cat("\n=== SUCCESS! ===\n")
print(cr_fit)
cat("\nThe fix is working correctly!\n")
