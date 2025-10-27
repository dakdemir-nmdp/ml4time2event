# Pipeline Quickstart --------------------------------------------------------

library(ml4time2event)
library(dplyr)

if (interactive() && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
}

# randomForestSRC can segfault on some systems when using multiple cores; force single core.
Sys.setenv(RF_CORES = "1", OMP_NUM_THREADS = "1")
options(rf.cores = 1)

# -----------------------------------------------------------------------------
# Survival pipeline
# -----------------------------------------------------------------------------

lung_df <- get_lung_survival_data()

lung_df %>%
  dplyr::select(time, status, age, sex, ph.ecog, performance_good) %>%
  dplyr::slice_head(n = 5)

# add missingness to demonstrate missing data handling (at random) in all variables including outcomes

set.seed(42)

for (col in names(lung_df)) {
  missing_indices <- sample(
    seq_len(nrow(lung_df)),
    size = floor(0.1 * nrow(lung_df))
  )
  lung_df[missing_indices, col] <- NA
}



surv_pipeline <- ml4t2e_fit_pipeline(
  data = lung_df,
  analysis_type = "survival",
  timevar = "time",
  eventvar = "status",
  models = c("glmnet", "coxph"),
  include_rf = FALSE,
  prediction_times = seq(0, 1000, length.out = 50)
)

print(surv_pipeline)

length(surv_pipeline$processed_expvars)
head(surv_pipeline$processed_expvars)

lung_holdout <- lung_df %>% slice(1:5)

surv_predictions <- predict(
  surv_pipeline,
  newdata = lung_holdout,
  new_times = seq(0, 730, length.out = 25)
)

names(surv_predictions$predictions)
head(surv_predictions$predictions$NewProbs)

# Persistence example (optional)
save_path <- tempfile("lung_pipeline_", fileext = ".rds")
ml4t2e_save_pipeline(surv_pipeline, save_path)
restored_pipeline <- ml4t2e_load_pipeline(save_path)
print(restored_pipeline)

restored_predictions <- predict(
  restored_pipeline,
  newdata = lung_holdout,
  new_times = seq(0, 730, length.out = 25)
)

all.equal(
  surv_predictions$predictions$NewProbs,
  restored_predictions$predictions$NewProbs
)

surv_plot_inputs <- c(
  lapply(
    surv_predictions$predictions$ModelPredictions,
    function(probs_mat) list(
      Probs = probs_mat,
      Times = surv_predictions$predictions$NewTimes
    )
  ),
  list(
    Ensemble = list(
      Probs = surv_predictions$predictions$NewProbs,
      Times = surv_predictions$predictions$NewTimes
    )
  )
)

surv_plot <- plot_survival_curves(
  predictions = surv_plot_inputs,
  patients_to_plot = 1:3,
  title = "Lung Cancer Survival Predictions"
)

print(surv_plot)

# -----------------------------------------------------------------------------
# Competing-risks pipeline
# -----------------------------------------------------------------------------

bmt_df <- get_bmt_competing_risks_data()

bmt_df %>%
  dplyr::select(ftime, status, age, sex, d, phase) %>%
  dplyr::slice_head(n = 5)

# add missingness to demonstrate missing data handling (at random) in all variables including outcomes
set.seed(42)
for (col in names(bmt_df)) {
  missing_indices <- sample(
    seq_len(nrow(bmt_df)),
    size = floor(0.1 * nrow(bmt_df))
  )
  bmt_df[missing_indices, col] <- NA
}


cr_pipeline <- ml4t2e_fit_pipeline(
  data = bmt_df,
  analysis_type = "competing_risks",
  timevar = "ftime",
  eventvar = "status",
  models = c("FG", "cox"),
  include_rf = FALSE,
  prediction_times = seq(0, 150, length.out = 40)
)

print(cr_pipeline)

bmt_holdout <- bmt_df %>% slice(1:4)

cr_predictions <- predict(
  cr_pipeline,
  newdata = bmt_holdout,
  ensemble_method = "average"
)

names(cr_predictions$predictions)
head(cr_predictions$predictions$NewProbs)

cr_plot_inputs <- c(
  lapply(
    cr_predictions$predictions$ModelPredictions,
    function(probs_mat) list(
      CIFs = probs_mat,
      Times = cr_predictions$predictions$NewTimes
    )
  ),
  list(
    Ensemble = list(
      CIFs = cr_predictions$predictions$NewProbs,
      Times = cr_predictions$predictions$NewTimes
    )
  )
)

cr_plot <- plot_cif_curves(
  predictions = cr_plot_inputs,
  patients_to_plot = 1:3,
  event_label = "Event of Interest",
  title = "Competing-Risks CIF Predictions"
)

print(cr_plot)


# model persistance for CR pipeline 
save_path_cr <- tempfile("bmt_pipeline_", fileext = ".rds")
ml4t2e_save_pipeline(cr_pipeline, save_path_cr)
restored_cr_pipeline <- ml4t2e_load_pipeline(save_path_cr)
print(restored_cr_pipeline)
restored_cr_predictions <- predict(
  restored_cr_pipeline,
  newdata = bmt_holdout,
  ensemble_method = "average"
)
all.equal(
  cr_predictions$predictions$NewProbs,
  restored_cr_predictions$predictions$NewProbs
)
