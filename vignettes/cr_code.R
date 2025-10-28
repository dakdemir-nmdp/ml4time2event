## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----load-libraries-----------------------------------------------------------
library(ml4time2event)
library(survival)
library(dplyr)
library(ggplot2)

if (!requireNamespace("cmprsk", quietly = TRUE)) {
  stop("Install the 'cmprsk' package to run this vignette.")
}
library(cmprsk)
if (interactive() && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all()
}


## ----load-data----------------------------------------------------------------
bmt_data <- ml4time2event::get_bmt_competing_risks_data()

glimpse(bmt_data)

summary(bmt_data[c("ftime", "status", "age")])


## ----explore-events-----------------------------------------------------------
status_counts <- dplyr::count(bmt_data, status)
status_counts

total_n <- nrow(bmt_data)
event_profile <- data.frame(
  status = status_counts$status,
  count = status_counts$n,
  proportion = round(status_counts$n / total_n, 3)
)
event_profile


## ----variable-profile---------------------------------------------------------
candidate_expvars <- c("sex", "d", "phase", "age", "source")
var_profile <- VariableProfile(bmt_data, candidate_expvars)
var_profile$Summary


## ----cif-phase----------------------------------------------------------------
cif_phase <- cuminc(
  ftime = bmt_data$ftime,
  fstatus = bmt_data$status,
  group = bmt_data$phase
)

timepoint_12mo <- 12
cif_summary <- data.frame()

for (phase_name in levels(bmt_data$phase)) {
  cif_name_relapse <- paste(phase_name, "1")
  cif_name_trm <- paste(phase_name, "2")

  if (cif_name_relapse %in% names(cif_phase)) {
    cif_obj <- cif_phase[[cif_name_relapse]]
    idx <- which.min(abs(cif_obj$time - timepoint_12mo))
    cif_summary <- rbind(
      cif_summary,
      data.frame(
        Phase = phase_name,
        Event = "Relapse",
        CIF_12_months = round(cif_obj$est[idx], 3)
      )
    )
  }

  if (cif_name_trm %in% names(cif_phase)) {
    cif_obj <- cif_phase[[cif_name_trm]]
    idx <- which.min(abs(cif_obj$time - timepoint_12mo))
    cif_summary <- rbind(
      cif_summary,
      data.frame(
        Phase = phase_name,
        Event = "TRM",
        CIF_12_months = round(cif_obj$est[idx], 3)
      )
    )
  }
}

cif_summary


## ----split-data---------------------------------------------------------------
set.seed(123)
split_data <- t2edata_split(bmt_data, prop = 0.7)
train_data <- split_data$Train
test_data <- split_data$Test

data.frame(
  dataset = c("Training", "Test"),
  observations = c(nrow(train_data), nrow(test_data))
)


## ----train-models-------------------------------------------------------------
timevar <- "ftime"
eventvar <- "status"
expvars <- c("sex", "d", "phase", "age", "source")
primary_event_code <- 1L
competing_event_codes <- 2L

data.frame(
  parameter = c(
    "Time variable",
    "Event variable",
    "Primary event",
    "Competing event(s)",
    "Predictors"
  ),
  value = c(
    timevar,
    eventvar,
    primary_event_code,
    paste(competing_event_codes, collapse = ", "),
    paste(expvars, collapse = ", ")
  )
)

models <- RunCRModels(
  datatrain = train_data,
  ExpVars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  models = c("cox", "FG", "xgboost", "gam", "bart", "rulefit", "survreg"),
  ntreeRF = 300,
  varsel = FALSE
)

# Individual model components now include:
# - For Cox model: ntimes, alpha (elastic net), nfolds (cross-validation)
# - For XGBoost: ntimes (prediction grid), verbose (progress reporting)
# - Consistent time grid specification across all competing risks models

models


## ----generate-predictions-----------------------------------------------------
# Note: With the updated API, time grids are now automatically managed by individual
# models. The ntimes parameter controls prediction grid resolution (default: 50).
# Individual model predictions include a 'Times' attribute for easy access.

ensemble_predictions <- PredictCRModels(
  models = models,
  newdata = test_data,
  new_times = seq(0, max(test_data$ftime), length.out = 50),
  ensemble_method = "average"
)

individual_predictions <- ensemble_predictions$ModelPredictions
ensemble_cif <- ensemble_predictions$NewProbs
common_times <- ensemble_predictions$NewTimes
models_used <- ensemble_predictions$models_used

data.frame(
  metric = c(
    "Subjects",
    "Prediction time points",
    "Models in ensemble"
  ),
  value = c(
    ncol(ensemble_cif),
    length(common_times),
    paste(models_used, collapse = ", ")
  )
)


## ----prepare-predictions------------------------------------------------------
pretty_cr_model_name <- function(x) {
  mapping <- c(
    RF_Model = "Random Forest",
    RF_Model2 = "Random Forest (Top Vars)",
    FG_Model = "Fine-Gray",
    BART_Model = "BART",
    Cox_Model = "Cox",
    rulefit_Model = "RuleFit",
    xgboost_Model = "XGBoost",
    gam_Model = "GAM",
    survreg_Model = "SurvReg"
  )
  if (x %in% names(mapping)) mapping[[x]] else tools::toTitleCase(gsub("_", " ", gsub("_Model$", "", x)))
}

model_prediction_objects <- lapply(individual_predictions, function(mat) {
  list(Times = common_times, CIFs = mat)
})

names(model_prediction_objects) <- vapply(names(individual_predictions), pretty_cr_model_name, character(1))

all_predictions <- c(
  model_prediction_objects,
  list(Ensemble = list(Times = common_times, CIFs = ensemble_cif))
)

names(all_predictions)


## ----evaluate-ensemble--------------------------------------------------------
actual_times <- test_data[[timevar]]
actual_events <- test_data[[eventvar]]
surv_obj <- Surv(actual_times, actual_events, type = "mstate")

eval_time <- median(actual_times[actual_events == primary_event_code], na.rm = TRUE)

etl_results <- list()

performance_summary <- do.call(rbind, lapply(names(all_predictions), function(model_name) {
  pred_obj <- all_predictions[[model_name]]

  concordance <- tryCatch({
    timedepConcordanceCR(
      SurvObj = surv_obj,
      Predictions = pred_obj$CIFs,
      time = eval_time,
      cause = primary_event_code,
      pred_times = pred_obj$Times,
      TestMat = test_data[, expvars]
    )
  }, error = function(e) NA_real_)

  # BrierScoreCR now uses 'eval_times' parameter (standardized across survival and CR)
  brier <- tryCatch({
    # CIF matrices retain the time x observation orientation for all metrics
    BrierScoreCR(
      SurvObj = surv_obj,
      Predictions = pred_obj$CIFs,
      eval_times = eval_time,  # New standardized parameter name
      cause = primary_event_code,
      TestMat = test_data[, expvars],
      pred_times = pred_obj$Times
    )
  }, error = function(e) NA_real_)

  etl_values <- tryCatch({
    etl_list <- CalculateExpectedTimeLost(
      PredictedCurves = list(list(NewProbs = pred_obj$CIFs)),
      modeltypes = c("CR"),
      times = pred_obj$Times,
      UL = max(actual_times)
    )
    etl_list[[1]]
  }, error = function(e) rep(NA_real_, ncol(pred_obj$CIFs)))

  etl_results[[model_name]] <<- etl_values

  data.frame(
    Model = model_name,
    Concordance = concordance,
    Brier_Score = brier,
    Mean_ETL = mean(etl_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

performance_summary


## ----expected-time-lost-------------------------------------------------------
max_time <- max(actual_times)
ensemble_etl <- etl_results[["Ensemble"]]
if (is.null(ensemble_etl)) {
  ensemble_etl <- rep(NA_real_, ncol(ensemble_cif))
}

etl_summary <- data.frame(
  statistic = c(
    "Maximum follow-up (months)",
    "Mean ETL (months)",
    "ETL range (months)"
  ),
  value = c(
    round(max_time, 2),
    round(mean(ensemble_etl, na.rm = TRUE), 2),
    if (all(is.na(ensemble_etl))) NA_character_ else paste(
      round(range(ensemble_etl, na.rm = TRUE), 2),
      collapse = " – "
    )
  )
)

etl_summary


## ----visualization, fig.width=12, fig.height=10-------------------------------
print(performance_summary)

ordered_models <- performance_summary[order(-performance_summary$Concordance), ]
ordered_models <- ordered_models[!is.na(ordered_models$Concordance), ]

# Plot all models with ensemble highlighted in black
plot_predictions <- all_predictions

p_cif <- plot_cif_curves(
  predictions = plot_predictions,
  patients_to_plot = seq_len(min(3, nrow(test_data))),
  highlight_ensemble = TRUE,
  title = "Predicted Cumulative Incidence Functions (CIF) - All Models (Ensemble in Black)"
)
print(p_cif)

ensemble_etl_clean <- ensemble_etl[!is.na(ensemble_etl)]
if (length(ensemble_etl_clean) > 0) {
  etl_df <- data.frame(
    Subject = seq_along(ensemble_etl_clean),
    ETL = ensemble_etl_clean
  )

  p_etl <- ggplot(etl_df, aes(x = ETL)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black") +
    labs(
      title = "Ensemble Expected Time Lost Distribution",
      x = "Expected Time Lost (months)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  print(p_etl)
}

test_cif_phase <- cuminc(
  ftime = test_data$ftime,
  fstatus = test_data$status,
  group = test_data$phase
)

cif_phase_plot_data <- data.frame()
for (cif_name in names(test_cif_phase)) {
  if (cif_name %in% c("Tests")) next

  cif_obj <- test_cif_phase[[cif_name]]
  parts <- strsplit(cif_name, " ")[[1]]

  if (length(parts) >= 2) {
    event_type <- as.numeric(parts[length(parts)])
    phase_name <- paste(parts[1:(length(parts) - 1)], collapse = " ")

    cif_data <- data.frame(
      Time = cif_obj$time,
      CIF = cif_obj$est,
      Phase = phase_name,
      Event = ifelse(event_type == 1, "Relapse", "TRM")
    )
    cif_phase_plot_data <- rbind(cif_phase_plot_data, cif_data)
  }
}

if (nrow(cif_phase_plot_data) > 0) {
  p_phase_relapse <- ggplot(
    cif_phase_plot_data[cif_phase_plot_data$Event == "Relapse", ],
    aes(x = Time, y = CIF, color = Phase)
  ) +
    geom_line(linewidth = 1.2) +
    labs(
      title = "Cumulative Incidence of Relapse by Disease Phase",
      subtitle = "Test Set Data",
      x = "Time (months)",
      y = "Cumulative Incidence of Relapse"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  p_phase_trm <- ggplot(
    cif_phase_plot_data[cif_phase_plot_data$Event == "TRM", ],
    aes(x = Time, y = CIF, color = Phase)
  ) +
    geom_line(linewidth = 1.2) +
    labs(
      title = "Cumulative Incidence of TRM by Disease Phase",
      subtitle = "Test Set Data",
      x = "Time (months)",
      y = "Cumulative Incidence of TRM"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  print(p_phase_relapse)
  print(p_phase_trm)
}


## ----save-models--------------------------------------------------------------
ensemble_file <- tempfile("cr_ensemble_model_", fileext = ".rds")
SaveEnsemble(models, ensemble_file)

loaded_ensemble <- LoadEnsemble(ensemble_file)

reloaded_predictions <- PredictCRModels(
  models = loaded_ensemble,
  newdata = test_data[seq_len(min(5, nrow(test_data))), ],
  new_times = common_times,
  ensemble_method = "average"
)

data.frame(
  object = c("Loaded ensemble class", "Predictions generated for subjects"),
  value = c(
    paste(class(loaded_ensemble), collapse = ", "),
    length(reloaded_predictions$NewProbs)
  )
)
