devtools::load_all()
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)
if (!require("cmprsk", quietly = TRUE)) install.packages("cmprsk")
library(cmprsk)

# Load data
data_path <- system.file("extdata", "bmtcrr_competing_risks.csv", package = "ml4time2event")
if (data_path == "") {
  data_path <- "inst/extdata/bmtcrr_competing_risks.csv"
}
bmt_data <- read.csv(data_path)

# Preprocessing
bmt_data$Sex <- factor(bmt_data$Sex)
bmt_data$D <- factor(bmt_data$D)
bmt_data$Phase <- factor(bmt_data$Phase)
bmt_data$Source <- factor(bmt_data$Source)

# Split data
set.seed(123)
train_indices <- sample(nrow(bmt_data), floor(0.7 * nrow(bmt_data)))
train_data <- bmt_data[train_indices, ]
test_data <- bmt_data[-train_indices, ]

# Train a simple model
timevar <- "ftime"
eventvar <- "Status"
expvars <- c("Sex", "D", "Phase", "Age", "Source")

cox_model <- CRModel_Cox(
  data = train_data,
  expvars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  failcode = 1
)

# Make predictions
pred <- Predict_CRModel_Cox(cox_model, test_data)
print("Prediction structure:")
print(str(pred))
print("Times range:")
print(range(pred$Times))
print("CIF dimensions:")
print(dim(pred$CIF))
print("CIF summary:")
print(summary(as.vector(pred$CIF)))

# Test concordance calculation
actual_times <- test_data[[timevar]]
actual_events <- test_data[[eventvar]]
surv_obj <- Surv(actual_times, actual_events, type = "mstate")

eval_time <- median(actual_times[actual_events == 1])
print(paste("Eval time:", eval_time))

time_idx <- which.min(abs(pred$Times - eval_time))
cif_at_eval <- pred$CIF[time_idx, ]

print("CIF at eval time summary:")
print(summary(cif_at_eval))
print("Any NA in CIF?")
print(any(is.na(cif_at_eval)))

print("Surv object summary:")
print(summary(surv_obj))

# Try concordance
tryCatch({
  concordance_val <- timedepConcordanceCR(
    SurvObj = surv_obj,
    Predictions = matrix(cif_at_eval, ncol = 1),
    time = eval_time,
    cause = 1,
    TestMat = test_data[, expvars]
  )
  print(paste("Concordance:", concordance_val))
}, error = function(e) {
  print(paste("Error:", e$message))
  print("Debugging...")
  print("obstimes summary:")
  print(summary(surv_obj[, "time"]))
  print("obsevents summary:")
  print(summary(surv_obj[, "status"]))
})