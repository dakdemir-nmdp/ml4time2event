set.seed(123)
library(ml4time2event)

# Load data
data <- get_bmt_competing_risks_data()
timevar <- "ftime"
eventvar <- "status"
expvars <- c("sex", "d", "phase", "age", "source")

# Split
split <- t2edata_split(data, prop = 0.7)
train_data <- split$Train
test_data <- split$Test

# Train models
cat("\n=== TRAINING MODELS ===\n")
models <- RunCRModels(
  datatrain = train_data,
  ExpVars = expvars,
  timevar = timevar,
  eventvar = eventvar,
  models = c("cox", "FG")  # Just test with 2 models for speed
)

cat("\n=== MODEL NAMES ===\n")
print(names(models))
cat("\nTotal models:", length(names(models)), "\n")

# Predict
cat("\n=== PREDICTING ===\n")
predictions <- PredictCRModels(
  models = models,
  newdata = test_data,
  new_times = seq(0, max(test_data[[timevar]]), length.out = 50),
  ensemble_method = "average"
)

cat("\n=== PREDICTION COMPONENTS ===\n")
cat("ModelPredictions names:", names(predictions$ModelPredictions), "\n")
cat("NewProbs is NULL:", is.null(predictions$NewProbs), "\n")
if (!is.null(predictions$NewProbs)) {
  cat("NewProbs dimensions:", dim(predictions$NewProbs), "\n")
}

# Evaluate
cat("\n=== EVALUATING ===\n")
leaderboard <- EvaluateCRModels(
  models = models,
  data = test_data,
  timevar = timevar,
  eventvar = eventvar,
  ensemble_method = "average",
  cause = 1L
)

cat("\n=== LEADERBOARD ===\n")
print(leaderboard[, c("model", "integrated_c", "integrated_brier", "rank")])
