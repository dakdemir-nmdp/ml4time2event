library(survival)

test_that("aalenJohansenFromCoxModels errors without fallback", {
  # 1. Setup: Train a Cox model
  set.seed(123)
  train_data <- data.frame(
    time = rexp(50),
    status = sample(0:2, 50, replace = TRUE),
    x = factor(sample(c("A", "B"), 50, replace = TRUE))
  )
  
  # Create two cause-specific models
  model1 <- coxph(Surv(time, status == 1) ~ x, data = train_data)
  model2 <- coxph(Surv(time, status == 2) ~ x, data = train_data)
  
  cox_models <- list("1" = model1, "2" = model2)
  
  # 2. Create newdata with a factor level not in the training data
  newdata <- data.frame(x = factor("C", levels = c("A", "B", "C")))
  times <- c(0.1, 0.5, 1.0)
  
  # 3. Test the current behavior (should not error, but returns a marginal prediction)
  # The `tryCatch` will hide the error from the new factor level.
  # We can't easily test the *value* of the fallback, but we can show it doesn't error out.
  # After the fix, this call should fail.
  
  # After the fix, this call should fail with a specific error message.
  expect_error(
    aalenJohansenFromCoxModels(cox_models, newdata, times, event_of_interest = "1"),
    "factor x has new level C" # survfit throws an error with this message
  )
})
