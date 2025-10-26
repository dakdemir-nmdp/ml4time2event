library(testthat)
library(here)
# library(data.table) # Removed
library(survival)
library(pracma) # For mocking Integrator

# Assuming the functions are available in the environment
# Use here::here for robustness
source(here::here("R/outcome_prediction.R"))
# Source dependencies needed by the functions under test
source(here::here("R/surv_ensemble.R"))
source(here::here("R/cr_ensemble.R"))
source(here::here("R/math_utils.R"))
# Need mocks or actual functions for PredictSurvModels, PredictCRModels, Integrator
# The mocks below will overwrite the sourced functions for testing purposes

context("Testing outcome_prediction functions")

# --- Test Data Setup ---
# Survival data
set.seed(789)
n_obs_surv <- 10
surv_data_pred <- data.frame(
  time = rexp(n_obs_surv, rate = 0.05),
  status = sample(0:1, n_obs_surv, replace = TRUE, prob = c(0.3, 0.7)),
  x1 = rnorm(n_obs_surv),
  x2 = factor(sample(c("C", "D"), n_obs_surv, replace = TRUE)),
  stringsAsFactors = FALSE
)

# Competing risks data
set.seed(123)
n_obs_cr <- 10
cr_data_pred <- data.frame(
  time = pmin(rexp(n_obs_cr, rate = 0.1), rexp(n_obs_cr, rate = 0.15)),
  status = sample(0:2, n_obs_cr, replace = TRUE, prob = c(0.2, 0.4, 0.4)),
  x1 = rnorm(n_obs_cr),
  x2 = factor(sample(c("A", "B"), n_obs_cr, replace = TRUE)),
  stringsAsFactors = FALSE
)

# Time points
times_pred <- c(5, 10, 15, 20)

# --- Mock Models and Predictions ---
# Mock model objects (simple lists)
mock_surv_model <- list(model_type = "SURV_MOCK", id = 1)
mock_cr_model <- list(model_type = "CR_MOCK", id = 2)
models_list_pred <- list(surv1 = mock_surv_model, cr1 = mock_cr_model, surv2 = mock_surv_model)
model_types_pred <- c("SURV", "CR", "SURV")


# --- Tests for PredictAllPossibleOutcomesSurvOrCifs ---

test_that("PredictAllPossibleOutcomesSurvOrCifs calls correct predict functions (mocked)", {
  skip_if_not_installed("mockery")

  # Mock the underlying prediction functions for this test
  mockery::stub(PredictAllPossibleOutcomesSurvOrCifs, 'PredictSurvModels', function(models, newdata, new_times, ...) {
      list(
        MockSurvModel = matrix(0.5, nrow=length(new_times), ncol=nrow(newdata)),
        NewProbs = matrix(0.5, nrow=length(new_times), ncol=nrow(newdata))
      )
  })
   mockery::stub(PredictAllPossibleOutcomesSurvOrCifs, 'PredictCRModels', function(models, newdata, new_times, ...) {
      list(
        MockCRModel = list(matrix(0.2, nrow=length(new_times), ncol=nrow(newdata))), # Only CIF1 needed for NewProbs
        NewProbs = matrix(0.2, nrow=length(new_times), ncol=nrow(newdata))
      )
  })

  predictions <- PredictAllPossibleOutcomesSurvOrCifs(data = surv_data_pred, # Use actual test data
                                                      modelslist = models_list_pred,
                                                      modeltypes = model_types_pred,
                                                      times = times_pred)

  expect_type(predictions, "list")
  expect_length(predictions, length(models_list_pred))

  # Check structure of SURV prediction output (based on mock)
  expect_type(predictions[[1]], "list")
  expect_named(predictions[[1]], c("MockSurvModel", "NewProbs"))
  expect_true(is.matrix(predictions[[1]]$NewProbs))
  expect_equal(dim(predictions[[1]]$NewProbs), c(length(times_pred), nrow(surv_data_pred))) # Use actual data dim
  expect_true(all(predictions[[1]]$NewProbs == 0.5)) # Check mock value

  # Check structure of CR prediction output (based on mock)
   expect_type(predictions[[2]], "list")
   expect_named(predictions[[2]], c("MockCRModel", "NewProbs"))
   expect_true(is.matrix(predictions[[2]]$NewProbs)) # Should contain CIF1
   expect_equal(dim(predictions[[2]]$NewProbs), c(length(times_pred), nrow(surv_data_pred))) # Use actual data dim
   expect_true(all(predictions[[2]]$NewProbs == 0.2)) # Check mock value

  # Check third model (SURV)
  expect_type(predictions[[3]], "list")
  expect_named(predictions[[3]], c("MockSurvModel", "NewProbs"))
  expect_true(is.matrix(predictions[[3]]$NewProbs))
  expect_true(all(predictions[[3]]$NewProbs == 0.5)) # Check mock value

  # Stubs created within test_that are typically reset automatically
})

test_that("PredictAllPossibleOutcomesSurvOrCifs handles NULL models and unsupported types", {
  skip_if_not_installed("mockery")
  models_list_bad <- list(surv1 = mock_surv_model, cr1 = NULL, other = list(), surv2 = mock_surv_model)
  model_types_bad <- c("SURV", "CR", "OTHER", "SURV")

  # Mock the underlying prediction functions for this test
  mockery::stub(PredictAllPossibleOutcomesSurvOrCifs, 'PredictSurvModels', function(...) list(NewProbs=matrix(0.5)))
  mockery::stub(PredictAllPossibleOutcomesSurvOrCifs, 'PredictCRModels', function(...) list(NewProbs=matrix(0.2)))

  # Combine warning checks into one expect_warning with a regex
  # The function call will generate both warnings, testthat captures them.
  # We check if *any* of the captured warnings match the regex.
  expect_warning(
      predictions <- PredictAllPossibleOutcomesSurvOrCifs(
          data = surv_data_pred,
          modelslist = models_list_bad,
          modeltypes = model_types_bad,
          times = times_pred
      ),
      regexp = "(Model object at index 2 is NULL)|(Unsupported model type 'OTHER' at index 3)"
  )

  # The predictions variable should now be assigned from the single call within expect_warning
  expect_length(predictions, 4)
  expect_true(!is.null(predictions[[1]]) && !is.na(predictions[[1]])) # Valid SURV (returns list)
  expect_true(is.na(predictions[[2]]))  # NULL model
  expect_true(is.na(predictions[[3]]))  # Unsupported type
  expect_true(!is.na(predictions[[4]])) # Valid SURV (returns list)

  # Stubs created within test_that are typically reset automatically
})

test_that("PredictAllPossibleOutcomesSurvOrCifs requires consistent lengths", {
   expect_error(PredictAllPossibleOutcomesSurvOrCifs(data = surv_data_pred, # Use actual test data
                                                      modelslist = models_list_pred,
                                                      modeltypes = model_types_pred[-1], # Mismatched length
                                                      times = times_pred),
                "Length of 'modelslist' and 'modeltypes' must be equal.")
})


# --- Tests for CalculateExpectedTimeLost ---

test_that("CalculateExpectedTimeLost calculates RMTL for SURV and CR", {
  # Generate mock prediction outputs (list with NewProbs matrix: times x obs)
  n_test_local <- nrow(surv_data_pred) # Define n_test locally
  n_times <- length(times_pred)
  mock_pred_surv <- list(
    NewProbs = matrix(
      seq(0.9, 0.5, length.out = n_test_local * n_times),
      nrow = n_times,
      ncol = n_test_local
    )
  )
  mock_pred_cr <- list(
    NewProbs = matrix(
      seq(0.1, 0.4, length.out = n_test_local * n_times),
      nrow = n_times,
      ncol = n_test_local
    )
  )
  predictions_list_calc <- list(mock_pred_surv, mock_pred_cr, mock_pred_surv)


  UL_time <- 15
  rmtl_list <- CalculateExpectedTimeLost(PredictedCurves = predictions_list_calc,
                                         modeltypes = model_types_pred,
                                         times = times_pred,
                                         UL = UL_time, LL = 0)

  expect_type(rmtl_list, "list")
  expect_length(rmtl_list, length(predictions_list_calc))

  # Check SURV RMTL (index 1)
  expect_type(rmtl_list[[1]], "double")
  expect_length(rmtl_list[[1]], n_test_local) # Use local n_test
  # Check calculation based on actual Integrator: integrate(1-ProbMatrix)
  event_probs1 <- 1 - predictions_list_calc[[1]]$NewProbs
  expected_rmtl1 <- apply(event_probs1, 2, function(scores) Integrator(times_pred, scores, c(0, UL_time), FALSE))
  expect_equal(rmtl_list[[1]], expected_rmtl1)

  # Check CR RMTL (index 2)
  expect_type(rmtl_list[[2]], "double")
  expect_length(rmtl_list[[2]], n_test_local) # Use local n_test
  # Check calculation based on actual Integrator: integrate(ProbMatrix)
  event_probs2 <- predictions_list_calc[[2]]$NewProbs
  expected_rmtl2 <- apply(event_probs2, 2, function(scores) Integrator(times_pred, scores, c(0, UL_time), FALSE))
  expect_equal(rmtl_list[[2]], expected_rmtl2)

   # Check SURV RMTL (index 3)
  expect_type(rmtl_list[[3]], "double")
  expect_length(rmtl_list[[3]], n_test_local) # Use local n_test
  event_probs3 <- 1 - predictions_list_calc[[3]]$NewProbs
  expected_rmtl3 <- apply(event_probs3, 2, function(scores) Integrator(times_pred, scores, c(0, UL_time), FALSE))
  expect_equal(rmtl_list[[3]], expected_rmtl3)
})

test_that("CalculateExpectedTimeLost handles invalid inputs", {
   # Generate valid mock predictions
   n_test_local <- nrow(surv_data_pred) # Define n_test locally
   n_times <- length(times_pred)
   mock_pred_surv <- list(
     NewProbs = matrix(
       seq(0.9, 0.5, length.out = n_test_local * n_times),
       nrow = n_times,
       ncol = n_test_local
     )
   )
   mock_pred_cr <- list(
     NewProbs = matrix(
       seq(0.1, 0.4, length.out = n_test_local * n_times),
       nrow = n_times,
       ncol = n_test_local
     )
   )
   predictions_list_calc <- list(mock_pred_surv, mock_pred_cr, mock_pred_surv)

   # Mismatched lengths
   expect_error(CalculateExpectedTimeLost(PredictedCurves = predictions_list_calc,
                                          modeltypes = model_types_pred[-1], # Mismatch
                                          times = times_pred, UL = 15),
                "Length of 'PredictedCurves' and 'modeltypes' must be equal.")

   # Invalid UL/LL
   expect_error(CalculateExpectedTimeLost(PredictedCurves = predictions_list_calc, modeltypes = model_types_pred, times = times_pred, UL = 5, LL = 10),
                "'UL' must be a single numeric value greater than 'LL'.")
   expect_error(CalculateExpectedTimeLost(PredictedCurves = predictions_list_calc, modeltypes = model_types_pred, times = times_pred, UL = c(10, 15)),
                "'UL' must be a single numeric value")
   expect_error(CalculateExpectedTimeLost(PredictedCurves = predictions_list_calc, modeltypes = model_types_pred, times = times_pred, UL = 15, LL = -1),
                 "'LL' must be a single non-negative numeric value")

   # Missing NewProbs or invalid prediction object
   predictions_list_bad <- predictions_list_calc
   predictions_list_bad[[1]] <- list() # Missing NewProbs
   predictions_list_bad[[2]] <- NA     # Invalid object
   expect_warning(rmtl_list_bad <- CalculateExpectedTimeLost(PredictedCurves = predictions_list_bad, modeltypes = model_types_pred, times = times_pred, UL = 15),
                  "Invalid prediction output or missing 'NewProbs' at index 1")
   expect_warning(rmtl_list_bad <- CalculateExpectedTimeLost(PredictedCurves = predictions_list_bad, modeltypes = model_types_pred, times = times_pred, UL = 15),
                  "Invalid prediction output or missing 'NewProbs' at index 2")
   expect_true(is.na(rmtl_list_bad[[1]]))
   expect_true(is.na(rmtl_list_bad[[2]]))
   expect_false(is.na(rmtl_list_bad[[3]][1])) # Third one should still calculate
})
