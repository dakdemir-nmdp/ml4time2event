# Tests for plotting functions
# TDD approach - write tests first, then ensure functions work

test_that("plot_survival_curves validates inputs correctly", {
  # Test NULL predictions
  expect_error(plot_survival_curves(NULL), "'predictions' cannot be NULL")
  
  # Test invalid prediction objects
  invalid_pred <- list(Probs = NULL, Times = NULL)
  expect_error(plot_survival_curves(invalid_pred), 
               "Each prediction object must have 'Probs' and 'Times' components")
  
  # Test non-matrix Probs
  invalid_pred2 <- list(Probs = c(1, 2, 3), Times = c(1, 2, 3))
  expect_error(plot_survival_curves(invalid_pred2), "'Probs' component must be a matrix")
  
  # Test non-numeric Times
  invalid_pred3 <- list(Probs = matrix(c(1, 0.8, 0.6), nrow = 3, ncol = 1), 
                        Times = c("a", "b", "c"))
  expect_error(plot_survival_curves(invalid_pred3), "'Times' component must be numeric")
})

test_that("plot_survival_curves handles single prediction object", {
  # Create valid single prediction object
  pred <- list(
    Probs = matrix(c(1, 0.8, 0.6, 1, 0.9, 0.7), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  # Should not error
  expect_no_error(plot_survival_curves(pred))
  
  # Should return ggplot object
  p <- plot_survival_curves(pred)
  expect_s3_class(p, "ggplot")
})

test_that("plot_survival_curves handles multiple prediction objects", {
  # Create valid prediction objects
  pred1 <- list(
    Probs = matrix(c(1, 0.8, 0.6, 1, 0.9, 0.7), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  pred2 <- list(
    Probs = matrix(c(1, 0.85, 0.65, 1, 0.95, 0.75), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  predictions <- list(Cox = pred1, RF = pred2)
  
  # Should not error
  expect_no_error(plot_survival_curves(predictions))
  
  # Should return ggplot object
  p <- plot_survival_curves(predictions)
  expect_s3_class(p, "ggplot")
})

test_that("plot_survival_curves validates patients_to_plot parameter", {
  pred <- list(
    Probs = matrix(c(1, 0.8, 0.6, 1, 0.9, 0.7), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  # Test patients_to_plot larger than available patients
  expect_error(plot_survival_curves(pred, patients_to_plot = 1:5),
               "'patients_to_plot' contains indices larger than number of patients")
  
  # Test negative patient indices  
  expect_error(plot_survival_curves(pred, patients_to_plot = c(-1, 0)),
               "'patients_to_plot' must contain positive integers")
  
  # Test valid patient indices
  expect_no_error(plot_survival_curves(pred, patients_to_plot = 1:2))
})

test_that("plot_cif_curves validates inputs correctly", {
  # Test NULL predictions
  expect_error(plot_cif_curves(NULL), "'predictions' cannot be NULL")
  
  # Test invalid prediction objects (missing CIFs)
  invalid_pred <- list(CIFs = NULL, Times = NULL)
  expect_error(plot_cif_curves(invalid_pred), 
               "Each prediction object must have 'CIFs' and 'Times' components")
  
  # Test non-matrix CIFs
  invalid_pred2 <- list(CIFs = c(0, 0.1, 0.2), Times = c(0, 6, 12))
  expect_error(plot_cif_curves(invalid_pred2), "'CIFs' component must be a matrix")
})

test_that("plot_cif_curves handles single prediction object", {
  # Create valid single prediction object for CIF
  pred <- list(
    CIFs = matrix(c(0, 0.1, 0.2, 0, 0.05, 0.15), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  # Should not error
  expect_no_error(plot_cif_curves(pred))
  
  # Should return ggplot object
  p <- plot_cif_curves(pred)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cif_curves handles multiple prediction objects", {
  # Create valid prediction objects
  pred1 <- list(
    CIFs = matrix(c(0, 0.1, 0.2, 0, 0.05, 0.15), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  pred2 <- list(
    CIFs = matrix(c(0, 0.12, 0.22, 0, 0.07, 0.17), nrow = 3, ncol = 2),
    Times = c(0, 6, 12)
  )
  
  predictions <- list(Cox = pred1, FineGray = pred2)
  
  # Should not error
  expect_no_error(plot_cif_curves(predictions))
  
  # Should return ggplot object
  p <- plot_cif_curves(predictions)
  expect_s3_class(p, "ggplot")
})

test_that("plot_model_performance validates inputs correctly", {
  # Test non-data.frame input
  expect_error(plot_model_performance(list()), "'performance_df' must be a data frame")
  
  # Test missing Model column
  df_no_model <- data.frame(Concordance = c(0.7, 0.8))
  expect_error(plot_model_performance(df_no_model), "'performance_df' must have a 'Model' column")
  
  # Test missing metric column
  df_no_metric <- data.frame(Model = c("Cox", "RF"))
  expect_error(plot_model_performance(df_no_metric, metric = "concordance"),
               "'performance_df' must have a Concordance column")
})

test_that("plot_model_performance creates valid plots", {
  # Create valid performance data
  performance_df <- data.frame(
    Model = c("Cox", "RF", "Ensemble"),
    Concordance = c(0.72, 0.75, 0.78),
    Brier_Score = c(0.15, 0.12, 0.10)
  )
  
  # Test concordance plot
  p1 <- plot_model_performance(performance_df, metric = "concordance")
  expect_s3_class(p1, "ggplot")
  
  # Test brier score plot
  p2 <- plot_model_performance(performance_df, metric = "brier")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_model_performance handles missing values", {
  # Create performance data with NAs
  performance_df <- data.frame(
    Model = c("Cox", "RF", "Ensemble"),
    Concordance = c(0.72, NA, 0.78)
  )
  
  # Should remove NA rows and still work
  p <- plot_model_performance(performance_df, metric = "concordance")
  expect_s3_class(p, "ggplot")
  
  # Test all NA values
  performance_df_all_na <- data.frame(
    Model = c("Cox", "RF"),
    Concordance = c(NA, NA)
  )
  
  expect_error(plot_model_performance(performance_df_all_na, metric = "concordance"),
               "No valid values found for Concordance")
})

test_that("plotting functions handle ensemble highlighting", {
  # Test survival curves with ensemble
  pred1 <- list(
    Probs = matrix(c(1, 0.8, 0.6), nrow = 3, ncol = 1),
    Times = c(0, 6, 12)
  )
  
  pred2 <- list(
    Probs = matrix(c(1, 0.85, 0.65), nrow = 3, ncol = 1),
    Times = c(0, 6, 12)
  )
  
  predictions <- list(Cox = pred1, Ensemble = pred2)
  
  # Should handle ensemble highlighting
  p1 <- plot_survival_curves(predictions, highlight_ensemble = TRUE)
  expect_s3_class(p1, "ggplot")
  
  # Test CIF curves with ensemble
  cif_pred1 <- list(
    CIFs = matrix(c(0, 0.1, 0.2), nrow = 3, ncol = 1),
    Times = c(0, 6, 12)
  )
  
  cif_pred2 <- list(
    CIFs = matrix(c(0, 0.12, 0.22), nrow = 3, ncol = 1),
    Times = c(0, 6, 12)
  )
  
  cif_predictions <- list(Cox = cif_pred1, Ensemble = cif_pred2)
  
  p2 <- plot_cif_curves(cif_predictions, highlight_ensemble = TRUE)
  expect_s3_class(p2, "ggplot")
  
  # Test performance plot with ensemble
  performance_df <- data.frame(
    Model = c("Cox", "RF", "Ensemble"),
    Concordance = c(0.72, 0.75, 0.78)
  )
  
  p3 <- plot_model_performance(performance_df, metric = "concordance", highlight_ensemble = TRUE)
  expect_s3_class(p3, "ggplot")
})