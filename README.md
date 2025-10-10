# ml4time2event

Machine learning for time to event analysis. This package provides a comprehensive collection of tools for predicting time-to-event outcomes using various statistical and machine learning methods, including support for both standard survival analysis and competing risks scenarios.

## Features

- **Multiple Model Types**: Cox regression, Random Forests, XGBoost, GAM, BART, DeepSurv, and more
- **Competing Risks Support**: Specialized models for competing risks analysis (Fine-Gray, cause-specific hazards)
- **Ensemble Methods**: Simple averaging, weighted averaging, and super learner stacking
- **Comprehensive Metrics**: C-index, Brier score, integrated metrics, and expected time lost
- **Robust Predictions**: Handles missing data, factor levels, and edge cases gracefully
- **Production Ready**: Extensive test coverage (1917 passing tests) and validated vignettes

## Installation

### Install from GitHub

To install the package from GitHub, run the following commands in R:

```r
# Install devtools if not already installed
# install.packages("devtools")

library(devtools)
install_github("dakdemir-nmdp/ml4time2event")
```

This command will install ml4time2event along with all required dependencies.

### System Requirements

- R >= 3.5.0
- Recommended: R >= 4.0.0 for optimal performance

### Dependencies

The package depends on several R packages that will be automatically installed:

**Required dependencies:**
- `magrittr` - Pipe operator for cleaner code
- `party` - Recursive partitioning methods
- `prodlim` - Survival analysis functions
- `survival` - Core survival analysis (Surv objects, Cox models)
- `labelled` - Data label management
- `xgboost` - Gradient boosting methods

**Suggested packages** (optional, for extended functionality):
- `randomForestSRC` - Random forests for survival
- `glmnet` - Penalized regression (Lasso, Ridge, Elastic Net)
- `mgcv` - Generalized additive models
- `gbm` - Gradient boosting machines
- `pec` - Prediction error curves
- `fastcmprsk` - Fast competing risks models
- `BART` - Bayesian additive regression trees
- `data.table`, `dplyr` - Data manipulation
- `ggplot2` - Visualization (for vignettes)

## Quick Start

Load the package in your R session:

```r
library(ml4time2event)
```

### Basic Survival Analysis Example

```r
# Load example data
library(survival)
data(veteran)

# Fit multiple survival models
models <- RunSurvModels(
  datatrain = veteran,
  ExpVars = c("trt", "karno", "diagtime", "age", "prior"),
  timevar = "time",
  eventvar = "status",
  models = c("coxph", "RF", "xgboost")
)

# Make predictions
predictions <- PredictSurvModels(
  models = models,
  newdata = veteran[1:10, ],
  newtimes = c(30, 90, 180),
  ensemble_method = "average"
)

# Access survival probabilities (rows = times, cols = observations)
head(predictions$NewProbs)
```

### Competing Risks Example

```r
# Fit competing risks models
cr_models <- RunCRModels(
  datatrain = your_data,
  ExpVars = c("age", "sex", "treatment"),
  timevar = "time",
  eventvar = "event",  # 0 = censored, 1 = event1, 2 = event2
  models = c("FineGray", "Cox", "GAM")
)

# Predict cumulative incidence functions
cr_predictions <- PredictCRModels(
  models = cr_models,
  newdata = test_data,
  newtimes = c(365, 730, 1095)
)
```

## Testing

The package includes comprehensive test coverage with 1917 passing tests.

### Run All Tests

```r
# Install test dependencies
install.packages("testthat")

# Run all tests
devtools::test()
```

Expected output:
```
FAIL 0 | WARN 44 | SKIP 1 | PASS 1917
```

### Run Specific Test Files

```r
# Test survival models
testthat::test_file("tests/testthat/test_surv_cox.R")

# Test competing risks models
testthat::test_file("tests/testthat/test_cr_xgboost.R")

# Test ensemble methods
testthat::test_file("tests/testthat/test_super_learner.R")

# Test metrics
testthat::test_file("tests/testthat/test_surv_metrics.R")
```

## Comprehensive Examples

The package includes two comprehensive vignettes demonstrating full analysis workflows:

### 1. Survival Analysis Vignette

Demonstrates the complete survival analysis workflow using the follicular lymphoma dataset:

```r
# Run the survival analysis vignette
source("vignettes/comprehensive_survival_analysis.R")
```

**What it covers:**
- Data loading and preprocessing
- Multiple survival models (Cox, Random Forest, XGBoost, GLMNet, GBM, GAM, SurvReg, BART, RuleFit, DeepSurv)
- Model evaluation and comparison (C-index, Brier score)
- Ensemble methods (averaging, weighted, super learner)
- Expected time lost calculations
- Comprehensive visualization
- Model persistence (saving/loading models)

**Output:**
- Model files: `.rds` files for each trained model
- Visualization: `follic_survival_analysis.pdf` with survival curves, C-index comparison, and more

### 2. Competing Risks Analysis Vignette

Demonstrates competing risks analysis using bone marrow transplant data:

```r
# Run the competing risks vignette
source("vignettes/comprehensive_competing_risks_analysis.R")
```

**What it covers:**
- Competing risks data preparation
- Multiple CR models (Cox cause-specific, Fine-Gray, RF, XGBoost, GAM, BART, DeepSurv, RuleFit, SurvReg)
- Cumulative Incidence Function (CIF) calculations
- Expected Time Lost for competing risks
- CR-specific metrics evaluation
- Ensemble methods for competing risks
- Comprehensive CIF and ETL visualization
- Model persistence

**Output:**
- Model files: 10 `.rds` files for trained models
- Visualization: `bmt_competing_risks_analysis.pdf` (4 pages) with CIF curves, model comparison, and ETL analysis

## Key Functions

### Survival Analysis

- `RunSurvModels()` - Train multiple survival models
- `PredictSurvModels()` - Generate predictions with ensemble options
- `SurvModel_Cox()` - Fit Cox proportional hazards model
- `SurvModel_RF()` - Fit random forest survival model
- `SurvModel_xgboost()` - Fit XGBoost survival model
- `timedepConcordance()` - Calculate time-dependent concordance
- `BrierScore()` - Calculate Brier score for survival predictions

### Competing Risks

- `RunCRModels()` - Train multiple competing risks models
- `PredictCRModels()` - Generate CIF predictions with ensemble options
- `CRModel_FineGray()` - Fit Fine-Gray subdistribution hazard model
- `CRModel_Cox()` - Fit cause-specific Cox models
- `CRModel_GAM()` - Fit generalized additive model for competing risks
- `CRModel_xgboost()` - Fit XGBoost for competing risks

### Ensemble Methods

- Simple averaging: `ensemble_method = "average"`
- Weighted averaging: `ensemble_method = "weighted"` (with `model_weights`)
- Super learner: `ensemble_method = "super_learner"` (with training data)
- `ComputeSuperLearnerWeights()` - Pre-compute optimal ensemble weights

### Utilities

- `survprobMatInterpolator()` - Interpolate survival probabilities to new time points
- `cifMatInterpolator()` - Interpolate cumulative incidence functions
- `buildObservedSurvivalMatrix()` - Create target matrix for super learner training
- `buildObservedCIFMatrix()` - Create CIF target matrix for super learner

## Contributing

Contributions, suggestions, and bug reports are welcome. Please use the GitHub repository to submit issues or pull requests.

### Development Setup

```r
# Clone the repository
# git clone https://github.com/dakdemir-nmdp/ml4time2event.git

# Install development dependencies
devtools::install_dev_deps()

# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Build documentation
devtools::document()

# Check package
devtools::check()
```

## Citation

If you use this package in your research, please cite:

```
Akdemir, D. (2025). ml4time2event: Machine Learning for Time-to-Event Analysis.
R package version 0.1.0. https://github.com/dakdemir-nmdp/ml4time2event
```

## License

This project is licensed under the GPL (>= 2) License.

## Contact

**Author and Maintainer:** Deniz Akdemir
**Email:** dakdemir@nmdp.org

## Acknowledgements

We thank all contributors and users for their support and feedback.

---

**Package Status:** Production Ready
**Test Coverage:** 1917 passing tests
**Last Updated:** 2025

Enjoy using ml4time2event!
