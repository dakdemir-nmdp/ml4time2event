# ml4time2event

Machine learning for time-to-event analysis. Provides tools for predicting survival outcomes and competing risks using statistical and machine learning methods.

## Installation

```r
# Install from CRAN (when available)
install.packages("ml4time2event")

# Or install from GitHub
# install.packages("devtools")
devtools::install_github("dakdemir-nmdp/ml4time2event")
```

## Quick Start

```r
library(ml4time2event)
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
  ensemble_method = "average"
)
```

## Features

- **Survival Models**: Cox, Random Forest, XGBoost, GAM, BART, DeepSurv, GLMNet, GBM, RuleFit
- **Competing Risks**: Fine-Gray, cause-specific Cox, and ML models for competing risks
- **Ensemble Methods**: Averaging, weighted averaging, super learner stacking
- **Metrics**: C-index, Brier score, integrated metrics, expected time lost
- **Comprehensive Testing**: 1900+ passing tests

## Vignettes

- `vignettes/comprehensive_survival_analysis.R` - Complete survival analysis workflow
- `vignettes/comprehensive_competing_risks_analysis.R` - Competing risks analysis

## License

GPL (>= 2)

## Author

Deniz Akdemir (dakdemir@nmdp.org)
