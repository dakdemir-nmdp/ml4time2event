# Veteran's Administration Lung Cancer Dataset

## Description

The Veteran's Administration Lung Cancer Trial dataset contains survival data from a randomized trial of two treatment regimens for lung cancer. This is a classic dataset used in survival analysis literature.

## Source

- **Package**: survival
- **Original Study**: Veterans' Administration Lung Cancer Trial
- **Downloaded**: October 2025
- **R Package Version**: survival 3.5+

## Data Characteristics

- **Size**: 137 observations
- **Variables**: 8
- **Event Type**: Death
- **Follow-up**: Survival time in days

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| trt          | Numeric | Treatment (1 = standard, 2 = test)                   |
| celltype     | Factor  | Cell type (squamous, smallcell, adeno, large)        |
| time         | Numeric | Survival time in days                                |
| status       | Numeric | Event indicator (0 = censored, 1 = death)            |
| karno        | Numeric | Karnofsky performance score (0-100)                  |
| diagtime     | Numeric | Months from diagnosis to randomization               |
| age          | Numeric | Age in years                                         |
| prior        | Numeric | Prior therapy (0 = no, 10 = yes)                     |

## Event Indicator

- **0**: Censored (subject was alive at last follow-up)
- **1**: Event (death)

## Summary Statistics

- Total observations: 137
- Events (deaths): 128
- Censored: 9
- Event rate: 93.4%
- Mean survival time: ~122 days

## Cell Type Distribution

- Squamous: Most common cell type
- Small cell: Typically associated with worse prognosis
- Adeno: Adenocarcinoma
- Large: Large cell carcinoma

## Usage

This dataset is particularly useful for:
- Demonstrating treatment comparisons in survival analysis
- Illustrating stratified analyses by cell type
- Teaching Cox proportional hazards models with categorical predictors
- Analyzing the effect of performance status on survival

## Processing Steps

1. Loaded from survival package using `data(veteran)`
2. Event indicator already standardized (0 = censored, 1 = event)
3. Saved as CSV in inst/extdata/veteran_dataset.csv

## References

Kalbfleisch, J. D., & Prentice, R. L. (2002). *The Statistical Analysis of Failure Time Data* (2nd ed.). Wiley.

## Example Usage in R

```r
# Load the dataset
veteran <- read.csv(system.file("extdata", "veteran_dataset.csv", package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = veteran$time, event = veteran$status)

# Kaplan-Meier curves by treatment
km_trt <- survfit(surv_obj ~ trt, data = veteran)
plot(km_trt, col = c("blue", "red"), xlab = "Time (days)",
     ylab = "Survival Probability")
legend("topright", legend = c("Standard", "Test"), col = c("blue", "red"), lty = 1)

# Cox model with treatment and cell type
cox_model <- coxph(surv_obj ~ trt + celltype + karno + age, data = veteran)
summary(cox_model)

# Stratified analysis by cell type
cox_stratified <- coxph(surv_obj ~ trt + strata(celltype) + karno + age,
                        data = veteran)
summary(cox_stratified)
```

## Key Features

- High event rate (93.4%) makes this dataset ideal for demonstrating survival analysis with minimal censoring
- Multiple prognostic factors (cell type, performance status, age)
- Randomized treatment assignment
- Well-studied dataset with known results
