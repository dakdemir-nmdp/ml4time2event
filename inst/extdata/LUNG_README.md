# Lung Cancer Dataset

## Description

The lung cancer dataset from the North Central Cancer Treatment Group (NCCTG) contains survival data from patients with advanced lung cancer. This is one of the most commonly used datasets for teaching survival analysis.

## Source

- **Package**: survival
- **Original Study**: NCCTG Lung Cancer Trial
- **Downloaded**: October 2025
- **R Package Version**: survival 3.5+

## Data Characteristics

- **Size**: 228 observations
- **Variables**: 10
- **Event Type**: Death
- **Follow-up**: Survival time in days

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| inst         | Numeric | Institution code                                     |
| time         | Numeric | Survival time in days                                |
| status       | Numeric | Event indicator (0 = censored, 1 = death)            |
| age          | Numeric | Age in years                                         |
| sex          | Numeric | Sex (1 = male, 2 = female)                           |
| ph.ecog      | Numeric | ECOG performance score (0-4)                         |
| ph.karno     | Numeric | Karnofsky performance score (physician-rated, 0-100) |
| pat.karno    | Numeric | Karnofsky performance score (patient-rated, 0-100)   |
| meal.cal     | Numeric | Calories consumed at meals                           |
| wt.loss      | Numeric | Weight loss in last six months (pounds)              |

## Event Indicator

- **0**: Censored (subject was alive at last follow-up)
- **1**: Event (death)

**Note**: The original status variable was recoded from the survival package format where 1 = censored and 2 = dead. We standardized it to 0 = censored and 1 = dead for consistency with other datasets.

## Summary Statistics

- Total observations: 228
- Events (deaths): 165
- Censored: 63
- Event rate: 72.4%
- Mean survival time: ~305 days
- Some variables contain missing values (NA)

## Usage

This dataset is particularly useful for:
- Teaching basic survival analysis concepts
- Demonstrating Kaplan-Meier curves
- Illustrating Cox proportional hazards models
- Comparing survival between groups (e.g., sex)
- Handling missing data in survival analysis

## Processing Steps

1. Loaded from survival package using `data(lung)`
2. Event indicator recoded: original status (1 = censored, 2 = dead) → standardized (0 = censored, 1 = dead)
3. Saved as CSV in inst/extdata/lung_dataset.csv

## References

Loprinzi, C. L., Laurie, J. A., Wieand, H. S., et al. (1994). Prospective evaluation of prognostic variables from patient-completed questionnaires. *Journal of Clinical Oncology*, 12(3), 601-607.

## Example Usage in R

```r
# Load the dataset
lung <- read.csv(system.file("extdata", "lung_dataset.csv", package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = lung$time, event = lung$status)

# Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ 1, data = lung)
plot(km_fit, xlab = "Time (days)", ylab = "Survival Probability")

# Compare survival by sex
km_sex <- survfit(surv_obj ~ sex, data = lung)
plot(km_sex, col = c("blue", "red"), xlab = "Time (days)",
     ylab = "Survival Probability")
legend("topright", legend = c("Male", "Female"), col = c("blue", "red"), lty = 1)

# Cox proportional hazards model
cox_model <- coxph(surv_obj ~ age + sex + ph.ecog, data = lung)
summary(cox_model)
```

## Missing Data

Several variables contain missing values:
- ph.ecog: 1 missing
- ph.karno: 1 missing
- pat.karno: 3 missing
- meal.cal: 47 missing
- wt.loss: 14 missing

Consider handling missing data appropriately in your analysis (e.g., complete case analysis, imputation).
