# Bone Marrow Transplant Dataset

## Description

This dataset contains data from bone marrow transplant (BMT) patients. It was originally designed for competing risks analysis but has been adapted for standard survival analysis by combining all event types into a single event indicator.

**Note**: This dataset is used as a substitute for the Framingham Heart Study dataset, which was not readily available in the requested format. The BMT dataset from the casebase package provides a medium-sized survival dataset suitable for demonstration purposes.

## Source

- **Package**: casebase
- **Original Dataset**: bmtcrr (Bone Marrow Transplant Competing Risks)
- **Downloaded**: October 2025
- **R Package Version**: casebase 0.10.6

## Data Characteristics

- **Size**: 177 observations
- **Variables**: 7
- **Event Type**: Any event (originally relapse or treatment-related mortality)
- **Follow-up**: Time in months

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| Sex          | Factor  | Sex (F = female, M = male)                           |
| D            | Factor  | Disease (ALL = acute lymphoblastic leukemia, AML = acute myeloid leukemia) |
| Phase        | Factor  | Disease phase at transplant (CR1, CR2, CR3, Relapse) |
| Age          | Integer | Age in years                                         |
| Status       | Integer | Event indicator (0 = censored, 1 = any event)        |
| Source       | Factor  | Stem cell source (BM+PB = bone marrow + peripheral blood, PB = peripheral blood) |
| ftime        | Numeric | Follow-up time in months                             |

## Event Indicator

- **0**: Censored (subject was event-free at last follow-up)
- **1**: Event (any event occurred - originally relapse or treatment-related mortality)

**Note**: The original dataset distinguished between relapse (Status = 1) and treatment-related mortality (Status = 2). For standard survival analysis without competing risks, we combined both event types into a single event indicator (Status > 0 → 1).

## Summary Statistics

- Total observations: 177
- Events: 131
- Censored: 46
- Event rate: 74.0%

## Disease Phase

- CR1: First complete remission
- CR2: Second complete remission
- CR3: Third complete remission
- Relapse: Active disease at transplant

## Usage

This dataset is particularly useful for:
- Demonstrating survival analysis in hematologic malignancies
- Analyzing the effect of disease phase on outcomes
- Comparing outcomes by stem cell source
- Teaching survival analysis with medium-sized datasets
- Illustrating stratified survival analysis

## Processing Steps

1. Loaded from casebase package using `data(bmtcrr)`
2. Event indicator recoded: original Status (0, 1, 2) → standardized (0 = censored, 1 = any event)
   - Original Status = 0 (censored) → 0 (censored)
   - Original Status = 1 (relapse) → 1 (event)
   - Original Status = 2 (treatment-related mortality) → 1 (event)
3. Saved as CSV in inst/extdata/framingham_dataset.csv

## References

Scrucca, L., Santucci, A., & Aversa, F. (2010). Regression modeling of competing risk using R: an in depth guide for clinicians. *Bone Marrow Transplantation*, 45(9), 1388-1395.

## Example Usage in R

```r
# Load the dataset
bmt <- read.csv(system.file("extdata", "framingham_dataset.csv", package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = bmt$ftime, event = bmt$Status)

# Kaplan-Meier curves by disease phase
km_phase <- survfit(surv_obj ~ Phase, data = bmt)
plot(km_phase, col = 1:4, xlab = "Time (months)",
     ylab = "Event-Free Survival Probability")
legend("topright", legend = levels(bmt$Phase), col = 1:4, lty = 1)

# Cox model with multiple predictors
cox_model <- coxph(surv_obj ~ Sex + D + Phase + Age + Source, data = bmt)
summary(cox_model)

# Compare by stem cell source
km_source <- survfit(surv_obj ~ Source, data = bmt)
plot(km_source, col = c("blue", "red"), xlab = "Time (months)",
     ylab = "Event-Free Survival Probability")
legend("topright", legend = c("BM+PB", "PB"), col = c("blue", "red"), lty = 1)
```

## Important Note

This dataset is stored as `framingham_dataset.csv` for organizational purposes, but it contains bone marrow transplant data, not Framingham Heart Study data. The Framingham dataset referenced in the original documentation was not readily available in a suitable format for this package.
