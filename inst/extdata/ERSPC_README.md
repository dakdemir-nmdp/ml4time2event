# ERSPC Dataset

## Description

The European Randomized Study of Prostate Cancer Screening (ERSPC) dataset contains survival data from a large multi-center randomized trial evaluating the effect of prostate-specific antigen (PSA) screening on prostate cancer mortality.

## Source

- **Package**: casebase
- **Original Study**: European Randomized Study of Prostate Cancer Screening (ERSPC)
- **Downloaded**: October 2025
- **R Package Version**: casebase 0.10.6

## Data Characteristics

- **Size**: 159,893 observations
- **Variables**: 3
- **Event Type**: Death from prostate cancer
- **Follow-up**: Variable follow-up time across participants

## Variables

| Variable Name    | Type    | Description                                          |
|-----------------|---------|------------------------------------------------------|
| ScrArm          | Factor  | Screening arm assignment (Control group, Screening group) |
| Follow.Up.Time  | Numeric | Follow-up time (in years)                            |
| DeadOfPrCa      | Integer | Event indicator (0 = censored, 1 = death from prostate cancer) |

## Event Indicator

- **0**: Censored (subject did not die from prostate cancer during follow-up)
- **1**: Event (death from prostate cancer)

## Summary Statistics

- Total observations: 159,893
- Events (deaths): 540
- Censored: 159,353
- Event rate: 0.34%

## Usage

This dataset is particularly useful for:
- Demonstrating flexible hazard function modeling
- Analyzing the effect of screening interventions on long-term survival
- Illustrating survival analysis with a single event type (no competing risks)
- Modeling time-to-event data with low event rates

## Processing Steps

1. Loaded from casebase package using `data(ERSPC)`
2. Event indicator already standardized (0 = censored, 1 = event)
3. Saved as CSV in inst/extdata/erspc_dataset.csv

## References

Schröder, F. H., Hugosson, J., Roobol, M. J., et al. (2009). Screening and prostate-cancer mortality in a randomized European study. *New England Journal of Medicine*, 360(13), 1320-1328.

## Example Usage in R

```r
# Load the dataset
erspc <- read.csv(system.file("extdata", "erspc_dataset.csv", package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = erspc$Follow.Up.Time, event = erspc$DeadOfPrCa)

# Fit Cox proportional hazards model
cox_model <- coxph(surv_obj ~ ScrArm, data = erspc)
summary(cox_model)
```
