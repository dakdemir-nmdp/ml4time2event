# Bone Marrow Transplant Competing Risks Dataset

## Description

The Bone Marrow Transplant Competing Risks (BMTCRR) dataset contains data from patients who received stem-cell transplants for acute leukemia. This dataset is specifically designed for competing risks analysis, where patients can experience one of two competing events: relapse or treatment-related mortality (TRM).

## Source

- **Package**: casebase
- **Original Study**: Bone Marrow Transplant study
- **Downloaded**: October 2025
- **R Package Version**: casebase 0.10.6

## Data Characteristics

- **Size**: 177 observations
- **Variables**: 7
- **Event Types**:
  - Event 1: Relapse
  - Event 2: Treatment-Related Mortality (TRM)
- **Follow-up**: Time in months

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| Sex          | Factor  | Sex (F = female, M = male)                           |
| D            | Factor  | Disease type (ALL = acute lymphoblastic leukemia, AML = acute myeloid leukemia) |
| Phase        | Factor  | Disease phase at transplant (CR1, CR2, CR3, Relapse) |
| Age          | Integer | Age in years at transplant                           |
| Status       | Integer | Event indicator (0 = censored, 1 = relapse, 2 = TRM) |
| Source       | Factor  | Stem cell source (BM+PB = bone marrow + peripheral blood, PB = peripheral blood only) |
| ftime        | Numeric | Follow-up time in months                             |

## Event Indicator (Status Variable)

- **0**: Censored (patient was event-free at last follow-up)
- **1**: Relapse (event of primary interest)
- **2**: Treatment-Related Mortality (TRM) - competing risk event

This is the standardized format for competing risks analysis where:
- Status = 0: No event occurred (censored)
- Status = 1: Event of interest occurred (relapse)
- Status = 2: Competing event occurred (TRM - death from causes other than relapse)

## Summary Statistics

- Total observations: 177
- Censored: 46 (26.0%)
- Relapse events: 56 (31.6%)
- TRM events: 75 (42.4%)
- Total events: 131 (74.0%)

## Disease Phase Categories

- **CR1**: First complete remission (best prognosis)
- **CR2**: Second complete remission
- **CR3**: Third complete remission
- **Relapse**: Active disease at time of transplant (worst prognosis)

## Stem Cell Source

- **BM+PB**: Bone marrow + peripheral blood stem cells
- **PB**: Peripheral blood stem cells only

## Usage

This dataset is particularly useful for:
- Demonstrating competing risks regression models
- Comparing Fine-Gray subdistribution hazard models vs. cause-specific hazard models
- Illustrating cumulative incidence functions (CIF)
- Analyzing the effect of disease phase and stem cell source on transplant outcomes
- Teaching competing risks concepts
- Flexible parametric modeling with case-base sampling

## Processing Steps

1. Loaded from casebase package using `data(bmtcrr)`
2. Event indicator already in standardized format (0 = censored, 1 = relapse, 2 = TRM)
3. Saved as CSV in inst/extdata/bmtcrr_competing_risks.csv

## References

Scrucca, L., Santucci, A., & Aversa, F. (2010). Regression modeling of competing risk using R: an in depth guide for clinicians. *Bone Marrow Transplantation*, 45(9), 1388-1395.

Klein, J. P., & Moeschberger, M. L. (2003). *Survival Analysis: Techniques for Censored and Truncated Data* (2nd ed.). Springer.

## Example Usage in R

```r
# Load the dataset
bmtcrr <- read.csv(system.file("extdata", "bmtcrr_competing_risks.csv",
                               package = "ml4time2event"))

# Install required packages
library(survival)
library(cmprsk)  # For competing risks analysis

# Create competing risks survival object
# For cmprsk package
cif_obj <- cuminc(ftime = bmtcrr$ftime,
                  fstatus = bmtcrr$Status,
                  group = bmtcrr$Phase)
plot(cif_obj, main = "Cumulative Incidence by Disease Phase")

# Fine-Gray subdistribution hazard model (for relapse)
library(survival)
# Create expanded dataset for Fine-Gray model
bmtcrr$Status_relapse <- ifelse(bmtcrr$Status == 1, 1, 0)
bmtcrr$Status_trm <- ifelse(bmtcrr$Status == 2, 1, 0)

# Fine-Gray model for relapse (event of interest = 1)
fg_fit <- finegray(Surv(ftime, Status, type = "mstate") ~ .,
                   data = bmtcrr, etype = 1)

# Cause-specific Cox model for relapse
cox_relapse <- coxph(Surv(ftime, Status == 1) ~ Sex + D + Phase + Age + Source,
                     data = bmtcrr)
summary(cox_relapse)

# Cause-specific Cox model for TRM
cox_trm <- coxph(Surv(ftime, Status == 2) ~ Sex + D + Phase + Age + Source,
                 data = bmtcrr)
summary(cox_trm)

# Compare Kaplan-Meier by disease phase
library(survminer)
fit <- survfit(Surv(ftime, Status > 0) ~ Phase, data = bmtcrr)
ggsurvplot(fit, data = bmtcrr, pval = TRUE,
           risk.table = TRUE, conf.int = TRUE)
```

## Important Notes for Competing Risks Analysis

1. **Cumulative Incidence Function (CIF)**: Use CIF instead of Kaplan-Meier to properly account for competing risks
2. **Fine-Gray Model**: Models the subdistribution hazard for the event of interest (relapse)
3. **Cause-Specific Models**: Can fit separate Cox models for each event type
4. **Interpretation**: Hazard ratios from Fine-Gray models are interpreted on the subdistribution hazard scale

## Key Clinical Insights

- Disease phase is a strong predictor of outcomes
- Patients in CR1 have the best prognosis
- Patients transplanted during relapse have the worst outcomes
- Both relapse and TRM are important competing endpoints
- Age and stem cell source may also influence outcomes
