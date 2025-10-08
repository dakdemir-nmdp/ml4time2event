# SUPPORT Dataset

## Description

The Study to Understand Prognoses Preferences Outcomes and Risks of Treatment (SUPPORT) dataset contains data from a study of hospitalized adults with serious illnesses. This dataset is commonly used to demonstrate variable selection and high-dimensional survival analysis.

## Source

- **Package**: casebase
- **Original Study**: Study to Understand Prognoses and Preferences for Outcomes and Risks of Treatments (SUPPORT)
- **Downloaded**: October 2025
- **R Package Version**: casebase 0.10.6

## Data Characteristics

- **Size**: 9,104 observations
- **Variables**: 34
- **Event Type**: Death
- **Follow-up**: Variable follow-up time (measured in days)

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| age          | Numeric | Age in years                                         |
| death        | Integer | Event indicator (0 = censored, 1 = death)            |
| sex          | Factor  | Sex (female, male)                                   |
| slos         | Integer | Days from study entry to discharge                   |
| d.time       | Integer | Days of follow-up                                    |
| dzgroup      | Factor  | Disease group (8 categories)                         |
| dzclass      | Factor  | Disease class (ARF/MOSF, Cancer, COPD/CHF/Cirrhosis, Coma) |
| num.co       | Integer | Number of comorbidities                              |
| edu          | Integer | Years of education                                   |
| scoma        | Integer | SUPPORT coma score                                   |
| avtisst      | Numeric | Average TISS score                                   |
| race         | Factor  | Race (asian, black, hispanic, other, white)          |
| hday         | Integer | Day in hospital at study admission                   |
| diabetes     | Integer | Diabetes indicator (0/1)                             |
| dementia     | Integer | Dementia indicator (0/1)                             |
| ca           | Factor  | Cancer status (metastatic, no, yes)                  |
| meanbp       | Numeric | Mean blood pressure                                  |
| wblc         | Numeric | White blood cell count                               |
| hrt          | Numeric | Heart rate                                           |
| resp         | Integer | Respiratory rate                                     |
| temp         | Numeric | Temperature                                          |
| pafi         | Numeric | PaO2/(FiO2*100)                                      |
| alb          | Numeric | Albumin                                              |
| bili         | Numeric | Bilirubin                                            |
| crea         | Numeric | Creatinine                                           |
| sod          | Integer | Sodium                                               |
| ph           | Numeric | pH                                                   |
| glucose      | Numeric | Glucose                                              |
| bun          | Numeric | Blood urea nitrogen                                  |
| urine        | Numeric | Urine output                                         |
| adlp         | Integer | ADL patient score                                    |
| adlsc        | Numeric | ADL surrogate score                                  |
| sps          | Numeric | SUPPORT physiology score                             |
| aps          | Integer | APACHE III score                                     |

## Event Indicator

- **0**: Censored (subject was alive at last follow-up)
- **1**: Event (death)

## Summary Statistics

- Total observations: 9,104
- Events (deaths): 6,200
- Censored: 2,904
- Event rate: 68.1%

## Usage

This dataset is particularly useful for:
- Variable selection in high-dimensional survival analysis
- Regularized hazard estimation
- Demonstrating penalized survival models
- Clinical prediction modeling

## Processing Steps

1. Loaded from casebase package using `data(support)`
2. Event indicator already standardized (0 = censored, 1 = event)
3. Some observations removed due to missingness (as noted in original data)
4. Saved as CSV in inst/extdata/support_dataset.csv

## References

Knaus, W. A., Harrell, F. E., Lynn, J., et al. (1995). The SUPPORT prognostic model: Objective estimates of survival for seriously ill hospitalized adults. *Annals of Internal Medicine*, 122(3), 191-203.

## Example Usage in R

```r
# Load the dataset
support <- read.csv(system.file("extdata", "support_dataset.csv", package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = support$d.time, event = support$death)

# Fit Cox model with multiple predictors
cox_model <- coxph(surv_obj ~ age + sex + dzclass + num.co + edu +
                   meanbp + hrt + resp + temp, data = support)
summary(cox_model)

# Penalized Cox model (lasso)
library(glmnet)
x <- model.matrix(~ age + sex + dzclass + num.co + edu + meanbp +
                  hrt + resp + temp - 1, data = support)
y <- Surv(support$d.time, support$death)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1)
```
