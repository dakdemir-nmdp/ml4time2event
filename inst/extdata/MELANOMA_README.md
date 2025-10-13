# Melanoma Competing Risks Dataset

## Description

The Melanoma dataset contains survival data from patients with malignant melanoma who had their tumors surgically removed at the Department of Plastic Surgery, University Hospital of Odense, Denmark during 1962-1977. This is a classic dataset for demonstrating competing risks analysis, where patients can either die from melanoma (event of interest) or die from other unrelated causes (competing risk).

## Source

- **Package**: boot
- **Original Study**: University Hospital of Odense, Denmark (1962-1977)
- **Downloaded**: October 2025
- **R Package Version**: boot 1.3+

## Data Characteristics

- **Size**: 205 observations
- **Variables**: 7
- **Event Types**:
  - Event 1: Death from melanoma (event of interest)
  - Event 2: Death from other causes (competing risk)
- **Follow-up**: Survival time in days since operation

## Variables

| Variable Name | Type    | Description                                          |
|--------------|---------|------------------------------------------------------|
| time         | Numeric | Survival time in days since the operation            |
| status       | Numeric | Event indicator (0 = censored/alive, 1 = died from melanoma, 2 = died from other causes) |
| sex          | Numeric | Sex (1 = male, 0 = female)                           |
| age          | Numeric | Age in years at time of operation                    |
| year         | Numeric | Year of operation                                    |
| thickness    | Numeric | Tumor thickness in millimeters                       |
| ulcer        | Numeric | Ulceration indicator (1 = present, 0 = absent)       |

## Event Indicator (Status Variable)

- **0**: Censored (patient was alive at end of study in 1977)
- **1**: Died from melanoma (event of interest)
- **2**: Died from other causes (competing risk event)

**Note**: The original status variable was recoded from the boot package format:
- Original: 1 = died from melanoma, 2 = alive, 3 = died from other causes
- Standardized: 0 = censored, 1 = died from melanoma, 2 = died from other causes

## Summary Statistics

- Total observations: 205
- Censored (alive): 134 (65.4%)
- Died from melanoma: 57 (27.8%)
- Died from other causes: 14 (6.8%)
- Total deaths: 71 (34.6%)

## Prognostic Variables

### Tumor Thickness
One of the most important prognostic factors in melanoma. Thicker tumors are associated with worse prognosis.

### Ulceration
Ulcerated tumors typically have a worse prognosis than non-ulcerated tumors.

### Age and Sex
Demographic factors that may influence survival outcomes.

## Usage

This dataset is particularly useful for:
- Teaching fundamental competing risks concepts
- Demonstrating Fine-Gray subdistribution hazard models
- Calculating and plotting cumulative incidence functions (CIF)
- Comparing cause-specific hazard models
- Analyzing the effect of tumor characteristics on survival
- Illustrating the difference between competing risks and standard survival analysis
- Clinical oncology research and teaching

## Processing Steps

1. Loaded from boot package using `data(melanoma)`
2. Event indicator recoded to standardized format:
   - Original status 2 (alive) → 0 (censored)
   - Original status 1 (died from melanoma) → 1 (event of interest)
   - Original status 3 (died from other causes) → 2 (competing risk)
3. Saved as CSV in inst/extdata/melanoma_competing_risks.csv

## References

Andersen, P. K., Borgan, O., Gill, R. D., & Keiding, N. (1993). *Statistical Models Based on Counting Processes*. Springer-Verlag.

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their Application*. Cambridge University Press.

Venables, W. N., & Ripley, B. D. (1994). *Modern Applied Statistics with S-Plus*. Springer-Verlag.

## Example Usage in R

```r
# Load the dataset
melanoma <- read.csv(system.file("extdata", "melanoma_competing_risks.csv",
                                 package = "ml4time2event"))

# Load required packages
library(survival)
library(cmprsk)

# ============================================
# 1. Cumulative Incidence Functions
# ============================================

# Calculate cumulative incidence by sex
cif_sex <- cuminc(ftime = melanoma$time,
                  fstatus = melanoma$status,
                  group = melanoma$sex)
plot(cif_sex, main = "Cumulative Incidence by Sex",
     xlab = "Time (days)", ylab = "Cumulative Incidence")

# Calculate cumulative incidence by ulceration status
cif_ulcer <- cuminc(ftime = melanoma$time,
                    fstatus = melanoma$status,
                    group = melanoma$ulcer)
plot(cif_ulcer, main = "Cumulative Incidence by Ulceration",
     xlab = "Time (days)", ylab = "Cumulative Incidence")

# ============================================
# 2. Fine-Gray Subdistribution Hazard Model
# ============================================

# Fine-Gray model for melanoma death (event type 1)
# This models the subdistribution hazard
library(survival)

# Model death from melanoma
fg_model <- crr(ftime = melanoma$time,
                fstatus = melanoma$status,
                cov1 = as.matrix(melanoma[, c("sex", "age", "thickness", "ulcer")]),
                failcode = 1)  # Event of interest is melanoma death
summary(fg_model)

# ============================================
# 3. Cause-Specific Cox Models
# ============================================

# Cause-specific model for melanoma death (status = 1)
cox_melanoma <- coxph(Surv(time, status == 1) ~ sex + age + thickness + ulcer,
                      data = melanoma)
summary(cox_melanoma)

# Cause-specific model for other causes death (status = 2)
cox_other <- coxph(Surv(time, status == 2) ~ sex + age + thickness + ulcer,
                   data = melanoma)
summary(cox_other)

# ============================================
# 4. Stratified Analysis by Tumor Thickness
# ============================================

# Create thickness categories
melanoma$thick_cat <- cut(melanoma$thickness,
                          breaks = c(0, 1, 4, Inf),
                          labels = c("Thin (<1mm)", "Medium (1-4mm)", "Thick (>4mm)"))

# Cumulative incidence by thickness category
cif_thick <- cuminc(ftime = melanoma$time,
                    fstatus = melanoma$status,
                    group = melanoma$thick_cat)
plot(cif_thick, main = "Cumulative Incidence by Tumor Thickness",
     xlab = "Time (days)", ylab = "Cumulative Incidence")

# ============================================
# 5. Compare with Standard Kaplan-Meier
# ============================================

# INCORRECT: Standard KM treating competing risks as censored
km_wrong <- survfit(Surv(time, status == 1) ~ 1, data = melanoma)

# CORRECT: Using competing risks CIF
cif_correct <- cuminc(ftime = melanoma$time,
                      fstatus = melanoma$status)

# Plot both for comparison (shows why competing risks matter)
plot(cif_correct, main = "Correct: Cumulative Incidence",
     xlab = "Time (days)", ylab = "Probability")
```

## Important Methodological Notes

### Why Competing Risks Matters

1. **Kaplan-Meier is biased**: Traditional KM estimates will overestimate the probability of melanoma death because it treats deaths from other causes as censored observations
2. **Cumulative Incidence Function (CIF)**: Properly accounts for competing risks and gives the actual probability of dying from melanoma
3. **Fine-Gray vs. Cox**:
   - Fine-Gray models the subdistribution hazard (useful for prediction)
   - Cause-specific Cox models the cause-specific hazard (useful for understanding risk factors)

### Interpretation of Results

- **Fine-Gray Model**: Hazard ratios represent the effect on the subdistribution hazard, which relates directly to cumulative incidence
- **Cause-Specific Models**: Hazard ratios represent the effect on the instantaneous risk of the specific event among those still at risk
- **Cumulative Incidence**: Gives the actual probability of experiencing each event type over time

## Key Clinical Findings

Based on the melanoma literature using this dataset:

1. **Tumor thickness** is the strongest predictor of melanoma death
2. **Ulceration** is associated with worse prognosis
3. **Age** may influence both melanoma death and death from other causes
4. **Competing risks** become more important with age (older patients more likely to die from other causes)

## Common Mistakes to Avoid

1. ❌ Using Kaplan-Meier and treating competing events as censored
2. ❌ Ignoring competing risks when they are common (>5-10% of events)
3. ❌ Confusing Fine-Gray and cause-specific hazard ratios
4. ❌ Not considering that risk factors may have different effects on different event types

## Additional Resources

For more information on competing risks analysis:
- Austin, P. C., et al. (2016). "Practical recommendations for reporting Fine-Gray model analyses for competing risk data." *Statistics in Medicine*, 35(8), 1083-1097.
- Putter, H., et al. (2007). "Tutorial in biostatistics: competing risks and multi-state models." *Statistics in Medicine*, 26(11), 2389-2430.
