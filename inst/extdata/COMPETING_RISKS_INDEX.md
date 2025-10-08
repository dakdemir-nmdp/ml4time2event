# Competing Risks Datasets Index

This document provides an overview of the competing risks datasets included in the ml4time2event package. All datasets have been standardized to use the format: 0 = censored, 1 = event of interest, 2 = competing risk event.

## What is Competing Risks?

Competing risks occur when subjects can experience one of several mutually exclusive events. For example:
- A bone marrow transplant patient can experience **relapse** OR **treatment-related mortality** (but not both)
- A melanoma patient can die from **melanoma** OR from **other causes** (but not both)

Traditional survival analysis methods (like Kaplan-Meier) can produce biased estimates when competing risks are present, because they treat competing events as if they were censored observations.

## Available Competing Risks Datasets

---

### 1. BMT Competing Risks Dataset (Medium)
- **File**: `bmtcrr_competing_risks.csv`
- **Documentation**: [BMTCRR_README.md](BMTCRR_README.md)
- **Size**: 177 observations
- **Variables**: 7
- **Event Types**:
  - Event 1 (Status = 1): Relapse (31.6%)
  - Event 2 (Status = 2): Treatment-Related Mortality (42.4%)
  - Censored (Status = 0): Event-free (26.0%)
- **Source**: Bone Marrow Transplant study (casebase package)
- **Key Variables**: Sex, D (disease type), Phase, Age, Status, Source, ftime

**Best for**:
- Teaching competing risks concepts
- Analyzing transplant outcomes
- Comparing Fine-Gray vs. cause-specific models
- Clinical hematology/oncology research

**Key Features**:
- High event rate (74% total events)
- Two competing events of similar clinical importance
- Multiple prognostic factors (disease phase, age, stem cell source)
- Well-balanced event types

---

### 2. Melanoma Competing Risks Dataset (Medium)
- **File**: `melanoma_competing_risks.csv`
- **Documentation**: [MELANOMA_README.md](MELANOMA_README.md)
- **Size**: 205 observations
- **Variables**: 7
- **Event Types**:
  - Event 1 (status = 1): Death from melanoma (27.8%)
  - Event 2 (status = 2): Death from other causes (6.8%)
  - Censored (status = 0): Alive (65.4%)
- **Source**: University Hospital of Odense, Denmark (boot package)
- **Key Variables**: time, status, sex, age, year, thickness, ulcer

**Best for**:
- Classic competing risks demonstrations
- Teaching cumulative incidence functions
- Analyzing tumor prognostic factors
- Clinical oncology research
- Long-term follow-up studies

**Key Features**:
- Classic competing risks dataset widely used in literature
- Important prognostic factors (tumor thickness, ulceration)
- Competing risk (other deaths) less common but clinically relevant
- Long follow-up period (up to 15 years)

---

### 3. SUPPORT Dataset (Large - Single Event)
- **File**: `support_dataset.csv`
- **Documentation**: [SUPPORT_README.md](SUPPORT_README.md)
- **Size**: 9,104 observations
- **Variables**: 34
- **Event Type**: Single event (death)
- **Source**: SUPPORT study (casebase package)

**Note**: While SUPPORT is primarily a single-event dataset, it is included here because:
1. It can be used with cause-specific models when stratifying by disease group
2. It's excellent for high-dimensional variable selection in survival/competing risks contexts
3. Different disease groups (dzclass) could theoretically be treated as competing risks in specialized analyses

**Best for**:
- High-dimensional survival analysis
- Variable selection with penalized methods
- Machine learning applications to survival data
- Clinical prediction modeling

---

## Dataset Selection Guide

### By Size
- **Large (>5,000)**: SUPPORT
- **Medium (150-250)**: Melanoma, BMT

### By Number of Events
- **High event rate (>70%)**: BMT (74% total events)
- **Moderate event rate (30-40%)**: Melanoma (34.6% total events)
- **Very high event rate**: SUPPORT (68.1% single event)

### By Competing Risk Balance
- **Well-balanced competing risks**: BMT (31.6% vs 42.4%)
- **Primary event dominant**: Melanoma (27.8% vs 6.8%)

### By Application
- **Teaching competing risks basics**: Melanoma (classic dataset)
- **Clinical hematology/transplant**: BMT
- **High-dimensional analysis**: SUPPORT
- **Tumor prognostics**: Melanoma
- **Treatment comparison**: BMT

### By Number of Variables
- **Simple (≤8 variables)**: BMT, Melanoma
- **Complex (>30 variables)**: SUPPORT

## Event Indicator Standardization

All competing risks datasets use the following standardized format:

```
Status = 0: Censored (no event occurred during follow-up)
Status = 1: Event of interest occurred (e.g., relapse, death from melanoma)
Status = 2: Competing risk event occurred (e.g., TRM, death from other causes)
Status = 3+: Additional competing risks (if applicable)
```

This standardization is critical for:
- Consistent analysis across datasets
- Proper use of competing risks methods
- Clear interpretation of results

## Original vs. Standardized Formats

### BMT Dataset
- **Original**: Already standardized (0, 1, 2)
- **Standardized**: No change needed

### Melanoma Dataset
- **Original**: 1 = died from melanoma, 2 = alive, 3 = died from other causes
- **Standardized**: 0 = censored/alive, 1 = died from melanoma, 2 = died from other causes

## Loading Datasets in R

```r
# Load BMT competing risks dataset
bmtcrr <- read.csv(
  system.file("extdata", "bmtcrr_competing_risks.csv", package = "ml4time2event")
)

# Load Melanoma competing risks dataset
melanoma <- read.csv(
  system.file("extdata", "melanoma_competing_risks.csv", package = "ml4time2event")
)

# Load SUPPORT dataset (for high-dimensional analysis)
support <- read.csv(
  system.file("extdata", "support_dataset.csv", package = "ml4time2event")
)
```

## Required R Packages for Competing Risks Analysis

```r
# Install required packages
install.packages(c("survival", "cmprsk", "mstate", "riskRegression"))

# For plotting
install.packages(c("survminer", "ggplot2"))

# For machine learning with competing risks
install.packages(c("randomForestSRC", "pec"))
```

## Common Competing Risks Methods

### 1. Cumulative Incidence Function (CIF)

```r
library(cmprsk)

# Calculate CIF
cif <- cuminc(ftime = bmtcrr$ftime,
              fstatus = bmtcrr$Status,
              group = bmtcrr$Phase)
plot(cif)
```

### 2. Fine-Gray Subdistribution Hazard Model

```r
library(cmprsk)

# Fine-Gray model for event type 1
fg_model <- crr(ftime = bmtcrr$ftime,
                fstatus = bmtcrr$Status,
                cov1 = model.matrix(~ Sex + Age + Phase, data = bmtcrr)[,-1],
                failcode = 1)
summary(fg_model)
```

### 3. Cause-Specific Cox Models

```r
library(survival)

# Model for event type 1
cox1 <- coxph(Surv(ftime, Status == 1) ~ Sex + Age + Phase, data = bmtcrr)

# Model for event type 2
cox2 <- coxph(Surv(ftime, Status == 2) ~ Sex + Age + Phase, data = bmtcrr)
```

### 4. Multi-State Models

```r
library(mstate)

# Define transitions
tmat <- transMat(x = list(c(2, 3), c(), c()), names = c("Initial", "Event1", "Event2"))

# Fit multi-state model
# (Requires data in long format)
```

## Key Differences: Competing Risks vs. Standard Survival

| Aspect | Standard Survival | Competing Risks |
|--------|------------------|-----------------|
| **Events** | One event type | Multiple mutually exclusive events |
| **Censoring** | Only administrative/loss to follow-up | Same, but competing events are NOT censored |
| **Primary Method** | Kaplan-Meier, Cox model | Cumulative Incidence, Fine-Gray, Cause-Specific Cox |
| **Interpretation** | Probability of event in hypothetical world with no competing risks | Actual probability accounting for competing risks |
| **Bias** | Overestimates probability when competing risks exist | Unbiased estimates |

## When to Use Competing Risks Methods

✅ Use competing risks methods when:
- Multiple mutually exclusive event types exist
- Competing events occur in >5-10% of subjects
- Accurate risk prediction is important
- Clinical interpretation requires accounting for competing events

❌ Standard survival methods may be adequate when:
- Only one event type of interest
- Competing events are very rare (<5%)
- Focus is solely on causal effects (but be cautious)

## Common Pitfalls to Avoid

1. ❌ **Using Kaplan-Meier for competing risks**: This will overestimate event probabilities
2. ❌ **Treating competing events as censored**: They are not censored; they prevent the event of interest
3. ❌ **Confusing Fine-Gray and cause-specific hazards**: They answer different questions
4. ❌ **Ignoring competing risks**: Can lead to biased risk predictions
5. ❌ **Wrong interpretation**: Fine-Gray models predict cumulative incidence; cause-specific models explain biological mechanisms

## Recommended Analysis Workflow

1. **Exploratory Analysis**
   - Plot cumulative incidence functions for each event type
   - Stratify by key covariates
   - Compare event rates across groups

2. **Model Selection**
   - Use Fine-Gray for prediction/risk stratification
   - Use cause-specific Cox for understanding risk factors
   - Consider both when possible

3. **Model Checking**
   - Check proportional subdistribution hazards assumption (Fine-Gray)
   - Check proportional hazards assumption (cause-specific Cox)
   - Validate predictions using proper scoring rules

4. **Reporting**
   - Report cumulative incidence curves
   - Present both subdistribution and cause-specific hazard ratios when relevant
   - Clearly state which model was used and why

## References

### Key Papers on Competing Risks

- Austin, P. C., et al. (2016). "Practical recommendations for reporting Fine-Gray model analyses for competing risk data." *Statistics in Medicine*, 35(8), 1083-1097.

- Putter, H., et al. (2007). "Tutorial in biostatistics: competing risks and multi-state models." *Statistics in Medicine*, 26(11), 2389-2430.

- Fine, J. P., & Gray, R. J. (1999). "A proportional hazards model for the subdistribution of a competing risk." *Journal of the American Statistical Association*, 94(446), 496-509.

- Andersen, P. K., et al. (2012). "Competing risks in epidemiology: possibilities and pitfalls." *International Journal of Epidemiology*, 41(3), 861-870.

### Books

- Pintilie, M. (2006). *Competing Risks: A Practical Perspective*. John Wiley & Sons.

- Klein, J. P., & Moeschberger, M. L. (2003). *Survival Analysis: Techniques for Censored and Truncated Data* (2nd ed.). Springer.

## Quality Assurance

All competing risks datasets have undergone:
1. ✅ Event indicators verified and standardized (0 = censored, 1 = event 1, 2 = event 2)
2. ✅ Time variables checked for non-negative values
3. ✅ CSV files verified for proper encoding
4. ✅ Comprehensive documentation created
5. ✅ Example code tested and verified

---

**Last Updated**: October 2025

For questions or issues, please file an issue on the package GitHub repository.
