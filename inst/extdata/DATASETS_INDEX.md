# Survival Datasets Index

This directory contains survival analysis datasets included in the ml4time2event package. All datasets have been standardized to use consistent event indicators (0 = censored, 1 = event).

## Available Datasets

### 1. ERSPC Dataset (Very Large)
- **File**: `erspc_dataset.csv`
- **Documentation**: [ERSPC_README.md](ERSPC_README.md)
- **Size**: 159,893 observations
- **Variables**: 3
- **Event Type**: Death from prostate cancer
- **Event Rate**: 0.34%
- **Source**: European Randomized Study of Prostate Cancer Screening
- **Key Variables**: ScrArm, Follow.Up.Time, DeadOfPrCa

**Best for**: Large-scale survival analysis, low event rate scenarios, screening trial analysis

---

### 2. SUPPORT Dataset (Large)
- **File**: `support_dataset.csv`
- **Documentation**: [SUPPORT_README.md](SUPPORT_README.md)
- **Size**: 9,104 observations
- **Variables**: 34
- **Event Type**: Death
- **Event Rate**: 68.1%
- **Source**: Study to Understand Prognoses and Preferences for Outcomes and Risks of Treatments
- **Key Variables**: age, death, sex, d.time, dzclass, plus 29 clinical/physiological variables

**Best for**: High-dimensional survival analysis, variable selection, clinical prediction models

---

### 3. Lung Cancer Dataset (Medium)
- **File**: `lung_dataset.csv`
- **Documentation**: [LUNG_README.md](LUNG_README.md)
- **Size**: 228 observations
- **Variables**: 10
- **Event Type**: Death
- **Event Rate**: 72.4%
- **Source**: NCCTG Lung Cancer Trial
- **Key Variables**: inst, time, status, age, sex, ph.ecog, ph.karno, pat.karno, meal.cal, wt.loss

**Best for**: Teaching survival analysis, demonstrating Kaplan-Meier curves, Cox models

---

### 4. Veteran Dataset (Small)
- **File**: `veteran_dataset.csv`
- **Documentation**: [VETERAN_README.md](VETERAN_README.md)
- **Size**: 137 observations
- **Variables**: 8
- **Event Type**: Death
- **Event Rate**: 93.4%
- **Source**: Veterans' Administration Lung Cancer Trial
- **Key Variables**: trt, celltype, time, status, karno, diagtime, age, prior

**Best for**: Treatment comparison, stratified analysis, high event rate scenarios

---

### 5. BMT Dataset (Medium)
- **File**: `framingham_dataset.csv`
- **Documentation**: [BMT_README.md](BMT_README.md)
- **Size**: 177 observations
- **Variables**: 7
- **Event Type**: Any event (relapse or treatment-related mortality)
- **Event Rate**: 74.0%
- **Source**: Bone Marrow Transplant data (casebase package)
- **Key Variables**: Sex, D, Phase, Age, Status, Source, ftime

**Best for**: Hematologic malignancy analysis, stratified survival analysis by disease phase

**Note**: This dataset is stored as `framingham_dataset.csv` but contains BMT data, not Framingham Heart Study data.

---

## Dataset Selection Guide

### By Size
- **Very Large (>100,000)**: ERSPC
- **Large (>5,000)**: SUPPORT
- **Medium (100-500)**: Lung, BMT
- **Small (<150)**: Veteran

### By Event Rate
- **Very Low (<1%)**: ERSPC
- **Moderate (60-75%)**: SUPPORT, Lung, BMT
- **Very High (>90%)**: Veteran

### By Number of Variables
- **Simple (≤8 variables)**: ERSPC, Veteran, BMT
- **Moderate (10 variables)**: Lung
- **Complex (>30 variables)**: SUPPORT

### By Application
- **Teaching/Demonstration**: Lung, Veteran
- **Variable Selection**: SUPPORT
- **Treatment Comparison**: Veteran
- **Large-scale Analysis**: ERSPC
- **Clinical Oncology**: Lung, Veteran, BMT
- **Screening Trials**: ERSPC

## Event Indicator Standardization

All datasets have been standardized to use the following event indicator convention:
- **0**: Censored (event did not occur during follow-up)
- **1**: Event (event occurred)

This standardization was necessary because different packages use different conventions:
- Some use (1, 2) where 1 = censored and 2 = event
- Some use (0, 1) where 0 = censored and 1 = event
- Some use (0, 1, 2) for competing risks where 0 = censored, 1 = event type 1, 2 = event type 2

## Loading Datasets in R

```r
# General pattern
dataset_name <- read.csv(
  system.file("extdata", "dataset_file.csv", package = "ml4time2event")
)

# Specific examples
erspc <- read.csv(system.file("extdata", "erspc_dataset.csv", package = "ml4time2event"))
support <- read.csv(system.file("extdata", "support_dataset.csv", package = "ml4time2event"))
lung <- read.csv(system.file("extdata", "lung_dataset.csv", package = "ml4time2event"))
veteran <- read.csv(system.file("extdata", "veteran_dataset.csv", package = "ml4time2event"))
bmt <- read.csv(system.file("extdata", "framingham_dataset.csv", package = "ml4time2event"))
```

## Creating Survival Objects

All datasets follow a consistent pattern for creating survival objects:

```r
library(survival)

# ERSPC
surv_erspc <- Surv(time = erspc$Follow.Up.Time, event = erspc$DeadOfPrCa)

# SUPPORT
surv_support <- Surv(time = support$d.time, event = support$death)

# Lung
surv_lung <- Surv(time = lung$time, event = lung$status)

# Veteran
surv_veteran <- Surv(time = veteran$time, event = veteran$status)

# BMT
surv_bmt <- Surv(time = bmt$ftime, event = bmt$Status)
```

## Quality Assurance

All datasets have undergone the following quality checks:
1. Event indicators verified and standardized to (0 = censored, 1 = event)
2. Time variables checked for non-negative values
3. CSV files verified for proper encoding and format
4. Documentation created with source information and processing steps

## Package Integration

These datasets are stored in `inst/extdata/` and can be accessed using `system.file()` as shown above. For package development, you can also document these datasets in `R/data.R` using roxygen2:

```r
#' ERSPC Dataset
#'
#' @description
#' European Randomized Study of Prostate Cancer Screening data
#'
#' @format A data frame with 159,893 rows and 3 variables
#' @source casebase package
"erspc"
```

## References

See individual dataset documentation files for specific references and citations.

## Contact

For questions or issues with these datasets, please file an issue on the package GitHub repository.

---

**Last Updated**: October 2025
