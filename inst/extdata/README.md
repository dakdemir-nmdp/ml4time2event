# ml4time2event Survival Datasets

This directory contains curated survival analysis datasets for the ml4time2event package. All datasets have been standardized for consistent use across survival and competing risks analyses.

## Quick Navigation

### 📊 Dataset Indexes
- **[Standard Survival Datasets Index](DATASETS_INDEX.md)** - Complete guide to single-event survival datasets
- **[Competing Risks Datasets Index](COMPETING_RISKS_INDEX.md)** - Complete guide to competing risks datasets

### 📖 Individual Dataset Documentation

#### Standard Survival Datasets (Single Event)
- [ERSPC Dataset](ERSPC_README.md) - Large prostate cancer screening study (159,893 obs)
- [SUPPORT Dataset](SUPPORT_README.md) - Large clinical study with 34 variables (9,104 obs)
- [Lung Cancer Dataset](LUNG_README.md) - Classic teaching dataset (228 obs)
- [Veteran Dataset](VETERAN_README.md) - VA lung cancer trial (137 obs)
- [BMT Dataset](BMT_README.md) - Bone marrow transplant, single-event version (177 obs)

#### Competing Risks Datasets (Multiple Events)
- [BMT Competing Risks](BMTCRR_README.md) - Relapse vs. TRM (177 obs)
- [Melanoma Competing Risks](MELANOMA_README.md) - Melanoma death vs. other causes (205 obs)

---

## Dataset Overview Table

### Standard Survival Datasets

| Dataset | Size | Variables | Event Rate | Time Unit | Best For |
|---------|------|-----------|------------|-----------|----------|
| **ERSPC** | 159,893 | 3 | 0.34% | Years | Large-scale analysis, low event rate |
| **SUPPORT** | 9,104 | 34 | 68.1% | Days | Variable selection, high-dimensional |
| **Lung** | 228 | 10 | 72.4% | Days | Teaching, demonstrations |
| **Veteran** | 137 | 8 | 93.4% | Days | Treatment comparison, high event rate |
| **BMT** | 177 | 7 | 74.0% | Months | Hematologic malignancies |

### Competing Risks Datasets

| Dataset | Size | Variables | Event 1 | Event 2 | Censored | Time Unit | Best For |
|---------|------|-----------|---------|---------|----------|-----------|----------|
| **BMT CR** | 177 | 7 | 31.6% | 42.4% | 26.0% | Months | Transplant outcomes, balanced events |
| **Melanoma** | 205 | 7 | 27.8% | 6.8% | 65.4% | Days | Classic CR teaching, tumor prognostics |

---

## Event Indicator Standardization

All datasets follow standardized coding conventions:

### Standard Survival (Single Event)
```
0 = Censored (event did not occur)
1 = Event occurred
```

### Competing Risks (Multiple Events)
```
0 = Censored (no event occurred)
1 = Event of interest (e.g., relapse, melanoma death)
2 = Competing risk event (e.g., TRM, other death)
3+ = Additional competing risks (if applicable)
```

---

## Quick Start Examples

### Standard Survival Analysis

```r
# Load a standard survival dataset
lung <- read.csv(system.file("extdata", "lung_dataset.csv",
                             package = "ml4time2event"))

# Create survival object
library(survival)
surv_obj <- Surv(time = lung$time, event = lung$status)

# Kaplan-Meier curve
km_fit <- survfit(surv_obj ~ 1, data = lung)
plot(km_fit)

# Cox proportional hazards model
cox_fit <- coxph(surv_obj ~ age + sex + ph.ecog, data = lung)
summary(cox_fit)
```

### Competing Risks Analysis

```r
# Load a competing risks dataset
bmtcrr <- read.csv(system.file("extdata", "bmtcrr_competing_risks.csv",
                               package = "ml4time2event"))

# Cumulative incidence function
library(cmprsk)
cif <- cuminc(ftime = bmtcrr$ftime,
              fstatus = bmtcrr$Status,
              group = bmtcrr$Phase)
plot(cif)

# Fine-Gray subdistribution hazard model
fg_fit <- crr(ftime = bmtcrr$ftime,
              fstatus = bmtcrr$Status,
              cov1 = model.matrix(~ Sex + Age + Phase, data = bmtcrr)[,-1],
              failcode = 1)
summary(fg_fit)
```

---

## Dataset Selection Decision Tree

```
Do you have competing events?
│
├─ NO → Use Standard Survival Datasets
│   │
│   ├─ Need large sample? → ERSPC (159K obs)
│   ├─ Need many variables? → SUPPORT (34 vars)
│   ├─ Teaching/learning? → Lung (228 obs)
│   ├─ High event rate? → Veteran (93% events)
│   └─ Clinical hematology? → BMT single event
│
└─ YES → Use Competing Risks Datasets
    │
    ├─ Well-balanced events? → BMT CR (relapse vs. TRM)
    ├─ Classic CR teaching? → Melanoma (melanoma vs. other)
    └─ High-dimensional? → SUPPORT (with stratification)
```

---

## File Organization

```
inst/extdata/
├── README.md                          # This file - main overview
├── DATASETS_INDEX.md                  # Standard survival datasets guide
├── COMPETING_RISKS_INDEX.md           # Competing risks datasets guide
│
├── Standard Survival Datasets
│   ├── erspc_dataset.csv              # ERSPC data (159,893 obs)
│   ├── ERSPC_README.md                # ERSPC documentation
│   ├── support_dataset.csv            # SUPPORT data (9,104 obs)
│   ├── SUPPORT_README.md              # SUPPORT documentation
│   ├── lung_dataset.csv               # Lung cancer data (228 obs)
│   ├── LUNG_README.md                 # Lung documentation
│   ├── veteran_dataset.csv            # Veteran data (137 obs)
│   ├── VETERAN_README.md              # Veteran documentation
│   ├── framingham_dataset.csv         # BMT data (177 obs)
│   └── BMT_README.md                  # BMT documentation
│
└── Competing Risks Datasets
    ├── bmtcrr_competing_risks.csv     # BMT CR data (177 obs)
    ├── BMTCRR_README.md               # BMT CR documentation
    ├── melanoma_competing_risks.csv   # Melanoma data (205 obs)
    └── MELANOMA_README.md             # Melanoma documentation
```

---

## Data Sources and Acknowledgments

All datasets are from well-established R packages:

- **casebase package**: ERSPC, SUPPORT, BMT/BMTCRR datasets
- **survival package**: Lung, Veteran datasets
- **boot package**: Melanoma dataset

We gratefully acknowledge the original investigators who collected these data and made them publicly available for research and education.

---

## Required R Packages

### For Standard Survival Analysis
```r
install.packages(c("survival", "survminer", "ggplot2"))
```

### For Competing Risks Analysis
```r
install.packages(c("survival", "cmprsk", "mstate", "riskRegression",
                   "survminer", "ggplot2"))
```

### For Machine Learning with Survival Data
```r
install.packages(c("randomForestSRC", "glmnet", "gbm", "xgboost",
                   "pec", "riskRegression"))
```

---

## Key References

### Survival Analysis
- Klein, J. P., & Moeschberger, M. L. (2003). *Survival Analysis: Techniques for Censored and Truncated Data* (2nd ed.). Springer.
- Therneau, T. M., & Grambsch, P. M. (2000). *Modeling Survival Data: Extending the Cox Model*. Springer.

### Competing Risks
- Pintilie, M. (2006). *Competing Risks: A Practical Perspective*. John Wiley & Sons.
- Austin, P. C., et al. (2016). "Practical recommendations for reporting Fine-Gray model analyses for competing risk data." *Statistics in Medicine*, 35(8), 1083-1097.

### Machine Learning for Survival
- Ishwaran, H., et al. (2008). "Random survival forests." *The Annals of Applied Statistics*, 2(3), 841-860.
- Katzman, J. L., et al. (2018). "DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network." *BMC Medical Research Methodology*, 18(1), 24.

---

## Quality Assurance

All datasets have been:
- ✅ Downloaded from reputable R packages
- ✅ Verified for data integrity
- ✅ Standardized for event indicators
- ✅ Documented with source information
- ✅ Tested with example code

---

## Usage in Package

To use these datasets in package vignettes or examples:

```r
# Standard survival
data_path <- system.file("extdata", "lung_dataset.csv",
                         package = "ml4time2event")
lung <- read.csv(data_path)

# Competing risks
data_path <- system.file("extdata", "bmtcrr_competing_risks.csv",
                         package = "ml4time2event")
bmtcrr <- read.csv(data_path)
```

---

## Contributing

If you have suggestions for additional datasets or find issues with existing datasets, please file an issue on the package GitHub repository.

---

## License

These datasets are included for educational and research purposes. Please cite the original sources when using these datasets in publications. See individual dataset documentation for specific references.

---

**Last Updated**: October 2025

**Package**: ml4time2event

**Maintainer**: [Package maintainer information]
