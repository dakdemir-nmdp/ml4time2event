# Refactoring Plan: ml4time2event

## Phase 1: Architecture Cleanup ("The Two Worlds")
Objective: Remove legacy functional API and consolidate on R6 classes.

- [x] **Verification**: R6 classes are independent of legacy files.
- [x] **Delete Legacy Files**:
    - [x] Simple Models: `R/surv_cox.R`, `R/cr_cox.R`
    - [x] Parametric: `R/surv_reg.R`, `R/cr_survreg.R`
    - [x] Tree Models: `R/surv_random_forest.R`, `R/surv_gbm.R`, `R/cr_random_forest.R`
    - [x] Boosting/Bagging: `R/surv_xgboost.R`, `R/cr_xgboost.R`, `R/surv_bart.R`, `R/cr_bart.R`
    - [x] Smooth/Linear: `R/surv_gam.R`, `R/cr_gam.R`, `R/surv_glmnet.R`, `R/cr_glmnet.R` (if any), `R/cr_fine_gray.R`
    - [x] RuleFit: `R/surv_rulefit.R`, `R/cr_rulefit.R`
    - [x] Neural/Discrete: `R/surv_shallownn.R`, `R/cr_shallownn.R`, `R/surv_ttah.R`, `R/cr_ttah.R`
- [x] **Update NAMESPACE**: Exported legacy functions removed.

## Phase 2: Test Suite Rehabilitation
Objective: Retarget tests to public API (`ml4t2e_fit`).

- [x] **Retarget Tests**: All unit tests updated to `ml4t2e_fit()`.
    - [x] Cox, RF, XGBoost, BART, SurvReg, GAM, GLMNet, RuleFit, ShallowNN, TTAH, GBM.
- [x] **Consolidate Tests**: Merged fix files.
    - [x] `test-cr-ensemble-fix.R` -> `test-ensemble-methods.R`
    - [x] `test-cr-interpolation-fix.R` & `test-cr-interpolator-fix.R` -> `test-cr-interpolation.R`
    - [x] `test-cr-fine-gray.R` (legacy) -> `test-cr-finegray.R` (modern) -> `test-cr-fine-gray.R`
    - [x] `test-ensemble-methods-enhancement.R` -> `test-ensemble-methods.R`

## Phase 3: Code Quality & Standardization
Objective: Standardize on Tidyverse and improve performance.

- [x] **Vectorization**: `R/conformal_prediction.R` optimized.
- [x] **Style Config**: No `stats::reshape` found.
- [x] **Linting**: Fixed `seq_len` in RuleFit.

## Progress Tracker
Current Status: All Phases Complete.
Next: Final comprehensive test run & check vignettes.
