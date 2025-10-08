# Aalen-Johansen CIF Calibration Verification Report

## Summary
✅ **VERIFICATION COMPLETE**: All cause-specific competing risks models are correctly implementing the Aalen-Johansen estimator instead of the biased `1 - S(t)` approximation.

## Detailed Review

### ✅ Models Using Correct Aalen-Johansen Implementation

#### 1. **Cox Model** (`cr_cox.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Fits separate Cox models for all competing events, uses `aalenJohansenCIF()`
- **Line 397**: `cif_matrix <- aalenJohansenCIF(cause_specific_survs, times, event_of_interest = failcode)`

#### 2. **GAM Model** (`cr_gam.R`) 
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Fits separate GAM models for all competing events, uses `aalenJohansenCIF()`
- **Line 450**: `cif_matrix <- aalenJohansenCIF(cause_specific_survs, times, event_of_interest = failcode)`

#### 3. **XGBoost Model** (`cr_xgboost.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Fits separate XGBoost models for all competing events, uses `aalenJohansenCIF()`
- **Line 371**: `cif_matrix <- aalenJohansenCIF(cause_specific_survs, times, event_of_interest = failcode)`

#### 4. **SurvReg Model** (`cr_survreg.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED  
- **Method**: Fits separate SurvReg models for all competing events, uses `aalenJohansenCIF()`
- **Line 450**: `cif_matrix <- aalenJohansenCIF(cause_specific_survs, times, event_of_interest = failcode)`

### ✅ Models Using Alternative Correct Methods

#### 5. **Fine-Gray Model** (`cr_fine_gray.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Uses subdistribution hazards directly - mathematically correct for Fine-Gray
- **Line 295**: `cif_vals <- 1 - exp(-cumsum(cumhaz))` (Correct for subdistribution hazards)

#### 6. **BART Model** (`cr_bart.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Uses BART package's internal competing risks implementation
- **Line 250**: Uses `BART::crisk.pre.bart()` and `predict()` - package handles CIF correctly

#### 7. **Random Forest Model** (`cr_random_forest.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Uses randomForestSRC's built-in CIF calculation
- **Line 218**: `cif_matrix <- pred_rf$cif[, , paste0("CIF.", modelout$failcode)]` - package handles CIF correctly

#### 8. **RuleFit Model** (`cr_rulefit.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Uses Fine-Gray internally via `Predict_CRModel_FineGray()`
- **Line 358**: Delegates to Fine-Gray which is mathematically correct

#### 9. **DeepSurv Model** (`cr_deepsurv.R`)
- **Status**: ✅ CORRECTLY IMPLEMENTED
- **Method**: Uses Fine-Gray style subdistribution hazards
- **Line 479**: `cif_matrix <- 1 - exp(-outer(baseline_cumsubhaz, exp(corrected_risk_scores)))` (Correct for subdistribution)

### ✅ Aalen-Johansen Implementation Details

The `aalenJohansenCIF()` function in `cr_interpolation.R` correctly implements the mathematical formula:

```
1. Extract cause-specific hazards from survival functions
2. Calculate overall survival: S_overall(t) = ∏_k S_k(t)  
3. Integrate: CIF_j(t) = ∫_0^t λ_j(s) S_overall(s^-) ds
```

**Key Features:**
- Properly handles competing risks by considering ALL causes
- Avoids the biased `1 - S_cause_specific(t)` approximation
- Uses numerical integration for practical computation
- Ensures proper probability bounds [0,1] and monotonicity

### ✅ Verification Tests

**Test Results:**
- All cause-specific models (Cox, GAM, XGBoost, SurvReg) pass tests
- Models produce different CIFs for different competing events
- CIFs are properly bounded and sum to ≤ 1 (proper calibration)
- Package loads successfully with all changes

## Mathematical Correctness

### ❌ **WRONG Approach** (Previously Used):
```r
# BIASED - inflates CIF curves
cif_matrix <- 1 - surv_probs_cause_specific
```

### ✅ **CORRECT Approach** (Now Implemented):
```r
# UNBIASED - proper Aalen-Johansen estimator
cif_matrix <- aalenJohansenCIF(
  cause_specific_survs = all_cause_survival_functions,
  times = times,
  event_of_interest = failcode
)
```

## Impact of Correction

### Before (Biased):
- CIF curves were inflated (too high)  
- Did not account for competing risks properly
- Sum of all CIFs could exceed 1

### After (Unbiased):
- CIF curves are properly calibrated
- Accounts for all competing events
- Sum of all CIFs ≤ 1 (mathematically correct)
- More accurate predictions when competing risks are significant

## Conclusion

✅ **ALL MODELS CORRECTLY IMPLEMENTED**: Every competing risks model in the package now uses mathematically correct methods for CIF calculation:

- **Cause-specific models** (Cox, GAM, XGBoost, SurvReg): Use Aalen-Johansen estimator
- **Subdistribution models** (Fine-Gray, DeepSurv): Use subdistribution hazards  
- **Package-based models** (BART, RF): Rely on package implementations
- **Hybrid models** (RuleFit): Use Fine-Gray internally

The package now provides properly calibrated CIF predictions that avoid the bias inherent in the `1 - S(t)` approximation.