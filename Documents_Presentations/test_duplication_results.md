# Test Data Duplication Investigation - Results

## Question
Is test data duplication necessary for `BART::crisk.pre.bart()` and competing risks predictions?

## Current Implementation (Before Fix)
```r
# Duplicated test data
test_duplicated <- rbind(x_test, x_test)

pre <- BART::crisk.pre.bart(
  time = modelout$times_train,
  delta = modelout$delta_train,
  x.train = BART::bartModelMatrix(modelout$x_train),
  x.test = BART::bartModelMatrix(test_duplicated),  # Duplicated
  x.train2 = BART::bartModelMatrix(modelout$x_train),
  x.test2 = BART::bartModelMatrix(test_duplicated), # Duplicated
  K = modelout$bart_model$K
)
```

## Test Results

### Test Setup
- Training data: 100 observations, 2 features
- Test data: 20 observations, 2 features
- K = 10 time points

### Test 1: WITHOUT Duplication
✅ **SUCCESS**
- `crisk.pre.bart()` completed without errors
- tx.test dimensions: 200 × 3 (20 obs × 10 times)
- tx.test2 dimensions: 200 × 3

### Test 2: WITH Duplication
✅ **SUCCESS**
- `crisk.pre.bart()` completed without errors
- tx.test dimensions: 400 × 3 (40 obs × 10 times)
- tx.test2 dimensions: 400 × 3

### Test 3: Prediction Comparison
- **Without duplication**: CIF length = 200 (N × K)
- **With duplication**: CIF length = 400 (2 × N × K)
- **Maximum difference**: 0 (predictions are IDENTICAL)

## Conclusion

✅ **TEST DATA DUPLICATION IS NOT NECESSARY**

The duplication only affects the output length (doubles it) but does not change the actual prediction values. The first N × K values are identical whether or not duplication is used.

## Implementation Changes

### Before (Lines 235-271 in cr_bart.R)
```r
# Note: BART crisk.pre.bart expects duplicated test data for some reason
# (this seems to be a BART package requirement)
test_duplicated <- rbind(x_test, x_test)

pre <- BART::crisk.pre.bart(
  x.test = BART::bartModelMatrix(test_duplicated),
  x.test2 = BART::bartModelMatrix(test_duplicated),
  ...
)

# Handle potential duplication in output
total_len <- length(pred$cif.test.mean)
expected_len <- N * K
if (total_len < expected_len) {
  stop(...)
} else if (total_len > expected_len) {
  duplication_factor <- total_len / expected_len
  if (!duplication_factor %in% c(2, 1)) {
    warning(...)
  }
}
cif_vector <- pred$cif.test.mean[seq_len(expected_len)]
cif_matrix <- matrix(cif_vector, nrow = N, ncol = K, byrow = TRUE)
```

### After (Simplified)
```r
pre <- BART::crisk.pre.bart(
  x.test = BART::bartModelMatrix(x_test),
  x.test2 = BART::bartModelMatrix(x_test),
  ...
)

# Simple validation and reshaping
total_len <- length(pred$cif.test.mean)
expected_len <- N * K
if (total_len != expected_len) {
  stop("Unexpected length of BART CIF predictions. Expected ", expected_len, " values, got ", total_len, ".")
}
cif_matrix <- matrix(pred$cif.test.mean, nrow = N, ncol = K, byrow = TRUE)
```

## Benefits of the Fix

1. **Simpler code**: Removed 10 lines of unnecessary duplication and handling logic
2. **Better performance**: No need to duplicate test data (saves memory and computation)
3. **Clearer intent**: Code now directly reflects what BART actually requires
4. **Easier maintenance**: Less complex logic to understand and maintain

## Testing Verification

All existing unit tests pass after the change:
```
✔ | F W  S  OK | Context
✔ |         37 | Testing cr_bart functions [1.3s]

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 37 ]
```

## Files Modified

1. **`R/cr_bart.R`**: Removed test data duplication (lines 235-271)
   - Removed `test_duplicated <- rbind(x_test, x_test)`
   - Updated `crisk.pre.bart()` calls to use `x_test` directly
   - Simplified prediction extraction logic

## Documentation Updates Needed

The following documents should be updated to reflect this finding:

1. **`BART_Implementation_Details.Rmd`**
   - Remove "Issue 2: Test Data Duplication" from Known Issues section
   - Update competing risks prediction section
   - Remove question about duplication necessity

2. **`BART_Implementation_Summary.Rmd`**
   - Remove duplication from implementation details
   - Remove from known issues list

## Date
December 11, 2025

## Test Script
`Documents_Presentations/test_crisk_duplication.R`
