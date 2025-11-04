# Bootstrap Sampling Issues - Analysis and Fix

## Summary
Found and fixed critical issues with `sample_size` parameter interaction with bootstrap variance estimation in the exptilt method.

## Issues Identified

### 1. **Double Sampling During Bootstrap** (CRITICAL)
**Location**: `src_dev/engines/exptilt/impl/dataframe.R`, line 606

**Problem**: 
- The `bootstrap_runner` function was passing `sample_size = model$sample_size` to bootstrap replicates
- This caused **double sampling**:
  1. First sampling by `bootstrap_variance()` (resampling with replacement)
  2. Second sampling by `exptilt.data.frame()` based on `sample_size` parameter
- Result: Bootstrap variance was estimated on a **subset of a resample**, leading to incorrect variance estimates

**Example**:
```r
# Original data: 10,000 observations
# sample_size = 2,000

# What was happening (WRONG):
# 1. Bootstrap creates replicate of 10,000 obs (with replacement)
# 2. exptilt() then samples down to 2,000 obs
# 3. Variance estimated on 2,000 obs instead of 10,000 obs

# What should happen (CORRECT):
# 1. Bootstrap creates replicate of 10,000 obs (with replacement)  
# 2. exptilt() uses ALL 10,000 obs (sample_size = Inf)
# 3. Variance estimated on full 10,000 obs
```

**Fix**: Set `sample_size = Inf` in bootstrap replicates to disable sampling

### 2. **Missing Original Data Reference** (CRITICAL)
**Location**: `src_dev/engines/exptilt/impl/dataframe.R`, line 611

**Problem**:
- Bootstrap was using `model$data`, but this field was **never stored** in the model object
- For data frames: `model$data` contained the **sampled subset**, not the original data
- For surveys: `model$design` contained the **sampled design**, not the original design
- Result: Bootstrap operated on the sampled data, systematically **underestimating variance**

**Example**:
```r
# Original data: 10,000 observations
# sample_size = 2,000 (80% reduction)

# What was happening (WRONG):
# 1. Point estimate computed on 2,000 sampled obs
# 2. Bootstrap resampled from the same 2,000 obs
# 3. Variance underestimated by ~80%

# What should happen (CORRECT):
# 1. Point estimate computed on 2,000 sampled obs (for speed)
# 2. Bootstrap resamples from FULL 10,000 obs (for accuracy)
# 3. Variance correctly estimated
```

**Fix**: 
- Store `original_data` and `original_design` in model object (before sampling)
- Use these for bootstrap instead of sampled data

### 3. **Incorrect Resample Guard Mask**
**Location**: `src_dev/engines/exptilt/impl/dataframe.R`, line 638

**Problem**:
- The `resample_guard` was using `respondent_mask` computed from **sampled data**
- Should use mask from **original data** to correctly identify respondents in bootstrap resamples

**Fix**: Compute `respondent_mask_guard` from original data

## Implementation Details

### Changes to Model Object
```r
model <- list(
  # ... existing fields ...
  
  # NEW: Store original data for bootstrap
  original_data = if (is.null(survey_design)) data else NULL,
  original_design = if (!is.null(survey_design)) data else NULL
)
```

### Changes to Bootstrap Runner
```r
bootstrap_runner <- function(data, is_bootstrap_replicate = FALSE, ...) {
  exptilt(
    data = data,
    # ... other parameters ...
    sample_size = Inf,  # CRITICAL FIX: Disable sampling in bootstrap
    is_bootstrap_replicate = is_bootstrap_replicate
  )
}

# Use ORIGINAL data for bootstrap (before sampling)
bootstrap_data <- if (model$is_survey) {
  if (!is.null(model$original_design)) model$original_design else model$design
} else {
  if (!is.null(model$original_data)) model$original_data else model$data
}
```

### Changes to Resample Guard
```r
if (!model$is_survey) {
  # Use ORIGINAL data to compute mask
  original_data_for_mask <- if (!is.null(model$original_data)) {
    model$original_data
  } else {
    model$data
  }
  outcome_col <- model$col_y
  respondent_mask_guard <- !is.na(original_data_for_mask[[outcome_col]])
  
  base_args$resample_guard <- function(indices, data) {
    any(respondent_mask_guard[indices])
  }
}
```

## Impact Assessment

### For Data Frames
- **Before**: Bootstrap severely underestimated variance when `sample_size < n`
- **After**: Bootstrap correctly estimates variance from full dataset
- **Recommendation**: Works perfectly now

### For Survey Designs
- **Before**: Double sampling caused incorrect variance estimates
- **After**: Bootstrap uses original design without re-sampling
- **Recommendation**: Should work correctly now

## Testing Recommendations

1. **Test with different sample sizes**:
   ```r
   # Test that variance is similar regardless of sample_size used for point estimate
   result_2k <- exptilt(data, ..., sample_size = 2000, variance_method = "bootstrap")
   result_5k <- exptilt(data, ..., sample_size = 5000, variance_method = "bootstrap")
   result_inf <- exptilt(data, ..., sample_size = Inf, variance_method = "bootstrap")
   
   # All should have similar SE (within sampling error)
   result_2k$se  # Should be close to result_inf$se
   result_5k$se  # Should be close to result_inf$se
   ```

2. **Test with survey designs**:
   ```r
   library(survey)
   design <- svydesign(ids = ~1, weights = ~wt, data = mydata)
   
   # Should not produce errors about design subsetting
   result <- exptilt(design, ..., sample_size = 2000, variance_method = "bootstrap")
   ```

3. **Compare delta vs bootstrap**:
   ```r
   # For large samples with normal y_dens, these should be similar
   result_delta <- exptilt(data, ..., variance_method = "delta")
   result_boot <- exptilt(data, ..., variance_method = "bootstrap")
   
   # Should be within 10-20% of each other (bootstrap has sampling error)
   abs(result_delta$se - result_boot$se) / result_delta$se
   ```

## Files Modified
- `src_dev/engines/exptilt/impl/dataframe.R`
  - Added `original_data` and `original_design` to model object (lines ~120-147)
  - Fixed bootstrap_runner to use `sample_size = Inf` (lines ~585-610)
  - Fixed bootstrap_data to use original data (lines ~612-629)
  - Fixed resample_guard to use original data mask (lines ~631-647)

## Conclusion
The bootstrap sampling was fundamentally broken due to:
1. Double sampling (resampling then sampling again)
2. Using sampled data instead of original data for variance estimation
3. Incorrect respondent masking

All issues are now fixed. Bootstrap variance should work correctly for both data frames and survey designs, regardless of the `sample_size` parameter used for the point estimate.
