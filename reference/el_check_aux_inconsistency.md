# Check auxiliary means consistency against respondents' sample support

Check auxiliary means consistency against respondents' sample support

## Usage

``` r
el_check_aux_inconsistency(
  respondent_df,
  aux_formula,
  provided_means = NULL,
  threshold = 8
)
```

## Arguments

- respondent_df:

  data.frame of respondents (no NAs in outcome indicator)

- aux_formula:

  RHS-only formula for auxiliaries (no intercept)

- provided_means:

  optional named numeric vector of auxiliary means on the same columns
  as model matrix

- threshold:

  numeric, z-score threshold for flagging

## Value

list(max_z = numeric(1) or NA, cols = character())
