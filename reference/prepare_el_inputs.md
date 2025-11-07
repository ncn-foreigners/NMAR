# Prepare inputs for EL estimation

Prepare inputs for EL estimation

## Usage

``` r
prepare_el_inputs(formula, data, require_na = TRUE)
```

## Details

Validates the two-sided outcome formula and constructs three internal
formulas: outcome (~ outcome_var), response (for the missingness model),
and auxiliary (RHS only, no intercept). Response-only predictors may
include variables not on the outcome RHS; such variables enter only the
response model (no auxiliary moment constraint). Only variables on the
outcome RHS are treated as auxiliaries and, when provided, must match
the names in \`auxiliary_means\`. See Qin, Leung and Shao (2002) for the
EL formulation.
