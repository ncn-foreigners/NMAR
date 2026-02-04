# Input preprocessing

Parses the two-part Formula, constructs EL design matrices, injects the
respondent delta indicator, attaches weights and survey metadata, and
returns the pieces needed by the EL core.

## Usage

``` r
el_prepare_inputs(
  formula,
  data,
  weights = NULL,
  n_total = NULL,
  design_object = NULL
)
```

## Details

Enforeces the following format required by the rest of el code:

- LHS references exactly one outcome source variable in `data`; any
  transforms are applied via the formula environment and must be defined
  for all respondent rows.

- The outcome is never allowed to appear on RHS1 (auxiliaries) or RHS2
  (missingness predictors), either explicitly in the formula or
  implicitly via dot (`.`) expansion. The missingness model uses the
  evaluated LHS expression as a dedicated predictor column instead.

- RHS1 always yields an intercept-free auxiliary design matrix with k-1
  coding for factor auxiliaries, regardless of user `+0`/`-1` syntax or
  custom contrasts. Auxiliary columns are validated to be fully observed
  and non-constant among respondents.

- RHS2 always yields a missingness-design matrix for respondents that
  includes an intercept column and zero-variance predictors emit
  warnings. NA among respondents is rejected.

- `respondent_mask` is defined from the raw outcome in `data`, not from
  the transformed LHS. An injected `..nmar_delta..` indicator in
  `analysis_data` must match this mask.

- `N_pop` is the analysis-scale population size: for IID it is
  `nrow(data)` unless overridden by `n_total`. For survey designs it is
  `sum(weights)` or `n_total` when supplied.
