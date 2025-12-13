# Prepare EL inputs for IID and survey designs

Parses the two-part Formula, constructs EL design matrices, injects the
respondent delta indicator, attaches weights and (optionally) survey
metadata, and returns the pieces needed by the EL core. The outcome
enters the missingness design only through the evaluated LHS expression;
any explicit use of the outcome variable on RHS2 is rejected.

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

Invariants enforced here (relied on by all downstream EL code):

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
  includes an intercept column and zero-variance predictors only emit
  warnings (not errors); NA among respondents is rejected.

- `respondent_mask` is defined from the raw outcome in `data`, not from
  the transformed LHS; an injected `..nmar_delta..` indicator in
  `analysis_data` matches this mask exactly.

- `N_pop` is the analysis-scale population size used in the EL system:
  for IID it is `nrow(data)` unless overridden by `n_total`; for survey
  designs it is `sum(weights)` or `n_total` when supplied.
