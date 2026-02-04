# Core of the empirical likelihood estimator

Core of the empirical likelihood estimator

## Usage

``` r
el_estimator_core(
  missingness_design,
  aux_matrix,
  aux_means,
  respondent_weights,
  analysis_data,
  outcome_expr,
  N_pop,
  formula,
  standardize,
  trim_cap,
  control,
  on_failure,
  family = logit_family(),
  variance_method,
  bootstrap_reps,
  start = NULL,
  trace_level = 0,
  auxiliary_means = NULL
)
```

## Arguments

- missingness_design:

  Respondent-side missingness model design matrix (intercept +
  predictors).

- aux_matrix:

  Auxiliary design matrix on respondents (may have zero columns).

- aux_means:

  Named numeric vector of auxiliary population means (aligned to columns
  of `aux_matrix`).

- respondent_weights:

  Numeric vector of respondent weights aligned with `missingness_design`
  rows.

- analysis_data:

  Data object used for logging and variance.

- outcome_expr:

  Character string identifying the outcome expression displayed in
  outputs.

- N_pop:

  Population size on the analysis scale.

- formula:

  Original model formula used for estimation.

- standardize:

  Logical. Whether to standardize predictors during estimation.

- trim_cap:

  Numeric. Upper bound for empirical likelihood weight trimming.

- control:

  List of control parameters for the nonlinear equation solver.

- on_failure:

  Character. Action when solver fails.

- family:

  List. Link function specification.

- variance_method:

  Character. Variance estimation method.

- bootstrap_reps:

  Integer. Number of bootstrap replications.

- auxiliary_means:

  Named numeric vector of known population means supplied by the user.

## Value

List containing estimation results, diagnostics, and metadata.
