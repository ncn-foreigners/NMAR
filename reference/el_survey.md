# Empirical likelihood for survey designs (NMAR)

Internal method dispatched by \`el()\` when \`data\` is a
\`survey.design\`.

## Usage

``` r
# S3 method for class 'survey.design'
el(
  data,
  formula,
  auxiliary_means = NULL,
  standardize = TRUE,
  strata_augmentation = TRUE,
  trim_cap = Inf,
  control = list(),
  on_failure = c("return", "error"),
  variance_method = c("bootstrap", "none"),
  bootstrap_reps = 500,
  n_total = NULL,
  start = NULL,
  trace_level = 0,
  family = logit_family(),
  ...
)
```

## Arguments

- data:

  A \`survey.design\` created with \[survey::svydesign()\].

- formula:

  Two-sided formula with an NA-valued outcome on the LHS; auxiliaries on
  the first RHS and, optionally, missingness predictors on the second
  RHS partition.

- auxiliary_means:

  Named numeric vector of population means for auxiliary design columns.
  Names must match the materialized \`model.matrix\` columns on the
  first RHS (after formula expansion), including factor indicators and
  transformed terms. The intercept is always excluded.

- standardize:

  Logical; standardize predictors.

- strata_augmentation:

  Logical; when `TRUE` (default), augment the auxiliary design with
  stratum indicators and stratum shares when a strata structure is
  present in the survey design.

- trim_cap:

  Numeric; cap for EL weights (Inf = no trimming).

- control:

  List; solver control for \`nleqslv(control=...)\`.

- on_failure:

  Character; "return" or "error" on solver failure.

- variance_method:

  Character; "bootstrap" or "none".

- bootstrap_reps:

  Integer; reps when \`variance_method = "bootstrap"\`.

- n_total:

  Optional analysis-scale population size `N_pop`; required for
  respondents-only designs.

- start:

  Optional list of starting values passed to solver helpers.

- trace_level:

  Integer 0-3 controlling estimator logging detail.

- family:

  Missingness (response) model family specification (defaults to logit).

- ...:

  Passed to solver.

## Value

\`c('nmar_result_el','nmar_result')\`.

## Details

Implements the empirical likelihood estimator with design weights. If
`n_total` is supplied, it is treated as the analysis-scale population
size `N_pop` used in the design-weighted QLS system. If `n_total` is not
supplied, `sum(weights(design))` is used as `N_pop`. Design weights are
not rescaled internally; the EL equations use respondent weights and
`N_pop` via `T0 = N_pop - sum(d_i)` in the linkage equation. When
respondents-only designs are used (no NA in the outcome), `n_total` must
be provided; if auxiliaries are requested you must also provide
population auxiliary means via `auxiliary_means`. Result weights are the
unnormalized EL masses `d_i/D_i(theta)` on this analysis scale;
`weights(result, scale = "population")` sums to `N_pop`.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.
