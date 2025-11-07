# Empirical likelihood for survey designs (NMAR)

Internal method dispatched by \`el()\` when \`data\` is a
\`survey.design\`. Variance via bootstrap is supported. Analytical delta
variance for EL is not implemented and returns NA when requested.

## Usage

``` r
# S3 method for class 'survey.design'
el(
  data,
  formula,
  auxiliary_means = NULL,
  standardize = TRUE,
  trim_cap = Inf,
  control = list(),
  on_failure = c("return", "error"),
  variance_method = c("delta", "bootstrap", "none"),
  bootstrap_reps = 500,
  n_total = NULL,
  start = NULL,
  trace_level = 0,
  ...
)
```

## Arguments

- data:

  A \`survey.design\` created with \[survey::svydesign()\].

- formula:

  Two-sided formula: NA-valued outcome on LHS; auxiliaries on RHS.

- auxiliary_means:

  Named numeric vector of population means for auxiliaries.

- standardize:

  Logical; standardize predictors.

- trim_cap:

  Numeric; cap for EL weights (Inf = no trimming).

- control:

  List; solver control for \`nleqslv(control=...)\`.

- on_failure:

  Character; "return" or "error" on solver failure.

- variance_method:

  Character; "delta" or "bootstrap".

- bootstrap_reps:

  Integer; reps when \`variance_method = "bootstrap"\`.

- ...:

  Passed to solver.

## Value

\`c('nmar_result_el','nmar_result')\`.

## Details

Implements the empirical likelihood estimator with design weights. If
`n_total` is supplied, design weights are rescaled internally to ensure
`sum(weights(design))` and `n_total` are on the same scale; this
guarantees the response-multiplier formula uses consistent totals. If
`n_total` is not supplied, `sum(weights(design))` is used as the
population total `N_pop`. When respondents-only designs are used (no NA
in the outcome), `n_total` must be provided; if auxiliaries are
requested you must also provide population auxiliary means via
`auxiliary_means`. Result weights are the unnormalized EL masses
`d_i/D_i(theta)` on this design
scale;`weights(result, scale = "population")` sums to `N_pop`.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.

Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using
complete auxiliary information from survey data. Journal of the American
Statistical Association, 96(453), 185-193.
