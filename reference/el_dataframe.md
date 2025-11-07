# Empirical likelihood for data frames (NMAR)

Internal method dispatched by \`el()\` when \`data\` is a
\`data.frame\`. Returns \`c('nmar_result_el','nmar_result')\` with the
point estimate, optional bootstrap SE, weights, coefficients,
diagnostics, and metadata.

## Usage

``` r
# S3 method for class 'data.frame'
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

  A \`data.frame\` where the outcome column contains \`NA\` for
  nonrespondents.

- formula:

  Two-sided formula \`Y_miss ~ auxiliaries\`.

- auxiliary_means:

  Named numeric vector of population means for auxiliary variables
  (names must match RHS of outcome formula).

- standardize:

  Logical; whether to standardize predictors prior to estimation.

- trim_cap:

  Numeric; cap for EL weights (\`Inf\` = no trimming).

- control:

  List; optional solver control parameters for \`nleqslv(control=...)\`.

- on_failure:

  Character; one of \`"return"\` or \`"error"\` on solver failure.

- variance_method:

  Character; one of \`"delta"\`, \`"bootstrap"\`, or \`"none"\`.

- bootstrap_reps:

  Integer; number of bootstrap reps if \`variance_method =
  "bootstrap"\`.

- ...:

  Additional arguments passed to the solver.

## Details

Implements the empirical likelihood estimator for IID data with optional
auxiliary moment constraints. The response-model score is the Bernoulli
derivative with respect to the linear predictor, supporting logit and
probit links. When respondents-only data is supplied (no NA in the
outcome), set `n_total` to the total number of sampled units; otherwise
the total is taken as `nrow(data)`. If respondents-only data is used and
auxiliaries are requested, you must also provide population auxiliary
means via `auxiliary_means`. Result weights are the unnormalized EL
masses `d_i/D_i(theta)` on this analysis scale.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.

Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using
complete auxiliary information from survey data. Journal of the American
Statistical Association, 96(453), 185-193.
