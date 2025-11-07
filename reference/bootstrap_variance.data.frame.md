# Bootstrap for i.i.d. data.frames

Bootstrap for i.i.d. data.frames

## Usage

``` r
bootstrap_variance.data.frame(
  data,
  estimator_func,
  point_estimate,
  bootstrap_reps = 500,
  ...
)
```

## Arguments

- data:

  a \`data.frame\` or a \`survey.design\`.

- estimator_func:

  function that returns an S3 result object; the primary estimate is
  extracted via \`\$y_hat\` and convergence via \`\$converged\`.

- point_estimate:

  numeric; point estimate used for some survey variance formulas.

- bootstrap_reps:

  integer; number of resamples.

- ...:

  passed through to \`estimator_func\`.

## Value

a list with \`se\`, \`variance\`, and the vector of \`replicates\`.
