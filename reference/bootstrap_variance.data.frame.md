# Bootstrap for IID data frames

Bootstrap for IID data frames

## Usage

``` r
# S3 method for class 'data.frame'
bootstrap_variance(
  data,
  estimator_func,
  point_estimate,
  bootstrap_reps = 500,
  ...
)
```

## Arguments

- data:

  A `data.frame`.

- estimator_func:

  Function returning an object with a numeric scalar component `y_hat`
  and an optional logical component `converged`.

- point_estimate:

  Unused for IID bootstrap, included for signature consistency.

- bootstrap_reps:

  integer; number of resamples.

- ...:

  Additional arguments. Some are consumed by
  [`bootstrap_variance()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/bootstrap_variance.md)
  itself (for example `resample_guard` for IID bootstrap or
  `bootstrap_settings`/`bootstrap_options`/`bootstrap_type`/`bootstrap_mse`
  or survey bootstrap). Remaining arguments are forwarded to
  `estimator_func`.

## Value

A list with components `se`, `variance`, and `replicates`.
