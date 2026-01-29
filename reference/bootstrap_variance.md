# Shared bootstrap variance helpers

Internal helpers to estimate the variance of a scalar estimator via
bootstrap resampling (IID data) or bootstrap replicate weights (survey
designs). Designed to be reused across NMAR engines.

## Usage

``` r
bootstrap_variance(data, estimator_func, point_estimate, ...)
```

## Arguments

- data:

  A `data.frame` or a `survey.design`.

- estimator_func:

  Function returning an object with a numeric scalar component `y_hat`
  and an optional logical component `converged`.

- point_estimate:

  Numeric scalar; used for survey bootstrap variance (passed to
  [`survey::svrVar()`](https://rdrr.io/pkg/survey/man/svrVar.html) as
  `coef`).

- ...:

  Additional arguments. Some are consumed by `bootstrap_variance()`
  itself (for example `resample_guard` for IID bootstrap or
  `bootstrap_settings`/`bootstrap_options`/`bootstrap_type`/`bootstrap_mse`
  for survey bootstrap); remaining arguments are forwarded to
  `estimator_func`.

## Details

- For `data.frame` inputs, performs IID bootstrap by resampling rows and
  rerunning `estimator_func` on each resample, then computing the
  empirical variance of the replicate estimates.

- For `survey.design` inputs, converts the design to a bootstrap
  replicate-weight design with
  [`svrep::as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html),
  evaluates `estimator_func` on each replicate weight vector (by
  injecting the replicate analysis weights into a copy of the input
  design), and passes the resulting replicate estimates and replicate
  scaling factors to
  [`survey::svrVar()`](https://rdrr.io/pkg/survey/man/svrVar.html).

`estimator_func` is typically an engine-level estimator (for example the
EL engine) and is called with the same arguments used for the point
estimate, except that the `data` argument is replaced by the resampled
data (IID) or a replicate-weighted `survey.design` (survey). Arguments
reserved for the bootstrap implementation are stripped from `...` before
forwarding.

## Bootstrap-specific options

- `resample_guard`:

  IID bootstrap only. A function `function(indices, data)` that returns
  `TRUE` to accept a resample and `FALSE` to reject it.

- `bootstrap_settings`:

  Survey bootstrap only. A list of arguments forwarded to
  [`svrep::as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html).

- `bootstrap_options`:

  Alias for `bootstrap_settings`.

- `bootstrap_type`:

  Shortcut for the `type` argument to
  [`svrep::as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html).

- `bootstrap_mse`:

  Shortcut for the `mse` argument to
  [`svrep::as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html).

## Progress Reporting

If the optional `progressr` package is installed, bootstrap calls signal
progress via a
[`progressr::progressor`](https://progressr.futureverse.org/reference/progressor.html)
inside
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html).
Users control whether progress is shown (and how) by registering
handlers with
[`progressr::handlers()`](https://progressr.futureverse.org/reference/handlers.html).
When `progressr` is not installed or no handlers are active, bootstrap
runs silently. Progress reporting is compatible with all future
backends.

## Reproducibility

For reproducible bootstrap results, always set a seed before calling the
estimation function:

      set.seed(123)  # Set seed for reproducibility
      result <- nmar(Y ~ X, data = df,
                     engine = el_engine(variance_method = "bootstrap",
                                        bootstrap_reps = 500))
      

## Parallelization

By default, bootstrap replicate evaluation runs sequentially via
[`base::lapply()`](https://rdrr.io/r/base/lapply.html) for both IID
resampling and survey replicate-weight bootstrap. If the optional
`future.apply` package is installed, bootstrap can use
`future.apply::future_lapply(future.seed = TRUE)` when the user has set
a parallel
[`future::plan()`](https://future.futureverse.org/reference/plan.html).

The backend is controlled by the package option `nmar.bootstrap_apply`:

- `"auto"`:

  (default) Use [`base::lapply()`](https://rdrr.io/r/base/lapply.html)
  unless the current future plan has more than one worker, in which case
  use
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
  if available.

- `"base"`:

  Always use [`base::lapply()`](https://rdrr.io/r/base/lapply.html)
  (never use `future.apply`, even if installed).

- `"future"`:

  Always use
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
  (requires `future.apply` to be installed).

When `future.apply` is used, random-number streams are parallel-safe and
backend-independent under the `future` framework. When
[`base::lapply()`](https://rdrr.io/r/base/lapply.html) is used, results
are reproducible under
[`set.seed()`](https://rdrr.io/r/base/Random.html) but will not match
the `future.seed` streams.
