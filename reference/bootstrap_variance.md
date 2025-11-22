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

  a \`data.frame\` or a \`survey.design\`.

- estimator_func:

  function that returns an S3 result object; the primary estimate is
  extracted via \`\$y_hat\` and convergence via \`\$converged\`.

- point_estimate:

  numeric; point estimate used for some survey variance formulas.

- ...:

  passed through to \`estimator_func\`.

## Details

\- For \`data.frame\` inputs, performs IID bootstrap by resampling rows
and rerunning \`estimator_func\` on each resample, then computing the
empirical variance of the replicate estimates. - For \`survey.design\`
inputs, converts the design to a bootstrap replicate-weight design with
\`svrep::as_bootstrap_design()\`, reconstructs the original sampling
design for each replicate weight vector, and passes the resulting
replicate estimates and replicate scaling factors to
\`survey::svrVar()\`.

\`estimator_func\` is typically an engine-level estimator (for example
the EL engine) and is called with the same arguments used for the point
estimate, except that the \`data\` argument is replaced by the resampled
data (IID) or the replicate \`survey.design\` (survey).

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
      

The `future` framework (via `future.seed = TRUE` in
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html))
ensures that each bootstrap replicate uses an independent L'Ecuyer-CMRG
random number stream derived from this seed. This gives reproducible
results across supported future backends (sequential, multisession,
cluster, and so on).
