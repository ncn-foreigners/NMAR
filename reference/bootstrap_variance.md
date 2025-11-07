# Shared bootstrap variance helpers

S3 generic + methods to estimate the variance of an estimator via
resampling (IID) or replicate weights (survey). Designed to be reused
across NMAR engines.

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

\- For \`data.frame\` inputs, performs i.i.d. bootstrap by resampling
rows and rerunning \`estimator_func\`. - For \`survey.design\` inputs,
converts to a bootstrap replicate-weight design
(\`svrep::as_bootstrap_design\`), computes replicate estimates by
rebuilding the original design with each replicate weight vector, and
then computes variance with \`survey::svrVar\` using the replicate
scales. \`estimator_func\` is typically an engine method (e.g.,
\`el()\`), and is called with the same arguments used for the point
estimate, except that the \`data\` argument is replaced by the resampled
data or replicate design.

## Progress Reporting

If the `progressr` package is installed, progress reporting is
available. Enable it by setting handlers before calling the bootstrap:

[`library(progressr)`](https://progressr.futureverse.org)

`handlers(global = TRUE)`

`handlers("txtprogressbar") # or "progress", "cli", etc.`

To disable progress in simulations or batch jobs, use
`handlers("void")`. If progressr is not installed or no handlers are
set, bootstrap runs silently (default behavior). Progress reporting
works with all future backends (sequential, multisession, cluster, etc.)
and does not affect reproducibility.

## Reproducibility

For reproducible bootstrap results, always set a seed before calling the
estimation function:

      set.seed(123)  # Set seed for reproducibility
      result <- nmar(Y ~ X, data = df,
                     engine = el_engine(variance_method = "bootstrap",
                                        bootstrap_reps = 500))
      

The `future` package (via `future.seed = TRUE`) ensures each bootstrap
replicate uses an independent L'Ecuyer-CMRG random number stream derived
from this seed, guaranteeing reproducibility across all future backends
(sequential, multisession, cluster, etc.).
