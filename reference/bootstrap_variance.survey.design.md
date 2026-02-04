# Bootstrap for survey designs

Bootstrap for survey designs

## Usage

``` r
# S3 method for class 'survey.design'
bootstrap_variance(
  data,
  estimator_func,
  point_estimate,
  bootstrap_reps = 500,
  survey_na_policy = c("strict", "omit"),
  ...
)
```

## Arguments

- data:

  A `survey.design`.

- estimator_func:

  Function returning an object with a numeric scalar component `y_hat`
  and an optional logical component `converged`.

- point_estimate:

  Numeric scalar; used for survey bootstrap variance (passed to
  [`survey::svrVar()`](https://rdrr.io/pkg/survey/man/svrVar.html) as
  `coef`).

- bootstrap_reps:

  integer; number of bootstrap replicates.

- survey_na_policy:

  Character string specifying how to handle replicates that fail to
  produce estimates. Options:

  `"strict"`

  :   (default) Any failed replicate causes an error. This is a
      conservative default that makes instability explicit.

  `"omit"`

  :   Failed replicates are omitted. The corresponding `rscales` are
      also omitted to maintain correct variance scaling. Use with
      caution: if failures are non-random, variance may be biased.

- ...:

  Additional arguments. Some are consumed by
  [`bootstrap_variance()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/bootstrap_variance.md)
  itself (for example `resample_guard` for IID bootstrap or
  `bootstrap_settings`/`bootstrap_options`/`bootstrap_type`/`bootstrap_mse`
  or survey bootstrap). Remaining arguments are forwarded to
  `estimator_func`.

## Value

A list with components `se`, `variance`, and `replicates`.

## Details

This path constructs a replicate-weight design using
[`svrep::as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html)
and evaluates the estimator on each set of bootstrap replicate analysis
weights. Replicate evaluation starts from a shallow template copy of the
input survey design (including its ids/strata/fpc structure) and injects
each replicate's analysis weights by updating the design's probability
slots (`prob`/`allprob`) so that `weights(design)` returns the desired
replicate weights. This avoids replaying or reconstructing a
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html) call and
therefore supports designs created via
[`subset()`](https://rdrr.io/r/base/subset.html) and
[`update()`](https://rdrr.io/r/stats/update.html). **NA policy:** By
default, survey bootstrap uses a strict NA policy: if any replicate
fails to produce a finite estimate, the entire bootstrap fails with an
error. Setting `survey_na_policy = "omit"` drops failed replicates and
proceeds with the remaining replicates.

## Limitations

**Calibrated/post-stratified designs:** Post-hoc adjustments applied via
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html),
[`survey::postStratify()`](https://rdrr.io/pkg/survey/man/postStratify.html),
or [`survey::rake()`](https://rdrr.io/pkg/survey/man/rake.html) are not
supported here and will cause the function to error. These adjustments
are not recomputed when replicate weights are injected, so the replicate
designs would not reflect the intended calibrated/post-stratified
analysis.
