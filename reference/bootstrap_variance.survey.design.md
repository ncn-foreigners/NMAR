# Bootstrap for survey designs via replicate weights

Bootstrap for survey designs via replicate weights

## Usage

``` r
bootstrap_variance.survey.design(
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

  a \`data.frame\` or a \`survey.design\`.

- estimator_func:

  function that returns an S3 result object; the primary estimate is
  extracted via \`\$y_hat\` and convergence via \`\$converged\`.

- point_estimate:

  numeric; point estimate used for some survey variance formulas.

- bootstrap_reps:

  integer; number of bootstrap replicates.

- survey_na_policy:

  Character string specifying how to handle replicates that fail to
  produce estimates. Options:

  `"strict"`

  :   (default) Any failed replicate causes an error. This ensures the
      full replicate design structure is maintained and is required for
      proper calibration-based variance.

  `"omit"`

  :   Failed replicates are omitted. The corresponding `rscales` are
      also omitted to maintain correct variance scaling. Use with
      caution: if failures are non-random, variance may be biased.

- ...:

  passed through to \`estimator_func\`.

## Value

a list with \`se\`, \`variance\`, and the vector of \`replicates\`.

## Details

This path constructs a replicate-weight design using
\[svrep::as_bootstrap_design()\] and rebuilds the original sampling
design for each replicate weight vector. The supplied design must have
been created directly with \[survey::svydesign()\].

**NA Policy:** Survey bootstrap uses a strict NA policy - if any
replicate fails to produce a finite estimate, the entire bootstrap fails
with an error. This ensures the replicate design structure is maintained
for design-calibrated variance via
[`survey::svrVar()`](https://rdrr.io/pkg/survey/man/svrVar.html). In
contrast, IID bootstrap allows up to 10% failures before warning, as it
uses uncalibrated [`stats::var()`](https://rdrr.io/r/stats/cor.html).

## Limitations

**Design Reconstruction:** Survey bootstrap currently supports only
designs created directly with
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).
Post-hoc adjustments applied via
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html),
[`survey::postStratify()`](https://rdrr.io/pkg/survey/man/postStratify.html),
or [`survey::rake()`](https://rdrr.io/pkg/survey/man/rake.html) cannot
be reconstructed across bootstrap replicates and will cause the function
to error.

Calibrated or post-stratified designs are not supported by this
bootstrap path. Start from the original \`survey::svydesign()\` object
prior to calibration/post-stratification.

**Supported Design Features:** The following `svydesign()` parameters
are preserved during reconstruction:

- `ids` (or `id`): Sampling unit identifiers

- `strata`: Stratification variables

- `fpc`: Finite population correction

- `nest`: Nested vs non-nested strata

The following are NOT preserved (they conflict with replicate weights):

- `probs`: Sampling probabilities (incompatible with direct weights)

- `pps`: PPS sampling specification (incompatible with direct weights)

**Rationale:** When reconstructing designs for each replicate, we
replace the original weights with bootstrap replicate weights.
Specifying both `weights` and `probs`/`pps` simultaneously is undefined
behavior in
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).
The structural design parameters (ids, strata, fpc, nest) define the
sampling topology and are preserved; the weights define the analysis
weights and are replaced.
