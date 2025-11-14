# Confidence intervals for coefficient table (summary objects)

Returns Wald-style confidence intervals for missingness-model
coefficients from a \`summary_nmar_result\*\` object. Uses t-quantiles
when finite degrees of freedom are available, otherwise normal
quantiles.

## Usage

``` r
# S3 method for class 'summary_nmar_result'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class \`summary_nmar_result\` (or subclass).

- parm:

  A specification of which coefficients are to be given confidence
  intervals, either a vector of names or a vector of indices; by
  default, all coefficients are considered.

- level:

  The confidence level required.

- ...:

  Ignored.

## Value

A numeric matrix with columns giving lower and upper confidence limits
for each parameter. Row names correspond to coefficient names. Returns
\`NULL\` if coefficients are unavailable.
