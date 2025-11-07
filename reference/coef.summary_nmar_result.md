# Coefficient table for summary objects

Returns a coefficients table (Estimate, Std. Error, statistic, p-value)
from a \`summary_nmar_result\*\` object when response-model coefficients
and a variance matrix are available. If the summary does not carry
response-model coefficients, returns \`NULL\`.

## Usage

``` r
# S3 method for class 'summary_nmar_result'
coef(object, ...)
```

## Arguments

- object:

  An object of class \`summary_nmar_result\` (or subclass).

- ...:

  Ignored.

## Value

A data.frame with rows named by coefficient, or \`NULL\` if not
available.

## Details

The statistic column is labelled "t value" when finite degrees of
freedom are available (e.g., survey designs); otherwise, it is labelled
"z value".
