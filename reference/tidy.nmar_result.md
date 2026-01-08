# Tidy summary for NMAR results

Return a data frame with the primary estimate and (if available)
missingness-model coefficients.

## Usage

``` r
# S3 method for class 'nmar_result'
tidy(x, conf.level = 0.95, ...)
```

## Arguments

- x:

  An object of class \`nmar_result\`.

- conf.level:

  Confidence level for the primary estimate.

- ...:

  Ignored.

## Value

A data frame with one row for the primary estimate and, when available,
additional rows for the response-model coefficients.
