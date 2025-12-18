# Extract weights from an \`nmar_result\`

Return analysis weights stored in an \`nmar_result\` as either
probability-scale (summing to 1) or population-scale (summing to
\`sample\$n_total\`). The function normalizes stored masses and attaches
informative attributes.

## Usage

``` r
# S3 method for class 'nmar_result'
weights(object, scale = c("probability", "population"), ...)
```

## Arguments

- object:

  An \`nmar_result\` object.

- scale:

  One of \`"probability"\` (default) or \`"population"\`.

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of weights with length equal to the number of
respondents.
