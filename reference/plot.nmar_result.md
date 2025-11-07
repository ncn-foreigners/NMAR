# Base plotting for NMAR results

Quick base plots for weights, fitted probabilities, constraints or
diagnostics.

## Usage

``` r
# S3 method for class 'nmar_result'
plot(x, which = c("weights", "fitted", "constraints", "diagnostics"), ...)
```

## Arguments

- x:

  An object of class \`nmar_result\`.

- which:

  Which plot: one of \`"weights"\`, \`"fitted"\`, \`"constraints"\`,
  \`"diagnostics"\`.

- ...:

  Ignored.
