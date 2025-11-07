# Numeric gradient helper

Thin wrapper around numDeriv::grad for internal use across engines.

## Usage

``` r
grad_numeric(x, func)
```

## Arguments

- x:

  Numeric vector at which to evaluate the gradient.

- func:

  Function mapping numeric vector to scalar; receives \`x\`.
