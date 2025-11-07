# Build denominator and floor pack

Build denominator and floor pack

## Usage

``` r
el_denominator(lambda_W, W, Xc_lambda, p_i, floor)
```

## Arguments

- lambda_W:

  numeric scalar

- W:

  numeric scalar in (0,1)

- Xc_lambda:

  numeric vector (X_centered %\*% lambda_x) or 0

- p_i:

  numeric vector of response probabilities

- floor:

  numeric scalar \> 0, denominator floor

## Value

list with denom, active, inv, inv_sq
