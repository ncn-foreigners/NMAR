# Unscale regression coefficients and covariance

Unscale regression coefficients and covariance

## Usage

``` r
unscale_coefficients(scaled_coeffs, scaled_vcov, recipe)
```

## Arguments

- scaled_coeffs:

  named numeric vector of coefficients estimated on the scaled space.

- scaled_vcov:

  covariance matrix of \`scaled_coeffs\`.

- recipe:

  \`nmar_scaling_recipe\` produced when scaling was applied.

## Value

a list with unscaled \`coefficients\` and \`vcov\`.
