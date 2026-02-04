# Unscale coefficients and covariance

Unscale coefficients and covariance

## Usage

``` r
unscale_coefficients(scaled_coeffs, scaled_vcov, recipe)
```

## Arguments

- scaled_coeffs:

  named numeric vector of coefficients estimated on the scaled space.

- scaled_vcov:

  covariance matrix of `scaled_coeffs`.

- recipe:

  Scaling recipe of class `nmar_scaling_recipe`.

## Value

A list with components `coefficients` and `vcov`.
