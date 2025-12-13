# Map unscaled coefficients to scaled space

Map unscaled coefficients to scaled space

## Usage

``` r
scale_coefficients(beta_unscaled, recipe, columns)
```

## Arguments

- beta_unscaled:

  named numeric vector of coefficients for the response model on the
  original scale, including an intercept named `"(Intercept)"`.

- recipe:

  Scaling recipe of class `nmar_scaling_recipe`, or `NULL`.

- columns:

  character vector of column names (order) for the scaled design matrix
  (including intercept).

## Value

numeric vector of coefficients in the scaled space, ordered by
`columns`.
