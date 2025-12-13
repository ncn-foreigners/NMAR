# Map unscaled auxiliary multipliers to scaled space

Map unscaled auxiliary multipliers to scaled space

## Usage

``` r
scale_aux_multipliers(lambda_unscaled, recipe, columns)
```

## Arguments

- lambda_unscaled:

  named numeric vector of auxiliary multipliers aligned to auxiliary
  design columns (no intercept) on original scale.

- recipe:

  Scaling recipe of class `nmar_scaling_recipe`.

- columns:

  character vector of auxiliary column names (order) for the scaled
  design.

## Value

numeric vector of multipliers in the scaled space.
