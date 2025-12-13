# Apply scaling to a matrix using a recipe

Apply scaling to a matrix using a recipe

## Usage

``` r
apply_nmar_scaling(matrix_to_scale, recipe)
```

## Arguments

- matrix_to_scale:

  A numeric matrix with column names present in `recipe`.

- recipe:

  An object of class `nmar_scaling_recipe`.

## Value

A numeric matrix with each column centered and scaled using `recipe`.
