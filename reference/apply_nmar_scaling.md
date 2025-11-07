# Apply scaling to a matrix using a recipe

Apply scaling to a matrix using a recipe

## Usage

``` r
apply_nmar_scaling(matrix_to_scale, recipe)
```

## Arguments

- matrix_to_scale:

  a numeric matrix with column names present in \`recipe\`.

- recipe:

  an object of class \`nmar_scaling_recipe\`.

## Value

a matrix with each named column centered and scaled using the recipe
