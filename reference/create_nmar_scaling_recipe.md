# Build a scaling recipe from one or more design matrices

Build a scaling recipe from one or more design matrices

## Usage

``` r
create_nmar_scaling_recipe(
  ...,
  intercept_col = "(Intercept)",
  weights = NULL,
  weight_mask = NULL,
  tol_constant = 1e-08,
  warn_on_constant = TRUE
)
```

## Arguments

- ...:

  one or more matrices with named columns.

- intercept_col:

  Intercept column name that should remain unscaled.

- weights:

  Optional numeric vector of weights used to compute weighted
  means/standard deviations.

- weight_mask:

  Optional logical/ numeric mask applied to \`weights\` before computing
  moments (useful for respondents-only scaling).

- tol_constant:

  Numeric tolerance below which columns are treated as constant and left
  unscaled.

- warn_on_constant:

  Logical; emit a warning when a column is treated as constant.
