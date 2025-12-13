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

  One or more numeric matrices with column names.

- intercept_col:

  Name of an intercept column that should remain unscaled.

- weights:

  Optional nonnegative numeric vector used to compute weighted means and
  standard deviations.

- weight_mask:

  Optional logical mask or nonnegative numeric multipliers applied to
  `weights` before computing moments (useful for respondents-only
  scaling). If `weights` is `NULL`, `weight_mask` is treated as weights.

- tol_constant:

  Numeric tolerance below which columns are treated as constant and left
  unscaled.

- warn_on_constant:

  Logical; warn when a column is treated as constant.
