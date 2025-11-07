# Prepare scaled matrices and moments (low-level)

Prepare scaled matrices and moments (low-level)

## Usage

``` r
prepare_nmar_scaling(
  Z_un,
  X_un,
  mu_x_un,
  standardize,
  weights = NULL,
  weight_mask = NULL
)
```

## Arguments

- Z_un:

  response model matrix (with intercept column).

- X_un:

  auxiliary model matrix (no intercept), or NULL.

- mu_x_un:

  named numeric vector of auxiliary means on the original scale (names
  must match \`colnames(X_un)\`), or NULL.

- standardize:

  logical; apply standardization if TRUE.

- weights:

  Optional numeric vector used for weighted scaling.

- weight_mask:

  Optional logical/numeric mask applied to \`weights\`.

## Value

a list with \`Z\`, \`X\`, \`mu_x\`, and \`recipe\`.
