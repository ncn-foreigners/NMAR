# Validate and apply scaling for engines

Validate and apply scaling for engines

## Usage

``` r
validate_and_apply_nmar_scaling(
  standardize,
  has_aux,
  response_model_matrix_unscaled,
  aux_matrix_unscaled,
  mu_x_unscaled,
  weights = NULL,
  weight_mask = NULL
)
```

## Arguments

- standardize:

  logical; apply standardization if TRUE.

- has_aux:

  logical; whether the engine uses auxiliary constraints.

- response_model_matrix_unscaled:

  response model matrix (with intercept).

- aux_matrix_unscaled:

  auxiliary matrix (no intercept) or an empty matrix.

- mu_x_unscaled:

  named auxiliary means on original scale, or NULL.

- weights:

  Optional numeric vector used for weighted scaling.

- weight_mask:

  Optional logical mask or nonnegative numeric multipliers applied to
  `weights`.

## Value

A list with components `nmar_scaling_recipe`,
`response_model_matrix_scaled`, `auxiliary_matrix_scaled`, and
`mu_x_scaled`.
