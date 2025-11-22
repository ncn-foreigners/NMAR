# Shared scaling infrastructure for NMAR engines

Centralized feature scaling and parameter unscaling routines used by
NMAR estimation engines to ensure consistent, numerically stable
behavior.

## Goals

- Provide an engine-agnostic API for standardizing design matrices and
  auxiliary moments before solving.

- Return a minimal "recipe" (per-column mean and standard deviation) for
  unscaling coefficients and covariance matrices after solving.

## Inputs/Outputs

- Inputs:

  `Z_un` (response model matrix with intercept), optional `X_un`
  (auxiliary model matrix, no intercept), optional named `mu_x_un`
  (auxiliary means on the original scale), and a logical `standardize`
  flag.

- Outputs:

  Scaled matrices `Z`, `X`, and `mu_x`, plus an `nmar_scaling_recipe`
  used later for unscaling.

## Integration pattern

1.  Before solving: call
    [`validate_and_apply_nmar_scaling()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/validate_and_apply_nmar_scaling.md)
    (engine-level) or
    [`prepare_nmar_scaling()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/prepare_nmar_scaling.md)
    (low-level) to obtain scaled matrices and recipe.

2.  Solve in the scaled space.

3.  After solving: call
    [`unscale_coefficients()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/unscale_coefficients.md)
    to unscale coefficients and their covariance matrices.

4.  Store the `nmar_scaling_recipe` in results for diagnostics and
    reproducibility.

## Notes

- The intercept column is never scaled.

- Columns with near-zero variance are centered but assigned `sd = 1` so
  that the corresponding parameter is not inflated by division by a very
  small standard deviation.

- Engines may use design-weighted scaling via the `weights` and
  `weight_mask` arguments.
