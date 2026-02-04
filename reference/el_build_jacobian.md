# Empirical likelihood analytical jacobian for srs

Builds the block Jacobian \\A = \partial F/\partial \theta\\ for the EL
system with \\\theta = (\beta, z, \lambda_x)\\ and \\z =
\operatorname{logit}(W)\\. Blocks follow QLS equations 7-10.

## Usage

``` r
el_build_jacobian(
  family,
  missingness_model_matrix,
  auxiliary_matrix,
  respondent_weights,
  N_pop,
  n_resp_weighted,
  mu_x_scaled
)
```
