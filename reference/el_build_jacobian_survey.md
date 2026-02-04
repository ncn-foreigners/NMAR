# Empirical likelihood analytical jacobian for survey designs

Empirical likelihood analytical jacobian for survey designs

## Usage

``` r
el_build_jacobian_survey(
  family,
  missingness_model_matrix,
  auxiliary_matrix,
  respondent_weights,
  N_pop,
  n_resp_weighted,
  mu_x_scaled
)
```

## Details

Builds the block Jacobian \\A = \partial g/\partial \theta\\ for the
survey EL system with \\\theta = (\beta, z, \lambda_W, \lambda_x)\\ and
\\z = \operatorname{logit}(W)\\.
