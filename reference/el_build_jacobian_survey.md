# Analytical Jacobian for survey EL system (design-weighted QLS analogue)

Analytical Jacobian for survey EL system (design-weighted QLS analogue)

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
\\z = \operatorname{logit}(W)\\. Blocks follow the design-weighted
analogue of Qin, Leung, and Shao (2002) used in
[`el_build_equation_system_survey()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_equation_system_survey.md).
Guarding policy matches the IID Jacobian:

- cap eta: `eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())`

- compute `w <- family$linkinv(eta)` and clip to `[1e-12, 1-1e-12]` when
  used in ratios

- denominator floor: `Di <- pmax(Di_raw, nmar_get_el_denom_floor())`;
  multiply terms depending on `d(1/Di)/d(.)` by
  `active = 1(Di_raw > floor)`

The Jacobian uses the same score and second-derivative machinery as
[`el_build_jacobian()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_jacobian.md);
when `family$d2mu.deta2` is missing, this function returns `NULL` and
the solver falls back to numeric/Broyden Jacobians.
