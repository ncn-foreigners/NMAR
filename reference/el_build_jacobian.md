# Analytical Jacobian for empirical likelihood

Analytical Jacobian for empirical likelihood

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

## Details

Builds the block Jacobian \\A = \partial F/\partial \theta\\ for the EL
system with \\\theta = (\beta, z, \lambda_x)\\ and \\z =
\operatorname{logit}(W)\\. Blocks follow Qin, Leung, and Shao (2002,
Eqs. 7-10). The derivative with respect to the linear predictor for the
missingness (response) model uses the Bernoulli score form
\\\partial/\partial\eta\\ \log w(\eta) = \mu.\eta(\eta)/w(\eta)\\ with
link-inverse clipping. Denominator guards are applied consistently when
forming terms depending on \\D_i(\theta)\\.

Guarding policy (must remain consistent across
equations/Jacobian/post): - Cap eta: eta \<- pmax(pmin(eta,
get_eta_cap()), -get_eta_cap()) - Compute w \<- family\$linkinv(eta);
clip to \[1e-12, 1-1e-12\] when used in ratios - Denominator floor: Di
\<- pmax(Di_raw, nmar_get_el_denom_floor()); multiply terms that depend
on d(1/Di)/d(.) by active = 1(Di_raw \> floor)

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.
