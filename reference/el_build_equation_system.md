# Empirical likelihood estimating equations

Empirical likelihood estimating equations

## Usage

``` r
el_build_equation_system(
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

Returns a function that evaluates the stacked EL system for \\\theta =
(\beta, z, \lambda_x)\\ with \\z = \operatorname{logit}(W)\\. Blocks
correspond to: (i) missingness (response) model score equations in
\\\beta\\, (ii) the response-rate equation in \\W\\, and (iii) auxiliary
moment constraints in \\\lambda_x\\. When no auxiliaries are present the
last block is omitted. The system matches Qin, Leung, and Shao (2002,
Eqs. 7-10) with empirical masses \\m_i = d_i/D_i(\theta)\\, \\D_i\\ as
in the paper. We cap \\\eta\\, clip \\p\\, and guard \\D_i\\ away from
zero to ensure numerical stability; these safeguards are applied
consistently in equations, Jacobian, and post-solution weights.

Guarding policy (must remain consistent across
equations/Jacobian/post): - Cap eta: eta \<- pmax(pmin(eta,
get_eta_cap()), -get_eta_cap()) - Compute w \<- family\$linkinv(eta);
clip to \[1e-12, 1-1e-12\] when used in ratios - Denominator floor: Di
\<- pmax(Di_raw, nmar_get_el_denom_floor()); in the Jacobian, multiply
terms that depend on d(1/Di)/d(.) by active = 1(Di_raw \> floor)

The score with respect to the linear predictor uses the Bernoulli form
\\s\_{\eta,i}(\beta) = \partial \log w_i / \partial \eta_i =
\mu.\eta(\eta_i)/w_i\\, which is valid for both logit and probit links
when \\w_i\\ is clipped.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.
