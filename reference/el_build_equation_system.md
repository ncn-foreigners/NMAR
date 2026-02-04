# Empirical likelihood estimating equations for SRS

Returns a function that evaluates the stacked EL system for \\\theta =
(\beta, z, \lambda_x)\\ with \\z = \operatorname{logit}(W)\\. Blocks
correspond to:

1.  missingness model score equations in \\\beta\\,

2.  the response-rate equation in \\W\\,

3.  auxiliary moment constraints in \\\lambda_x\\.

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

When no auxiliaries are present the last block is omitted. The system
matches QLS equations 7-10. We cap \\\eta\\, clip \\w_i\\ in ratios, and
guard \\D_i\\ away from zero to ensure numerical stability.

**Guarding policy:**

- Cap \\\eta\\: `eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())`.

- Compute `w <- family$linkinv(eta)` and clip to `[1e-12, 1 - 1e-12]`
  when used in ratios.

- Denominator floor: `Di <- pmax(Di_raw, nmar_get_el_denom_floor())`. In
  the Jacobian, terms that depend on `d(1/Di)/d(.)` are multiplied by
  `active = 1(Di_raw > floor)` to match the clamped equations.
