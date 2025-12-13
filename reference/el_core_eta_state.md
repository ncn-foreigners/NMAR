# Core eta-state computation for EL engines

Computes the capped linear predictor, response probabilities,
derivatives, and stable scores with respect to the linear predictor for
a given family. This helper centralizes the numerically delicate pieces
(capping, clipping, Mills ratios, and score derivatives) and is used
consistently across the EL equation system and analytical Jacobians for
both IID and survey designs.

## Usage

``` r
el_core_eta_state(family, eta_raw, eta_cap)
```

## Arguments

- family:

  List-like response family bundle (see
  [`logit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/logit_family.md)
  and
  [`probit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/probit_family.md)).

- eta_raw:

  Numeric vector of unconstrained linear predictors.

- eta_cap:

  Scalar cap applied symmetrically to `eta_raw`.

## Value

A list with components:

- `eta`:

  Capped linear predictor.

- `w`:

  Mean function `family$linkinv(eta)`.

- `w_clipped`:

  `w` clipped to `[1e-12, 1-1e-12]` for use in ratios.

- `mu_eta`:

  Derivative `family$mu.eta(eta)`.

- `d2mu`:

  Second derivative `family$d2mu.deta2(eta)` when available, otherwise
  `NULL`.

- `s_eta`:

  Score with respect to `eta`, using stable logit/probit forms where
  possible.

- `ds_eta_deta`:

  Derivative of `s_eta` with respect to `eta` when `d2mu` is available,
  otherwise `NULL`.
