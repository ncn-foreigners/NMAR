# Core Empirical Likelihood Estimator

Implements the core computational engine for empirical likelihood
estimation under nonignorable nonresponse, including parameter solving,
variance calculation, and diagnostic computation.

## Usage

``` r
el_estimator_core(
  full_data,
  respondent_data,
  respondent_weights,
  N_pop,
  internal_formula,
  auxiliary_means,
  standardize,
  trim_cap,
  control,
  on_failure,
  family = logit_family(),
  variance_method,
  bootstrap_reps,
  user_args,
  start = NULL,
  trace_level = 0,
  ...
)
```

## Arguments

- full_data:

  Data frame or survey design object containing all units.

- respondent_data:

  Data frame containing only responding units.

- respondent_weights:

  Numeric vector of base sampling weights for respondents.

- N_pop:

  Numeric. Total population size (weighted if survey design).

- internal_formula:

  List of internal formulas for outcome, response, and auxiliary models.

- auxiliary_means:

  Named numeric vector of known population means.

- standardize:

  Logical. Whether to standardize predictors during estimation.

- trim_cap:

  Numeric. Upper bound for empirical likelihood weight trimming.

- control:

  List of control parameters for the nonlinear equation solver.

- on_failure:

  Character. Action when solver fails: "return" or "error".

- family:

  List. Link function specification (typically logit).

- variance_method:

  Character. Variance estimation method.

- bootstrap_reps:

  Integer. Number of bootstrap replications.

- user_args:

  List. Original user arguments for bootstrap replication.

- ...:

  Additional arguments passed to the solver.

## Value

List containing estimation results, diagnostics, and metadata.

## Details

Orchestrates EL estimation for NMAR following Qin, Leung, and Shao
(2002). The stacked system in \\(\beta, z, \lambda_x)\\ with \\z =
\logit(W)\\ is solved by `nleqslv` using an analytic Jacobian. Numerical
safeguards are applied consistently across equations, Jacobian, and
post-solution weights: bounded linear predictors, probability clipping
in ratios, and a small floor on denominators `D_i()` with an active-set
mask in derivatives. After solving, unnormalized masses `d_i/D_i()` are
formed, optional trimming may be applied (with normalization only for
reporting), and optional variance is computed via bootstrap. The
analytical delta method for EL is not implemented.
