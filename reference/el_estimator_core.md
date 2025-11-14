# Core Empirical Likelihood Estimator

Implements the core computational engine for empirical likelihood
estimation under nonignorable nonresponse, including parameter solving,
variance calculation, and diagnostic computation.

## Usage

``` r
el_estimator_core(
  missingness_design,
  auxiliary_matrix,
  mu_x,
  respondent_weights,
  full_data,
  outcome_var,
  N_pop,
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
  auxiliary_means = NULL
)
```

## Arguments

- missingness_design:

  Respondent-side missingness (response) model design matrix
  (intercept + predictors).

- auxiliary_matrix:

  Auxiliary design matrix on respondents (may have zero columns).

- mu_x:

  Named numeric vector of auxiliary population means (aligned to columns
  of \`auxiliary_matrix\`).

- respondent_weights:

  Numeric vector of respondent weights aligned with
  \`missingness_design\` rows.

- full_data:

  Data object used for logging (survey designs supply the design
  object).

- outcome_var:

  Character string identifying the outcome expression displayed in
  outputs.

- N_pop:

  Population size on the analysis scale.

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

- auxiliary_means:

  Named numeric vector of known population means supplied by the user
  (optional; used for diagnostics).

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
