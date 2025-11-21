# Core Empirical Likelihood Estimator

Implements the core computational engine for empirical likelihood
estimation under nonignorable nonresponse, including parameter solving,
variance calculation, and diagnostic computation.

## Usage

``` r
el_estimator_core(
  missingness_design,
  aux_matrix,
  aux_means,
  respondent_weights,
  analysis_data,
  outcome_expr,
  N_pop,
  formula,
  standardize,
  trim_cap,
  control,
  on_failure,
  family = logit_family(),
  variance_method,
  bootstrap_reps,
  start = NULL,
  trace_level = 0,
  auxiliary_means = NULL
)
```

## Arguments

- missingness_design:

  Respondent-side missingness (response) model design matrix
  (intercept + predictors).

- aux_matrix:

  Auxiliary design matrix on respondents (may have zero columns).

- aux_means:

  Named numeric vector of auxiliary population means (aligned to columns
  of \`aux_matrix\`).

- respondent_weights:

  Numeric vector of respondent weights aligned with
  \`missingness_design\` rows.

- analysis_data:

  Data object used for logging and variance (survey designs supply the
  design object).

- outcome_expr:

  Character string identifying the outcome expression displayed in
  outputs.

- N_pop:

  Population size on the analysis scale.

- formula:

  Original model formula used for estimation.

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

- auxiliary_means:

  Named numeric vector of known population means supplied by the user
  (optional; used for diagnostics).

## Value

List containing estimation results, diagnostics, and metadata.

## Details

Orchestrates EL estimation for NMAR following Qin, Leung, and Shao
(2002). For `data.frame` inputs (IID setting) the stacked system in
\\(\beta, z, \lambda_x)\\ with \\z = \logit(W)\\ is solved by `nleqslv`
using an analytic Jacobian. For `survey.design` inputs a design-weighted
analogue in \\(\beta, z, \lambda_W, \lambda_x)\\ is solved with an
analytic Jacobian when the response family supplies second derivatives,
or with numeric/Broyden Jacobians otherwise. Numerical safeguards are
applied consistently across equations, Jacobian, and post-solution
weights: bounded linear predictors, probability clipping in ratios, and
a small floor on denominators \\D_i(\theta)\\ with an active-set mask in
derivatives. After solving, unnormalized masses \\d_i/D_i(\theta)\\ are
formed, optional trimming may be applied (with normalization only for
reporting), and optional variance is computed via bootstrap. The
analytical delta method for EL is not implemented.
