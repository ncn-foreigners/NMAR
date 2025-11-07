# Post-solution: denominators, masses, and mean (EL)

Computes EL denominators and masses and the mean estimate after the
solver converges. The mean is calculated from normalized probability
masses per Qin, Leung, and Shao (2002, Eq. 11). When trimming is active,
the masses are capped and renormalized before computing the mean;
diagnostics still report constraint sums based on untrimmed masses.
Denominator floors and the guarding policy used in equations/Jacobian
are applied consistently here for diagnostic coherence.

## Usage

``` r
el_post_solution(
  estimates,
  response_model_matrix_scaled,
  response_model_matrix_unscaled,
  auxiliary_matrix_scaled,
  mu_x_scaled,
  respondent_data,
  outcome_var,
  family,
  N_pop,
  respondent_weights,
  K_beta,
  K_aux,
  nmar_scaling_recipe,
  standardize,
  trim_cap,
  X_centered = NULL
)
```

## Arguments

- estimates:

  Numeric vector (beta, z, lambda) at the solution.

- response_model_matrix_scaled:

  Scaled design matrix for the response model.

- response_model_matrix_unscaled:

  Unscaled design matrix for the response model.

- auxiliary_matrix_scaled:

  Scaled auxiliary matrix (or empty matrix).

- mu_x_scaled:

  Vector of population means for scaled auxiliaries (or NULL).

- respondent_data:

  Data frame of respondents.

- outcome_var:

  Character; outcome column name in respondent_data.

- family:

  Family object with linkinv and mu.eta.

- N_pop:

  Numeric; population total on the analysis scale.

- respondent_weights:

  Base weights for respondents.

- K_beta, K_aux:

  Integers; sizes of beta and lambda.

- nmar_scaling_recipe:

  Scaling recipe object for unscaling.

- standardize:

  Logical; whether coefficients need unscaling.

- trim_cap:

  Numeric; weight trimming cap (Inf = no trimming).
