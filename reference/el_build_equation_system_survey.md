# Empirical likelihood equations for survey designs (design-weighted QLS system)

Empirical likelihood equations for survey designs (design-weighted QLS
system)

## Usage

``` r
el_build_equation_system_survey(
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

Returns a function that evaluates the stacked EL system for complex
survey designs using design weights. Unknowns are \\\theta = (\beta, z,
\lambda_W, \lambda_x)\\ with \\z = \operatorname{logit}(W)\\. Blocks
correspond to:

- response-model score equations in \\\beta\\,

- the response-rate equation in \\W\\ based on \\\sum d_i (w_i - W)/D_i
  = 0\\,

- auxiliary moment constraints \\\sum d_i (X_i - \mu_x)/D_i = 0\\,

- and the design-based linkage between \\\lambda_W\\ and the
  nonrespondent total: \\T_0/(1-W) - \lambda_W \sum d_i / D_i = 0\\,
  where \\T_0 = N\_{\mathrm{pop}} - \sum d_i\\ on the analysis scale.

When all design weights are equal and \\N\_{\mathrm{pop}}\\ and the
respondent count match the simple random sampling setup, this system
reduces to the Qin, Leung, and Shao (2002) equations (6)-(10).
