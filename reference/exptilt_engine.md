# Exponential tilting (ET) engine for NMAR estimation

Constructs a configuration for the exponential tilting estimator under
nonignorable nonresponse (NMAR). The estimator solves
\\S_2(\boldsymbol{\phi}, \hat{\boldsymbol{\gamma}}) = 0,\\ using nleqslv
to apply EM algorithm.

## Usage

``` r
exptilt_engine(
  standardize = FALSE,
  on_failure = c("return", "error"),
  variance_method = c("delta", "bootstrap", "none"),
  bootstrap_reps = 10,
  supress_warnings = FALSE,
  auxiliary_means = NULL,
  control = list(),
  family = c("logit", "probit"),
  y_dens = c("normal", "lognormal", "exponential"),
  stopping_threshold = 1,
  sample_size = 2000
)
```

## Arguments

- standardize:

  logical; standardize predictors. Default `TRUE`.

- on_failure:

  character; `"return"` or `"error"` on solver failure

- variance_method:

  character; one of `"delta"`, `"bootstrap"`, or `"none"`.

- bootstrap_reps:

  integer; number of bootstrap replicates when
  `variance_method = "bootstrap"`.

- supress_warnings:

  Logical; suppress variance-related warnings.

- auxiliary_means:

  Optional named numeric vector of population moments for auxiliary
  covariates.

- control:

  Named list of control parameters passed to
  [`nleqslv::nleqslv`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html).
  Common parameters include:

  - `maxit`: Maximum number of iterations (default: 100)

  - `method`: Solver method - `"Newton"` or `"Broyden"` (default:
    `"Newton"`)

  - `global`: Global strategy - `"dbldog"`, `"pwldog"`, `"qline"`,
    `"gline"`, `"hook"`, or `"none"` (default: `"dbldog"`)

  - `xtol`: Tolerance for relative error in solution (default: 1e-8)

  - `ftol`: Tolerance for function value (default: 1e-8)

  - `btol`: Tolerance for backtracking (default: 0.01)

  - `allowSingular`: Allow singular Jacobians (default: `TRUE`)

  See
  [`?nleqslv::nleqslv`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html)
  for full details.

- family:

  character; response model family, either `"logit"` or `"probit"`, or a
  family object created by
  [`logit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/logit_family.md)
  /
  [`probit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/probit_family.md).

- y_dens:

  Outcome density model (`"auto"`, `"normal"`, `"lognormal"`, or
  `"exponential"`).

- stopping_threshold:

  Numeric; early stopping threshold. If the maximum absolute value of
  the score function falls below this threshold, the algorithm stops
  early (default: 1).

- sample_size:

  Integer; maximum sample size for stratified random sampling (default:
  2000). When the dataset exceeds this size, a stratified random sample
  is drawn to optimize memory usage. The sampling preserves the ratio of
  respondents to non-respondents in the original data.

## Value

An engine object of class `c("nmar_engine_exptilt","nmar_engine")`. This
is a configuration list; it is not a fit. Pass it to
[nmar](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).

## Details

The method is a robust Propensity-Score Adjustment (PSA) approach for
Not Missing at Random (NMAR). It uses Maximum Likelihood Estimation
(MLE), basing the likelihood on the observed part of the sample
(\\f(\boldsymbol{Y}\_i \| \delta_i = 1, \boldsymbol{X}\_i)\\), making it
robust against outcome model misspecification. The propensity score is
estimated by assuming an instrumental variable (\\X_2\\) that is
independent of the response status given other covariates and the study
variable. Estimator calculates fractional imputation weights \\w_i\\.
The final estimator is a weighted average, where the weights are the
inverse of the estimated response probabilities \\\hat{\pi}\_i\\,
satisfying the estimating equation: \$\$ \sum\_{i \in \mathcal{R}}
\frac{\boldsymbol{g}(\boldsymbol{Y}\_i, \boldsymbol{X}\_i ;
\boldsymbol{\theta})}{\hat{\pi}\_i} = 0, \$\$ where \\\mathcal{R}\\ is
the set of observed respondents.

## References

Minsun Kim Riddles, Jae Kwang Kim, Jongho Im A
Propensity-score-adjustment Method for Nonignorable Nonresponse *Journal
of Survey Statistics and Methodology*, Volume 4, Issue 2, June 2016,
Pages 215–245.

## Examples

``` r
# \donttest{
generate_test_data <- function(
  n_rows = 500,
  n_cols = 1,
  case = 1,
  x_var = 0.5,
  eps_var = 0.9,
  a = 0.8,
  b = -0.2
) {
# Generate X variables - fixed to match comparison
  X <- as.data.frame(replicate(n_cols, rnorm(n_rows, 0, sqrt(x_var))))
  colnames(X) <- paste0("x", 1:n_cols)

# Generate Y - fixed coefficients to match comparison
  eps <- rnorm(n_rows, 0, sqrt(eps_var))
  if (case == 1) {
# Use fixed coefficient of 1 for all x variables to match: y = -1 + x1 + epsilon
    X$Y <- as.vector(-1 + as.matrix(X) %*% rep(1, n_cols) + eps)
  }
  else if (case == 2) {
    X$Y <- -2 + 0.5 * exp(as.matrix(X) %*% rep(1, n_cols)) + eps
  }
  else if (case == 3) {
    X$Y <- -1 + sin(2 * as.matrix(X) %*% rep(1, n_cols)) + eps
  }
  else if (case == 4) {
    X$Y <- -1 + 0.4 * as.matrix(X)^3 %*% rep(1, n_cols) + eps
  }

  Y_original <- X$Y

# Missingness mechanism - identical to comparison
  pi_obs <- 1 / (1 + exp(-(a + b * X$Y)))

# Create missing values
  mask <- runif(nrow(X)) > pi_obs
  mask[1] <- FALSE # Ensure at least one observation is not missing
  X$Y[mask] <- NA

  return(list(X = X, Y_original = Y_original))
}
res_test_data <- generate_test_data(n_rows = 500, n_cols = 1, case = 1)
x <- res_test_data$X

exptilt_config <- exptilt_engine(
  y_dens = 'normal',
  control = list(maxit = 10),
  stopping_threshold = 0.01,
  standardize = FALSE,
  family = 'logit',
  bootstrap_reps = 50
)
formula = Y ~ x1
res <- nmar(formula = formula, data = x, engine = exptilt_config, trace_level = 1)
#> [STEP-1] ============================================================ 
#> [STEP-1]   EXPTILT ESTIMATION STARTED 
#> [STEP-1] ============================================================ 
#> [INFO] Running with trace_level = 1 | For more detailed output, use trace_level = 2. Available trace_level = c(1,2,3) 
#> [INFO] Formula: Y ~ x1 
#> [INFO]  
#> [INFO] -- DATA SUMMARY -- 
#> [INFO]   Total observations:      500 
#> [INFO]   Respondents:             355 (29.0%) 
#> [INFO]   Non-respondents:         145 (71.0%) 
#> [INFO]  
#> [INFO] -- CONDITIONAL DENSITY ESTIMATION -- 
#> [INFO]   Selected distribution:   normal 
#> [INFO]  
#> [INFO] -- NONLINEAR SOLVER (nleqslv) -- 
#> [INFO]   Solving... 
#> [INFO]  
#> [INFO]   OK Converged 
#> [INFO]   Iterations:               5 
#> [INFO]  
#> [INFO] -- VARIANCE ESTIMATION (Delta Method) -- 
#> Warning: Delta variance failed to evaluate; using bootstrap instead.
#> [INFO]  
#> [INFO] -- VARIANCE ESTIMATION (Bootstrap) -- 
#> [INFO]   Bootstrap replications:   50 
#> [INFO]   OK Bootstrap complete 
#> [INFO]   Standard error:           0.000003 
#> [INFO]  
#> [RESULT] ============================================================ 
#> [RESULT]   ESTIMATION COMPLETE 
#> [RESULT] ============================================================ 
#> [RESULT]   Mean estimate:            -1.007480 
#> [RESULT]   Standard error:           0.000003 
#> [RESULT]   95% CI:                   [-1.007485, -1.007476] 
#> [RESULT] ============================================================ 
summary(res)
#> NMAR Model Summary
#> =================
#> Y mean: -1.007480 (0.000003)
#> 95% CI: (-1.007485, -1.007476)
#> Converged: TRUE 
#> Variance method: delta 
# }
```
