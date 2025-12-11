# Empirical likelihood (EL) engine for NMAR

Constructs a configuration object for the empirical likelihood estimator
under nonignorable nonresponse (NMAR) with optional auxiliary moment
constraints. For `data.frame` inputs (IID setting) the estimator solves
a stacked system in \\\theta = (\beta, z, \lambda_x)\\ with \\z =
\operatorname{logit}(W)\\ using a Newton method with an analytic
Jacobian and globalization via
[nleqslv](https://rdrr.io/pkg/nleqslv/man/nleqslv.html). For
`survey.design` inputs it solves a design-weighted analogue in \\\theta
= (\beta, z, \lambda_W, \lambda_x)\\. When the response family supplies
second derivatives (logit and probit) an analytic Jacobian is used;
otherwise the solver falls back to numeric/Broyden Jacobians. Numerical
safeguards (bounded linear predictor, link-inverse clipping, denominator
floors, and stable linear algebra) improve robustness. Pass the engine
to
[nmar](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
together with a formula and data.

## Usage

``` r
el_engine(
  standardize = TRUE,
  trim_cap = Inf,
  on_failure = c("return", "error"),
  variance_method = c("bootstrap", "none"),
  bootstrap_reps = 500,
  auxiliary_means = NULL,
  control = list(),
  strata_augmentation = TRUE,
  n_total = NULL,
  start = NULL,
  family = c("logit", "probit")
)
```

## Arguments

- standardize:

  logical; standardize predictors. Default `TRUE`.

- trim_cap:

  numeric; cap for EL weights (`Inf` = no trimming).

- on_failure:

  character; `"return"` or `"error"` on solver failure.

- variance_method:

  character; one of `"bootstrap"` or `"none"`.

- bootstrap_reps:

  integer; number of bootstrap replicates when
  `variance_method = "bootstrap"`.

- auxiliary_means:

  named numeric vector; population means for auxiliary design columns.
  Names must match the materialized model.matrix column names on the
  first RHS (after formula expansion), e.g., factor indicators like
  \`F_b\` or transformed terms \`I(X^2)\`. Auxiliary intercepts are
  always dropped automatically, so do not supply \`(Intercept)\`.
  Optional.

- control:

  list; optional solver control for
  [`nleqslv::nleqslv()`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html).
  Recognized fields (defaults in parentheses):

  - Top-level: `global` = `"qline"` (quadratic line search) or one of
    `"dbldog"`, `"pwldog"`, `"cline"`, `"gline"`, `"hook"`, `"none"`;
    `xscalm` = `"auto"` or `"fixed"`

  - In `control=`: `xtol`, `ftol`, `btol`, `maxit`, `trace`, `stepmax`,
    `delta`, `allowSing`

  Unknown names are ignored. For `data.frame` inputs the EL system is
  solved by Newton with an analytic Jacobian; for `survey.design` inputs
  a design-weighted analogue is solved with an analytic Jacobian when
  available or numeric/Broyden Jacobians otherwise.

- strata_augmentation:

  logical; when `TRUE` (default), survey designs with an identifiable
  strata structure are augmented with stratum indicators and
  corresponding population shares in the auxiliary block (Wu-style
  strata augmentation). Has no effect for `data.frame` inputs or survey
  designs without strata.

- n_total:

  numeric; optional when supplying respondents-only data (no `NA` in the
  outcome). For `data.frame` inputs, set to the total number of sampled
  units before filtering to respondents. For `survey.design` inputs, set
  to the total design weight or known population total. If omitted and
  the outcome contains no NAs, the estimator errors, requesting
  `n_total`.

- start:

  list; optional starting point for the solver. Fields:

  - `beta`: named numeric vector of missingness-model coefficients on
    the original (unscaled) scale, including `(Intercept)`.

  - `W` or `z`: starting value for population response rate
    (`0 < W < 1`) or its logit (`z`). If both are provided, `z` takes
    precedence.

  - `lambda`: named numeric vector of auxiliary multipliers on the
    original scale (names must match auxiliary design columns; no
    intercept). Values are mapped to the scaled space internally.

- family:

  character; missingness (response) model family, either `"logit"` or
  `"probit"`, or a family object created by
  [`logit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/logit_family.md)
  /
  [`probit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/probit_family.md).

## Value

A list of class `"nmar_engine_el"` (also inheriting from
`"nmar_engine"`) containing configuration fields to be supplied to
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).
Users rarely access fields directly; instead, pass the engine to
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
together with a formula and data.

## Details

This engine implements an empirical likelihood estimator for NMAR
response based on Qin, Leung and Shao (2002) for IID data, and a
design-weighted analogue for complex survey designs inspired by Chen and
Sitter (1999) and Wu (2005). For `data.frame` inputs the unknowns are
`(beta, z, lambda_x)` with `z = logit(W)`, and the QLS closed-form
identity is used to profile out the multiplier `lambda_W`. For
`survey.design` inputs the system is extended to
`(beta, z, lambda_W, lambda_x)` and solved with design weights and, when
present, Wu-style strata augmentation in the auxiliary block. Numerical
guards (capped linear predictors, clipped response probabilities,
denominator floors) are applied consistently in equations and Jacobians.

**Formula syntax**:
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
supports a partitioned right-hand side `y_miss ~ aux1 + aux2 | z1 + z2`.
Variables left of `|` are auxiliaries (used in EL moment constraints);
variables right of `|` are missingness-model predictors only. The
outcome appears on the left-hand side and is included as a response
predictor by default. Auxiliary design matrices are constructed with an
intercept dropped automatically; missingness models always include an
intercept even if the formula uses `-1` or `+0`.

**Variance**: The EL engine supports bootstrap standard errors via
`variance_method = "bootstrap"` or can skip variance with
`variance_method = "none"`.

## Progress Reporting

When `variance_method = "bootstrap"`, progress reporting is available
via the `progressr` package. To enable it:

    library(progressr)
    library(future)

    # Enable progress reporting
    handlers(global = TRUE)
    handlers("txtprogressbar")  # or "progress", "cli", etc.

    # Set parallel backend (optional)
    plan(multisession, workers = 4)

    # Always set seed for reproducibility
    set.seed(123)

    # Run with progress bar
    result <- nmar(Y ~ X, data = df,
                   engine = el_engine(variance_method = "bootstrap",
                                      bootstrap_reps = 500))

    # Reset to sequential
    plan(sequential)

To disable progress in simulations or batch jobs:

`handlers("void") # Silent`

If progressr is not installed or no handlers are set, bootstrap runs
silently (default behavior). Progress reporting works with all future
backends and does not affect reproducibility.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.

Chen, J., and Sitter, R. R. (1999). A pseudo empirical likelihood
approach for complex survey data. Biometrika, 86(2), 373-385.

Wu, C. (2005). Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling. Canadian Journal of Statistics,
33(3), 497-509.

## See also

\[nmar()\], \[weights.nmar_result()\], \[summary.nmar_result\]

## Examples

``` r
# \donttest{
set.seed(1)
n <- 200
X <- rnorm(n)
Z <- rnorm(n)
Y <- 2 + 0.5 * X + Z
p <- plogis(-0.7 + 0.4 * scale(Y)[, 1])
R <- runif(n) < p
df <- data.frame(Y_miss = Y, X = X)
df$Y_miss[!R] <- NA_real_
eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
fit <- nmar(Y_miss ~ X, data = df, engine = eng)
summary(fit)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.979080
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 200 
#> Respondents: 76 
#> Call: nmar(Y_miss ~ X, data = <data.frame: N=200>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -2.535301
#> Y_miss       0.987603

# Response-only predictors can be placed to the right of `|`:
df2 <- data.frame(Y_miss = Y, X = X, Z = Z)
df2$Y_miss[!R] <- NA_real_
eng2 <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
fit2 <- nmar(Y_miss ~ X | Z, data = df2, engine = eng2)
print(fit2)
#> Call: nmar(Y_miss ~ X | Z, data = <data.frame: N=200>, engine = empirical_likelihood)
#> 
#> NMAR Result
#> ------------
#> Y_miss mean: 2.372442
#> Converged: TRUE 
#> Variance method: none 
#> Estimator: empirical_likelihood 
#> Sample size: 200 (respondents: 76)
#> 
#> Method: empirical_likelihood
#> Max equation residual: 8.380e-09
#> Constraint sum (W): -8.380e-09
#> Constraint sums (aux):
#>             X 
#> -3.356769e-09 

# Survey design usage
if (requireNamespace("survey", quietly = TRUE)) {
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  eng3 <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  fit3 <- nmar(Y_miss ~ X, data = des, engine = eng3)
  summary(fit3)
}
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.979080
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 200 
#> Respondents: 76 
#> Call: nmar(Y_miss ~ X, data = <survey.design: N=200>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -2.535301
#> Y_miss       0.987603
# }
```
