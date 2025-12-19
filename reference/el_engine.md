# Empirical likelihood (EL) engine for NMAR

Constructs an engine specification for the empirical likelihood (EL)
estimator of a full-data mean under nonignorable nonresponse (NMAR).

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
  first RHS (after formula expansion), e.g., factor indicator columns
  created by
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) or
  transformed terms like `I(X^2)`. Auxiliary intercepts are always
  dropped automatically, so do not supply `(Intercept)`. If `NULL`
  (default) and the outcome contains at least one `NA`, auxiliary means
  are estimated from the full input (including nonrespondents): IID uses
  unweighted column means of the auxiliary design; survey designs use
  the design-weighted means based on `weights(design)`. This corresponds
  to the QLS case where \\\mu_x\\ is replaced by \\\bar X\\ (the
  full-sample mean) when auxiliary variables are observed for all
  sampled units.

- control:

  Optional solver configuration forwarded to
  [`nleqslv::nleqslv()`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html).
  Provide a single list that may include solver tolerances (e.g.,
  `xtol`, `ftol`, `maxit`) and, optionally, top-level entries `global`
  and `xscalm` for globalization and scaling. Example:
  `control = list(maxit = 500, xtol = 1e-10, ftol = 1e-10, global = "qline", xscalm = "auto")`.

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
  to the total design-weight total on the same analysis scale as
  `weights(design)` (default `sum(weights(design))`). If omitted and the
  outcome contains no NAs, the estimator errors, requesting `n_total`.

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

  Missingness (response) model family. Either `"logit"` (default) or
  `"probit"`, or a custom family object: a list with components `name`,
  `linkinv`, `mu.eta`, `score_eta`, and optionally `d2mu.deta2`. When
  `d2mu.deta2` is absent the solver uses Broyden/numeric Jacobians.

## Value

A list of class `"nmar_engine_el"` (also inheriting from
`"nmar_engine"`) containing configuration fields to be supplied to
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).
Users rarely access fields directly; instead, pass the engine to
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
together with a formula and data.

## Details

The implementation follows Qin, Leung, and Shao (2002): the response
mechanism is modeled as \\w(y, x; \beta) = P(R = 1 \mid Y = y, X = x)\\
and the joint law of \\(Y, X)\\ is represented nonparametrically by
respondent masses that satisfy empirical likelihood constraints. The
mean is estimated as a respondent weighted mean with weights
proportional to \\\tilde w_i = a_i / D_i(\beta, W, \lambda)\\, where
\\a_i\\ are base weights (\\a_i \equiv 1\\ for IID data and \\a_i =
d_i\\ for survey designs) and \\D_i\\ is the EL denominator.

For `data.frame` inputs the estimator solves the Qin-Leung-Shao (QLS)
estimating equations for \\(\beta, W, \lambda_x)\\ with \\W\\
reparameterized as \\z = \operatorname{logit}(W)\\, and profiles out the
response multiplier \\\lambda_W\\ using the closed-form QLS identity
(their Eq. 10). For `survey.design` inputs the estimator uses a
design-weighted analogue (Chen and Sitter 1999; Wu 2005) with an
explicit \\\lambda_W\\ and an additional linkage equation involving the
nonrespondent design-weight total \\T_0\\.

Numerical stability:

- \\W\\ is optimized on the logit scale so \\0 \< W \< 1\\.

- The response-model linear predictor is capped and EL denominators
  \\D_i\\ are floored at a small positive value; the analytic Jacobian
  is consistent with this guard via an active-set mask.

- Optional trimming (`trim_cap`) is applied only after solving, to the
  unnormalized masses \\\tilde w_i = a_i/D_i\\; this changes the
  returned weights and therefore the point estimate.

**Formula syntax and data constraints**:
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
accepts a partitioned right-hand side
`y_miss ~ auxiliaries | response_only`. Variables left of `|` enter
auxiliary moment constraints; variables right of `|` enter only the
response model. The outcome (LHS) is always included as a response-model
predictor through the evaluated LHS expression; explicit use of the
outcome on the RHS is rejected. The response model always includes an
intercept; the auxiliary block never includes an intercept.

To include a covariate in both the auxiliary constraints and the
response model, repeat it on both sides, e.g. `y_miss ~ X | X`.

**Auxiliary means**: If `auxiliary_means = NULL` (default) and the
outcome contains at least one `NA`, auxiliary means are estimated from
the full input and used as \\\bar X\\ in the QLS constraints. For
respondents-only data (no `NA` in the outcome), `n_total` must be
supplied; and if the auxiliary RHS is non-empty, `auxiliary_means` must
also be supplied. When `standardize = TRUE`, supply `auxiliary_means` on
the original data scale; the engine applies the same standardization
internally.

**Survey scale**: For `survey.design` inputs, `n_total` (if provided)
must be on the same analysis scale as `weights(design)`. The default is
`sum(weights(design))`.

**Convergence and identification**: the stacked EL system can have
multiple solutions. Adding response-only predictors (variables to the
right of `|`) can make the problem sensitive to starting values. Inspect
diagnostics such as `jacobian_condition_number` and consider supplying
`start = list(beta = ..., W = ...)` when needed.

**Variance**: The EL engine supports bootstrap standard errors via
`variance_method = "bootstrap"` or can skip variance with
`variance_method = "none"`. Set a seed for reproducible bootstrap
results.

Bootstrap requires suggested packages: for IID resampling it requires
`future.apply` (and `future`); for survey replicate-weight bootstrap it
requires `survey` and `svrep`.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. Journal of the
American Statistical Association, 97(457), 193-200.
[doi:10.1198/016214502753479338](https://doi.org/10.1198/016214502753479338)

Chen, J., and Sitter, R. R. (1999). A pseudo empirical likelihood
approach for the effective use of auxiliary information in complex
surveys. Statistica Sinica, 9, 385-406.

Wu, C. (2005). Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling. Survey Methodology, 31(2),
239-243.

## See also

[`nmar`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md),
[`weights.nmar_result`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/weights.nmar_result.md),
[`summary.nmar_result`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/summary.nmar_result.md)

## Examples

``` r
set.seed(1)
n <- 200
X <- rnorm(n)
Y <- 2 + 0.5 * X + rnorm(n)
p <- plogis(-0.7 + 0.4 * scale(Y)[, 1])
R <- runif(n) < p
if (all(R)) R[1] <- FALSE
df <- data.frame(Y_miss = Y, X = X)
df$Y_miss[!R] <- NA_real_

# Estimate auxiliary mean from full data (QLS "use Xbar" case)
eng <- el_engine(auxiliary_means = NULL, variance_method = "none")

# Put X in both the auxiliary block and the response model (QLS-like)
fit <- nmar(Y_miss ~ X | X, data = df, engine = eng)
summary(fit)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 2.255538
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 200 
#> Respondents: 76 
#> Call: nmar(Y_miss ~ X | X, data = <data.frame: N=200>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -1.410187
#> Y_miss       0.396554
#> X            0.138009

# \donttest{
# Response-only predictors can be placed to the right of |:
Z <- rnorm(n)
df2 <- data.frame(Y_miss = Y, X = X, Z = Z)
df2$Y_miss[!R] <- NA_real_
eng2 <- el_engine(auxiliary_means = NULL, variance_method = "none")
fit2 <- nmar(Y_miss ~ X | Z, data = df2, engine = eng2)
print(fit2)
#> Call: nmar(Y_miss ~ X | Z, data = <data.frame: N=200>, engine = empirical_likelihood)
#> 
#> NMAR Result
#> ------------
#> Y_miss mean: 2.107642
#> Converged: TRUE 
#> Variance method: none 
#> Estimator: empirical_likelihood 
#> Sample size: 200 (respondents: 76)
#> 
#> Method: empirical_likelihood
#> Max equation residual: 8.833e-09
#> Constraint sum (W): 8.833e-09
#> Constraint sums (aux):
#>            X 
#> 1.796749e-09 

# Survey design usage
if (requireNamespace("survey", quietly = TRUE)) {
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  eng3 <- el_engine(auxiliary_means = NULL, variance_method = "none")
  fit3 <- nmar(Y_miss ~ X, data = des, engine = eng3)
  summary(fit3)
}
#> NMAR Model Summary
#> =================
#> Y_miss mean: 2.099150
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 200 
#> Respondents: 76 
#> Call: nmar(Y_miss ~ X, data = <survey.design: N=200>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -2.113955
#> Y_miss       0.747676
# }
```
