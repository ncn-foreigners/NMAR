# Empirical likelihood (EL) engine for NMAR

Constructs a configuration object for the empirical likelihood estimator
under nonignorable nonresponse (NMAR) with optional auxiliary moment
constraints. The estimator solves the stacked system in \\\theta =
(\beta, z, \lambda_x)\\ with \\z = \operatorname{logit}(W)\\ using a
Newton method with analytic Jacobian and globalization via
[nleqslv](https://rdrr.io/pkg/nleqslv/man/nleqslv.html). Numerical
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
  variance_method = c("delta", "bootstrap", "none"),
  bootstrap_reps = 500,
  auxiliary_means = NULL,
  control = list(),
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

  character; one of `"delta"`, `"bootstrap"`, or `"none"`. The
  analytical delta method for EL is currently not implemented; when
  `"delta"` is supplied it is coerced to `"none"` with a warning.

- bootstrap_reps:

  integer; number of bootstrap replicates when
  `variance_method = "bootstrap"`.

- auxiliary_means:

  named numeric vector; population means for auxiliary design columns.
  Names must match the materialized model.matrix column names on the
  first RHS (after formula expansion), e.g., factor indicators like
  \`F_b\` or transformed terms \`I(X^2)\`. Intercept is always excluded.
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

  Unknown names are ignored. The method is Newton with an analytic
  Jacobian.

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

An engine object of class `c("nmar_engine_el","nmar_engine")`. This is a
configuration list; it is not a fit. Pass it to
[nmar](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).

## Details

Empirical likelihood assigns masses \\m_i = d_i / D_i(\theta)\\ to
respondents with base weights \\d_i\\. Following Qin, Leung, and Shao
(2002, JASA 97:193-200), the denominator is \$\$D_i(\theta) = 1 +
\lambda_W\\w_i(\beta) - W\\ + X\_{i\cdot}^{(c)}\\\lambda_x,\$\$ where
\\w_i(\beta) = g(\eta_i)\\ with link-inverse \\g\\, \\\eta_i =
x_i^\top\beta\\, \\W\\ is the (unknown) response rate, and \\X^{(c)} =
X - \mu_X\\ centers auxiliary columns at their population means. The
multiplier for the response-rate equation is \$\$\lambda_W =
\\N\_\mathrm{pop}/\sum_i d_i - 1\\/(1 - W),\$\$ which ensures scale
coherence for both IID data (\\d_i \equiv 1\\) and survey designs
(\\d_i\\ are design weights). The estimating equations impose \\\sum_i
m_i\\w_i(\beta) - W\\ = 0\\, \\\sum_i m_i\\s_i(\beta) = 0\\ with
\\s_i(\beta) = \partial \log w_i / \partial \eta_i\\, and, when present,
\\\sum_i m_i X\_{i\cdot}^{(c)} = 0\\.

The missingness-model score used in both equations and Jacobian is the
derivative of the Bernoulli log-likelihood with respect to the linear
predictor, i.e. `mu.eta(eta) / linkinv(eta)` (logit: `1 - w`; probit:
Mills ratio `phi/Phi`). We apply a consistent guarding policy (cap
\\\eta\\, clip \\w\\, floor denominators with an "active" mask in the
Jacobian) to ensure numerical stability and to make the analytic
Jacobian match the piecewise-smooth equations being solved. When
`variance_method = "delta"` is requested, the estimator returns `NA`
standard errors with a message; use `variance_method = "bootstrap"` for
SEs.

Survey designs are handled by replacing counts in QLS (2002) with
design-weighted totals. In particular, the response-rate multiplier
generalizes to \\\lambda_W = \\N\_\mathrm{pop}/\sum\_{i:\\R_i=1} d_i -
1\\/(1 - W)\\, which reduces to the QLS expression when \\d_i \equiv
1\\. Solver configuration uses `nleqslv` with an analytic Jacobian and
line-search globalization. Defaults are `global = "qline"` and
`xscalm = "auto"`; users can override via `control`. Invalid values are
coerced to these defaults with a warning.

**Formula syntax**:
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
supports a partitioned right-hand side `y_miss ~ aux1 + aux2 | z1 + z2`.
Variables left of `|` are auxiliaries (used in EL moment constraints);
variables right of `|` are missingness-model predictors only. The
outcome appears on the left-hand side and is included as a response
predictor by default.

**Weights in results**: Calling
[`weights()`](https://rdrr.io/r/stats/weights.html) on the returned
`nmar_result` gives respondent weights on either the probability scale
(sum to 1) or the population scale (sum to \\N\_\mathrm{pop}\\). The
reported masses come from the empirical likelihood construction
\\a_i/D_i(\theta)\\ and are normalized in
[`weights()`](https://rdrr.io/r/stats/weights.html).

**Variance**: Analytical delta variance for EL is not implemented.
Requesting `variance_method = "delta"` is coerced to `"none"` with a
warning. For standard errors in both IID and survey settings, use
`variance_method = "bootstrap"`.

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

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. *Journal of the
American Statistical Association*, 97(457), 193-200.

Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using
complete auxiliary information from survey data. *Journal of the
American Statistical Association*, 96(453), 185-193. Related to
design-based calibration; our EL approach balances auxiliary moments
through empirical likelihood constraints rather than calibration
adjustments to weights.

## See also

[nmar](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md);
see the vignette "Empirical Likelihood Theory for NMAR" for derivations.

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
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) -2.535301         NA      NA       NA
#> Y_miss       0.987603         NA      NA       NA

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
#>              Estimate Std. Error t value Pr(>|t|)
#> (Intercept) -2.535301         NA      NA       NA
#> Y_miss       0.987603         NA      NA       NA
# }
```
