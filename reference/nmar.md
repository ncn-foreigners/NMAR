# Not Missing at Random (NMAR) Estimation

High-level interface for NMAR estimation.

`nmar()` validates basic inputs and dispatches to an engine (for example
[`el_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md)).
The engine controls the estimation method and interprets `formula`; see
the engine documentation for model-specific requirements.

## Usage

``` r
nmar(formula, data, engine, trace_level = 0)
```

## Arguments

- formula:

  A two-sided formula. Many engines support a partitioned right-hand
  side via `|`, for example `y_miss ~ block1_vars | block2_vars`. The
  meaning of these blocks is engine-specific (see the engine
  documentation). In the common "missing values indicate nonresponse"
  workflow, the left-hand side is the outcome with `NA` values for
  nonrespondents.

- data:

  A `data.frame` or a `survey.design` containing the variables
  referenced by `formula`.

- engine:

  An NMAR engine configuration object, typically created by
  [`el_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md),
  [`exptilt_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/exptilt_engine.md),
  or
  [`exptilt_nonparam_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/exptilt_nonparam_engine.md).
  This object defines the estimation method and tuning parameters and
  must inherit from class `"nmar_engine"`.

- trace_level:

  Integer 0-3; controls verbosity during estimation (default `0`):

  - 0: no output (silent mode);

  - 1: major steps only (initialization, convergence, final results);

  - 2: iteration summaries and key diagnostics;

  - 3: full diagnostic output.

## Value

An object of class `"nmar_result"` with an engine-specific subclass (for
example `"nmar_result_el"`). Use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`se`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/se.md),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`weights()`](https://rdrr.io/r/stats/weights.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
to access estimates, standard errors, weights, and diagnostics.

## See also

[`el_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md),
[`exptilt_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/exptilt_engine.md),
[`exptilt_nonparam_engine`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/exptilt_nonparam_engine.md),
[`summary.nmar_result`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/summary.nmar_result.md),
[`weights.nmar_result`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/weights.nmar_result.md)

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n)
z1 <- rnorm(n)
y_true <- 0.5 + 0.3 * x1 + 0.2 * z1 + rnorm(n, sd = 0.3)
resp <- rbinom(n, 1, plogis(2 + 0.1 * y_true + 0.1 * z1))
if (all(resp == 1)) resp[sample.int(n, 1)] <- 0L
y_obs <- ifelse(resp == 1, y_true, NA_real_)

# Empirical likelihood engine
df_el <- data.frame(Y_miss = y_obs, X = x1, Z = z1)
eng_el <- el_engine(variance_method = "none")
fit_el <- nmar(Y_miss ~ X | Z, data = df_el, engine = eng_el)
summary(fit_el)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 0.527456
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 200 
#> Respondents: 179 
#> Call: nmar(Y_miss ~ X | Z, data = <data.frame: N=200>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept)  2.460382
#> Y_miss       0.118914
#> Z           -0.875952

# \donttest{
# Exponential tilting engine (illustrative)
dat_et <- data.frame(y = y_obs, x2 = z1, x1 = x1)
eng_et <- exptilt_engine(
  y_dens = "normal",
  family = "logit",
  variance_method = "none"
)
fit_et <- nmar(y ~ x2 | x1, data = dat_et, engine = eng_et)
summary(fit_et)
#> NMAR Model Summary (Exponential tilting)
#> =================================
#> y mean: 0.515218
#> Converged: TRUE 
#> Variance method: none 
#> Call: nmar(y ~ x2 | x1, data = <data.frame: N=?>, engine = exponential_tilting)
#> 
#> Response-model (theta) coefficients:
#>   (Intercept)          : 2.214475
#>   x1                   : 0.120835
#>   y                    : -0.161063

# Survey design example (same outcome, random weights)
if (requireNamespace("survey", quietly = TRUE)) {
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w,
                           data = data.frame(Y_miss = y_obs, X = x1, Z = z1))
  eng_svy <- el_engine(variance_method = "none")
  fit_svy <- nmar(Y_miss ~ X | Z, data = des, engine = eng_svy)
  summary(fit_svy)
}
#> NMAR Model Summary
#> =================
#> Y_miss mean: 0.516401
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 253.0348 
#> Respondents: 179 
#> Call: nmar(Y_miss ~ X | Z, data = <survey.design: N=253.035>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept)  1.938651
#> Y_miss       0.265915
#> Z           -0.046734

# Bootstrap variance usage
if (requireNamespace("future.apply", quietly = TRUE)) {
  set.seed(2)
  eng_boot <- el_engine(
    variance_method = "bootstrap",
    bootstrap_reps = 20
  )
  fit_boot <- nmar(Y_miss ~ X | Z, data = df_el, engine = eng_boot)
  se(fit_boot)
}
#> [1] 0.0298057
# }
```
