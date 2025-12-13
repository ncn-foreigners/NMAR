# Empirical Likelihood

## Overview

This vignette demonstrates the empirical likelihood (EL) estimator for
Not Missing At Random (NMAR) data in the NMAR package. The primary
estimand is the full-data mean of the outcome $Y$ under the QLS NMAR
model. The method implements the estimator of Qin, Leung, and Shao
(2002), using empirical likelihood weights that satisfy estimating
equations for the response mechanism and (optionally) auxiliary moment
constraints. For full derivations, the analytic Jacobian, and variance
discussion (bootstrap), see the companion article “Empirical Likelihood
Theory for NMAR”.

Key features:

- Supports `data.frame` (IID) and `survey.design` objects via the same
  [`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
  API.
- Variance via bootstrap (IID resampling or survey replicate weights).
- Optional standardization of predictors; weight trimming for
  robustness.
- Rich S3 surface: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), `tidy()`,
  `glance()`, [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

### Quick start

- Specify the model with a two-sided formula:
  `Y_miss ~ X1 + X2 | Z1 + Z2` where auxiliaries are left of `|` and
  response-only predictors are right of `|`. If `|` is omitted, RHS
  variables are treated only as auxiliaries; the response model uses the
  outcome (LHS) and any predictors explicitly placed to the right of
  `|`.
  - Variables on the outcome RHS (e.g., `X1 + X2`) are auxiliaries;
    supply their known population means via
    `auxiliary_means = c(X1 = ..., X2 = ...)`.
  - Predictors to the right of `|` enter only the response model (no
    auxiliary constraint) and do not need population means.
  - If you want a covariate to enter both the auxiliary constraints and
    the response model (as in the original QLS setup), include it on
    both sides, for example: `Y_miss ~ X | X`.
- Choose the engine: `el_engine(...)`, e.g.,
  `el_engine(auxiliary_means = c(X1 = 0), variance_method = "bootstrap", standardize = TRUE)`.
- Fit:
  `nmar(formula = Y_miss ~ X1 + X2 | Z1 + Z2, data = df_or_design, engine = el_engine(...))`.
- Inspect: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
  `fit$diagnostics`.

## Data-frame example (IID)

We simulate an NMAR mechanism where the response probability depends on
the unobserved outcome.

``` r
set.seed(123)
library(NMAR)

N <- 500
X <- rnorm(N)
Z <- rnorm(N)
Y <- 2 + 0.5 * X + Z

# NMAR response: depends on Y
p <- plogis(-1.0 + 0.4 * scale(Y)[, 1])
R <- runif(N) < p

dat <- data.frame(Y_miss = Y, X = X)
dat$Y_miss[!R] <- NA_real_
```

``` r
engine <- el_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE)
# Fit EL estimator (no variance for speed in vignette)
fit <- nmar(
  formula = Y_miss ~ X,
  data = dat,
  engine = engine
)

summary(fit)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.878660
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 500 
#> Respondents: 150 
#> Call: nmar(Y_miss ~ X, data = <data.frame: N=500>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -1.570694
#> Y_miss       0.366754
# For confidence intervals, use bootstrap variance (see example below).
```

Probit family (optional):

``` r
engine <- el_engine(auxiliary_means = c(X = 0), family = "probit", variance_method = "none", standardize = TRUE)

fit_probit <- nmar(
  formula = Y_miss ~ X,
  engine = engine,
  data = dat
)
summary(fit_probit)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.880128
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 500 
#> Respondents: 150 
#> Call: nmar(Y_miss ~ X, data = <data.frame: N=500>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -0.949678
#> Y_miss       0.217504
```

Tidy/glance summaries (via the `generics` package):

``` r
generics::tidy(fit)
#>          term   estimate std.error conf.low conf.high component statistic
#> 1      Y_miss  1.8786601        NA       NA        NA  estimand        NA
#> 2 (Intercept) -1.5706942        NA       NA        NA  response        NA
#> 3      Y_miss  0.3667545        NA       NA        NA  response        NA
#>   p.value
#> 1      NA
#> 2      NA
#> 3      NA
generics::glance(fit)
#>     y_hat std.error conf.low conf.high converged trimmed_fraction
#> 1 1.87866        NA       NA        NA      TRUE                0
#>   variance_method jacobian_condition_number max_equation_residual
#> 1            none                  36.83019           5.17808e-13
#>   min_denominator fraction_small_denominators nobs nobs_resp is_survey
#> 1       0.4613402                           0  500       150     FALSE
```

Outputs and diagnostics at a glance (probability-scale weights sum to 1;
population-scale weights sum to the analysis total $N_{\text{pop}}$):

``` r
head(weights(fit), 10)
#>  [1] 0.008384402 0.008937803 0.004381958 0.014450653 0.006959120 0.007160884
#>  [7] 0.005995129 0.006358571 0.007390444 0.007842733
head(weights(fit, scale = "population"), 10)
#>  [1] 4.192201 4.468901 2.190979 7.225327 3.479560 3.580442 2.997564 3.179285
#>  [9] 3.695222 3.921366
head(fitted(fit), 10)
#>  [1] 0.2385382 0.2237686 0.4564170 0.1384020 0.2873926 0.2792951 0.3336042
#>  [8] 0.3145361 0.2706197 0.2550132
str(fit$diagnostics)
#> List of 35
#>  $ convergence_code                  : int 1
#>  $ message                           : chr "Function criterion near zero"
#>  $ vcov_message                      : chr "Variance skipped (variance_method='none')"
#>  $ trimmed_fraction                  : num 0
#>  $ solver_method                     : chr "Newton"
#>  $ nleqslv_global                    : chr "qline"
#>  $ nleqslv_xscalm                    : chr "auto"
#>  $ solver_iterations                 : int 7
#>  $ solver_time                       : num 0.004
#>  $ variance_time                     : num 0
#>  $ reparam_W                         : chr "logit"
#>  $ max_equation_residual             : num 5.18e-13
#>  $ jacobian_condition_number         : num 36.8
#>  $ auxiliary_inconsistency_max_z     : num 0.114
#>  $ auxiliary_inconsistency_cols      : chr "X"
#>  $ min_denominator                   : num 0.461
#>  $ fraction_small_denominators       : num 0
#>  $ denom_q01                         : num 0.496
#>  $ denom_q05                         : num 0.664
#>  $ denom_median                      : num 1.05
#>  $ denom_count_lt_1e4                : int 0
#>  $ denom_floor                       : num 1e-08
#>  $ denom_floor_hits                  : num 0
#>  $ weight_max_share                  : num 0.0145
#>  $ weight_top5_share                 : num 0.0665
#>  $ weight_ess                        : num 137
#>  $ constraint_sum_W                  : num -3.33e-13
#>  $ constraint_sum_aux                : Named num -5.18e-13
#>   ..- attr(*, "names")= chr "X"
#>  $ constraint_sum_link               : num NA
#>  $ sum_respondent_weights            : num 150
#>  $ sum_unnormalized_weights_untrimmed: num 150
#>  $ normalization_ratio               : num 1
#>  $ max_constraint_residual           : num 5.18e-13
#>  $ auxiliary_means                   : Named num 0
#>   ..- attr(*, "names")= chr "X"
#>  $ auxiliary_matrix                  : num [1:150, 1] -0.56 -0.23 1.559 -1.265 -0.687 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:150] "1" "2" "3" "8" ...
#>   .. ..$ : chr "X"
```

Bootstrap variance (keep reps small for speed). This example requires
the optional future.apply package; the chunk is skipped if it is not
installed:

``` r
set.seed(123)
engine <- el_engine(
  auxiliary_means = c(X = 0),
  variance_method = "bootstrap",
  bootstrap_reps = 10,
  standardize = TRUE
)
fit_boot <- nmar(
  formula = Y_miss ~ X,
  engine = engine,
  data = dat
)
se(fit_boot)
#> [1] 0.2817672
confint(fit_boot)
#>           2.5 %   97.5 %
#> Y_miss 1.326407 2.430914
```

## Respondents-only data (n_total)

If you pass respondents-only data (the outcome contains no NA), provide
the total sample size via `n_total` in the engine so the estimator can
recover the population response rate:

``` r
set.seed(124)
N <- 300
X <- rnorm(N); Z <- rnorm(N); Y <- 1.5 + 0.4 * X + Z
p <- plogis(-0.5 + 0.4 * scale(Y)[, 1])
R <- runif(N) < p
df_resp <- subset(data.frame(Y_miss = Y, X = X), R == 1)
eng_resp <- el_engine(auxiliary_means = c(X = 0), variance_method = "none", n_total = N)
fit_resp <- nmar(Y_miss ~ X, data = df_resp, engine = eng_resp)
summary(fit_resp)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.509469
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 300 
#> Respondents: 102 
#> Call: nmar(Y_miss ~ X, data = <data.frame: N=300>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -0.886331
#> Y_miss       0.145580
```

## Response-only predictors

You can include predictors that enter only the response model (and are
not constrained as auxiliaries) by placing them to the right of `|` in
the formula.

``` r
set.seed(125)
N <- 400
X <- rnorm(N)
Z <- rnorm(N)
Y <- 1 + 0.6 * X + 0.3 * Z + rnorm(N)
p <- plogis(-0.6 + 0.5 * scale(Y)[, 1] + 0.4 * Z)
R <- runif(N) < p
df2 <- data.frame(Y_miss = Y, X = X, Z = Z)
df2$Y_miss[!R] <- NA_real_
engine <- el_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE)

# Use X as auxiliary (known population mean 0), and Z as response-only predictor
fit_resp_only <- nmar(
  formula = Y_miss ~ X | Z,
  data = df2,
  engine = engine
)
summary(fit_resp_only)
#> NMAR Model Summary
#> =================
#> Y_miss mean: 1.136991
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 400 
#> Respondents: 139 
#> Call: nmar(Y_miss ~ X | Z, data = <data.frame: N=400>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -1.019139
#> Y_miss       0.369372
#> Z           -0.232795
```

Auxiliary means and formulas:

- Names of `auxiliary_means` must match the `model.matrix` columns
  generated by the outcome RHS (for numeric variables this typically
  equals the variable names; for factors it corresponds to the indicator
  columns).
- When `standardize = TRUE`, the engine automatically transforms
  `auxiliary_means` to the standardized scale internally and reports
  coefficients on the original scale.
- Response-only predictors (to the right of `|`) do not need auxiliary
  means.

## Survey design example (optional)

The estimator supports complex surveys via
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).
This chunk runs only if the `survey` package is available. If the survey
weights were normalized (for example, rescaled to sum to the sample
size), pass an appropriate population total via
`el_engine(n_total = ...)` so that population-scale EL weights are on
the intended scale.

``` r
suppressPackageStartupMessages(library(survey))
data(api)

set.seed(42)
apiclus1$api00_miss <- apiclus1$api00
ystd <- scale(apiclus1$api00)[, 1]
prob <- plogis(-0.5 + 0.4 * ystd + 0.2 * scale(apiclus1$ell)[, 1])
miss <- runif(nrow(apiclus1)) > prob
apiclus1$api00_miss[miss] <- NA_real_

dclus1 <- svydesign(id = ~dnum, weights = ~pw, data = apiclus1, fpc = ~fpc)
# Let the engine infer auxiliary means from the full design (design-weighted).
# Alternatively, you can supply known population means via auxiliary_means.
engine <- el_engine(auxiliary_means = NULL, variance_method = "none", standardize = TRUE)

fit_svy <- nmar(
  formula = api00_miss ~ ell | ell,
  data = dclus1,
  engine = engine
)
summary(fit_svy)
#> NMAR Model Summary
#> =================
#> api00_miss mean: 667.708681
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 6194 
#> Respondents: 65 
#> Call: nmar(api00_miss ~ ell | ell, data = <survey.design: N=6194>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept) -1.246046
#> api00_miss   0.000168
#> ell          0.019037
```

## Practical guidance

- Variance method: Analytical delta variance for EL is currently
  disabled; use bootstrap for standard errors.
- Trimming: Use a finite `trim_cap` to improve robustness when large
  weights occur; prefer bootstrap variance when trimming.
- Solver control: set
  `control = list(xtol = 1e-10, ftol = 1e-10, maxit = 200)` for tighter
  tolerances if needed. Globalization details are managed internally by
  `nleqslv`.
- Standardization: `standardize = TRUE` typically improves numerical
  stability and comparability across predictors and auxiliary means.
- Diagnostics: Inspect `fit$diagnostics` (Jacobian condition number, max
  equation residuals, trimming fraction) to assess numerical health and
  identification strength.
- Response-only predictors: Variables to the right of `|` do not need to
  appear on the RHS of the outcome formula; they enter only the response
  model. Auxiliary means must be supplied only for variables on the
  outcome RHS.
- Inconsistent auxiliaries: If provided auxiliary means are grossly
  inconsistent with the respondent support, the engine will issue a
  warning via auxiliary inconsistency diagnostics and the solver may
  fail or yield highly concentrated weights. Consider revisiting the
  constraints, relaxing them, or using `trim_cap` and bootstrap
  variance.

Troubleshooting:

- Extreme weights: set a finite `trim_cap`; prefer
  `variance_method = "bootstrap"` for SE.
- Ill-conditioned Jacobian (large
  `fit$diagnostics$jacobian_condition_number`): prefer
  `variance_method = "bootstrap"`. You may also tighten solver
  tolerances via `control = list(xtol=..., ftol=..., maxit=...)`.
- Convergence issues: check `fit$diagnostics$max_equation_residual`,
  rescale predictors (`standardize = TRUE`), or reduce the number of
  constraints.

## Solver control and notes

- Control example: increase iterations and tighten tolerances via
  `control` (passed to `nleqslv`):

``` r
ctrl <- list(maxit = 200, xtol = 1e-10, ftol = 1e-10)
eng_ctrl <- el_engine(auxiliary_means = c(X = 0), variance_method = "none", control = ctrl)
invisible(nmar(Y_miss ~ X, data = dat, engine = eng_ctrl))
```

## References and further reading

- Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
  under nonignorable nonresponse or informative sampling. Journal of the
  American Statistical Association, 97(457), 193-200.
  <doi:10.1198/016214502753479338>
- Chen, J. and Sitter, R. R. (1999). A pseudo empirical likelihood
  approach to the effective use of auxiliary information in complex
  surveys. Statistica Sinica, 9, 385-406.
- Wu, C. (2005). Algorithms and R codes for the pseudo empirical
  likelihood method in survey sampling. Survey Methodology, 31(2),
  239-243.

## Families and numerical stability

- Family: `el_engine(family = "logit")` (default) or
  `family = "probit"`.
- Probit stability: the response-model score is evaluated as the Mills
  ratio $\phi(\eta)/\Phi(\eta)$ in the log domain for numerical
  stability; the logit score simplifies to $1 - {plogis}(\eta)$. We also
  cap the linear predictor and clip probabilities used in ratios.
- Theory mapping: see the companion article “Empirical Likelihood Theory
  for NMAR” for equations, Jacobian blocks, and variance details.

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] survey_4.4-8   survival_3.8-3 Matrix_1.7-4   future_1.68.0  NMAR_0.0.0.1  
#> 
#> loaded via a namespace (and not attached):
#>  [1] future.apply_1.20.1 jsonlite_2.0.0      compiler_4.5.2     
#>  [4] Rcpp_1.1.0          nleqslv_3.3.5       parallel_4.5.2     
#>  [7] jquerylib_0.1.4     globals_0.18.0      splines_4.5.2      
#> [10] systemfonts_1.3.1   textshaping_1.0.4   yaml_2.3.12        
#> [13] fastmap_1.2.0       lattice_0.22-7      R6_2.6.1           
#> [16] generics_0.1.4      Formula_1.2-5       knitr_1.50         
#> [19] htmlwidgets_1.6.4   desc_1.4.3          DBI_1.2.3          
#> [22] bslib_0.9.0         rlang_1.1.6         cachem_1.1.0       
#> [25] xfun_0.54           fs_1.6.6            sass_0.4.10        
#> [28] cli_3.6.5           progressr_0.18.0    pkgdown_2.2.0      
#> [31] digest_0.6.39       lifecycle_1.0.4     evaluate_1.0.5     
#> [34] listenv_0.10.0      codetools_0.2-20    ragg_1.5.0         
#> [37] mitools_2.4         parallelly_1.46.0   rmarkdown_2.30     
#> [40] tools_4.5.2         htmltools_0.5.9
```

## Notes on variance choices

- Analytical delta variance for EL is not implemented in this version.
  Use `variance_method = "bootstrap"` for SEs.
- For speed-critical bootstraps, you can set `variance_method = "none"`
  for point fits and let the bootstrap call use it internally for
  replicates (this is the default behavior when you request bootstrap
  SEs through `el_engine(...)`).
