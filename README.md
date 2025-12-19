
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NMAR <img src="logo.png" align="right" height="120" alt="NMAR" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/ncn-foreigners/NMAR/graph/badge.svg)](https://app.codecov.io/gh/ncn-foreigners/NMAR)
[![R-CMD-check](https://github.com/ncn-foreigners/NMAR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncn-foreigners/NMAR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

NMAR provides estimators for finite-population means when outcomes are
subject to nonignorable nonresponse (Not Missing at Random, NMAR). It
supports iid `data.frame` inputs and complex survey designs via
`survey::survey.design`, and exposes a unified interface through
`nmar()`.

## Methods

NMAR currently provides the following engines:

- `el_engine()`: empirical likelihood (Qin, Leung and Shao, 2002).
- `exptilt_engine()`: exponential tilting (Riddles, Kim and Im, 2016).
- `exptilt_nonparam_engine()`: nonparametric exponential tilting for
  aggregated categorical data (Riddles, Kim and Im, 2016, Appendix 2).

References:

- Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
  under nonignorable nonresponse or informative sampling. Journal of the
  American Statistical Association, 97(457), 193-200.
  <https://doi.org/10.1198/016214502753479338>
- Riddles, M. K., Kim, J. K., and Im, J. (2016). A
  propensity-score-adjustment method for nonignorable nonresponse.
  Journal of Survey Statistics and Methodology, 4(2), 215-245.
  <https://doi.org/10.1093/jssam/smv047>

See `browseVignettes("NMAR")` and the package website for worked
examples and engine-specific assumptions:
<https://ncn-foreigners.ue.poznan.pl/NMAR/>.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("ncn-foreigners/NMAR")
```

Or with remotes:

``` r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/NMAR")
```

To build vignettes locally:

``` r
remotes::install_github("ncn-foreigners/NMAR", build_vignettes = TRUE)
```

## Usage

### Formula interface

`nmar()` uses a two-sided formula with the outcome on the left-hand
side. In the common “missing values indicate nonresponse” workflow,
nonrespondents are encoded as `NA` in the outcome.

Many engines support a partitioned right-hand side via `|` (e.g.,
`y_miss ~ block1_vars | block2_vars`), but the interpretation of these
blocks is engine-specific. See `?el_engine` and `?exptilt_engine` for
details.

### Example

``` r
suppressPackageStartupMessages(library(NMAR))

data("riddles_case1", package = "NMAR")

# Empirical likelihood (EL)
fit_el <- nmar(
  y ~ x,
  data = riddles_case1,
  engine = el_engine(variance_method = "none")
)
summary(fit_el)
#> NMAR Model Summary
#> =================
#> y mean: -1.001986
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 500 
#> Respondents: 368 
#> Call: nmar(y ~ x, data = <data.frame: N=500>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept)  0.860555
#> y           -0.175376

# Exponential tilting (ET)
fit_et <- nmar(
  y ~ x,
  data = riddles_case1,
  engine = exptilt_engine(y_dens = "normal", family = "logit", variance_method = "none")
)
summary(fit_et)
#> NMAR Model Summary (Exponential tilting)
#> =================================
#> y mean: -1.003729
#> Converged: TRUE 
#> Variance method: none 
#> Call: nmar(y ~ x, data = <data.frame: N=?>, engine = exponential_tilting)
#> 
#> Response-model (theta) coefficients:
#>   (Intercept)          : 0.863720
#>   y                    : -0.170861
```

Result objects returned by `nmar()` support methods such as `summary()`,
`weights()`, `se()`, and `confint()`, and broom-style `tidy()` /
`glance()` via the `generics` package.

### Survey designs

``` r
if (requireNamespace("survey", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survey))
  set.seed(1)
  d <- riddles_case1
  d$w <- runif(nrow(d), 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = d)

  fit_svy <- nmar(y ~ x, data = des, engine = el_engine(variance_method = "none"))
  summary(fit_svy)
}
#> NMAR Model Summary
#> =================
#> y mean: -1.005961
#> Converged: TRUE 
#> Variance method: none 
#> Variance notes: Variance skipped (variance_method='none') 
#> Total units: 621.7412 
#> Respondents: 368 
#> Call: nmar(y ~ x, data = <survey.design: N=621.741>, engine = empirical_likelihood)
#> 
#> Missingness-model coefficients:
#>              Estimate
#> (Intercept)  0.893453
#> y           -0.129256
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/ncn-foreigners/NMAR/issues)

## Citation

If you use NMAR in academic work, please cite the package and the
relevant method paper(s):

``` r
citation("NMAR")
```

## Authors and acknowledgments

Research grant: [OPUS 20
\#2020/39/B/HS4/00941](https://ncn-foreigners.ue.poznan.pl/)

- [Maciej Beręsewicz](https://github.com/berenz)
- [Igor Kołodziej](https://github.com/IgorKolodziej)
- [Mateusz Iwaniuk](https://github.com/Iwaniukooo11)
