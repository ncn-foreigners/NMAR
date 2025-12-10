# NMAR

The goal of **NMAR** is to provide a set of functions to **estimate the
population mean** of data subject to Not Missing At Random (NMAR)
mechanisms using advanced statistical methods

## Overview

The **NMAR** library provides functions to estimate the mean of **Not
Missing At Random (NMAR)** data. The estimation functions are built upon
the following engine functions, each implementing a distinct method:

- `exptilt_engine`: Exponential Tilting Estimator  
  Minsun Kim Riddles, Jae Kwang Kim, Jongho Im  
  A Propensity-score-adjustment Method for Nonignorable Nonresponse  
  <https://doi.org/10.1093/jssam/smv047>

- `exptilt_nonparam_engine`: Nonparametric Exponential Tilting
  Estimator  
  Minsun Kim Riddles, Jae Kwang Kim, Jongho Im  
  A Propensity-score-adjustment Method for Nonignorable Nonresponse  
  Appendix 2  
  <https://doi.org/10.1093/jssam/smv047>

- `el_engine`: Empirical Likelihood Estimator  
  Jing Qin, Denis Leung, Jun Shao  
  Estimation With Survey Data Under Nonignorable Nonresponse or
  Informative Sampling  
  <http://dx.doi.org/10.1198/016214502753479338>

The main user-facing function is
[`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md),
which acts as a unified wrapper around these engine functions.

## Installation

You can install the development version of NMAR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ncn-foreigners/NMAR")
```

or

``` R
remotes::install_github("ncn-foreigners/NMAR@main", force = T, build = T, build_manual = T, build_vignettes = T)
```

### Project branches:

- `main` - The stable production branch. This is the version you should
  install for general use
- `package-pre_prod` - The pre-production development branch. This is a
  stable, installable version reflecting the latest features. (1:1 clone
  of `package-dev`)
- `package-dev` - The **internal development branch**. This branch uses
  a non-R package structure for internal development purposes. **DO
  NOT** download or install this branch. For contribution guidelines,
  please see the CONTRIBUTING.md file or the \[project webpage (TODO:
  add link to webpage)\]

## Usage

### General naming

- `outcome_var` (f.e Y): The outcome variable (with missing values).
- `covariates_for_outcome`(f.e x1,x2): Predictors of `outcome_var` value
  (used in the response model).
- `covariates_for_missingness`(f.e x3): Predictor of `outcome_var`
  missingness (used in the missingness model).

`formula` = outcome_var ~ covariates_for_outcome

### Example

If, in a salary survey, richer people were less likely to answer, `Y` is
salary, `x1` and `x2` are experience and education, and `x3` might be
gender (if we assume gender affects the likelihood of responding but not
the salary value itself).

``` r
library(NMAR, quietly = T)

exptilt_config <- exptilt_engine(
  y_dens = 'normal',
  family = 'probit', # or logit
  variance_method = 'bootstrap', # or delta
  bootstrap_reps = 10
)

formula = y ~ x
res <- nmar(formula = formula, data = riddles_case1, engine = exptilt_config, trace_level = 0)
print(coef(res))
#> (Intercept)           y 
#>   0.5331082  -0.1024071
print(res)
#> NMAR Result
#> ------------
#> y mean: -1.003112 (0.056375)
#> Converged: TRUE 
#> Variance method: bootstrap 
#> Estimator: exponential_tilting
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/ncn-foreigners/NMAR/issues)

## Authors and acknowledgments

Research grant: [OPUS 20
\#2020/39/B/HS4/00941](https://ncn-foreigners.ue.poznan.pl/)

- [Maciej Beręsewicz](https://github.com/berenz)
- [Igor Kołodziej](https://github.com/IgorKolodziej)
- [Mateusz Iwaniuk](https://github.com/Iwaniukooo11)
