# Changelog

## NMAR 0.1.2

- Bootstrap variance no longer hard-requires ‘future.apply’; it now
  falls back to sequential execution via base::lapply() when
  ‘future.apply’ is not installed (emits a one-time warning).
- When ‘future.apply’ is installed, bootstrap uses future-seeded RNG
  streams (future.seed = TRUE) for reproducibility across future
  backends.

## NMAR 0.1.1

CRAN release: 2026-01-16

- CRAN release-related fixes
- Fix `return` roxygen keyword in S3 Functions
- Add research doi references to DESCRIPTION file

## NMAR 0.1.0

### Initial CRAN Release

- First release of the **NMAR** package for estimating nonignorable
  nonresponse (NMAR) bias in survey data.

### Methods

- **Empirical Likelihood (EL):** Added
  [`el_engine()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md)
  implementing the estimator of Qin, Leung, and Shao (2002). This method
  uses empirical likelihood weights satisfying response mechanism
  equations and auxiliary moment constraints.
- **Exponential Tilting (Parametric & Nonparametric):** Included robust
  implementations for both microdata (`exptilt_engine`) and aggregated
  contingency tables (`exptilt_nonparam_engine`) based on Riddles, Kim,
  and Im (2016).

### Key Features

- **Unified API:** All estimators are accessible via a single,
  consistent
  [`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
  interface supporting standard formula syntax (e.g., `Y ~ X | Z`).
- **Complex Survey Support:** Seamless integration with the `survey`
  package.
  [`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
  accepts `survey.design` objects, automatically handling weights and
  stratification.
- **Variance Estimation:** Robust bootstrapping (S3) implementation for
  standard errors and confidence intervals across all engines.
- **Diagnostics:** Rich return objects including convergence statistics,
  Jacobian condition numbers, and detailed weight summaries.

### Major Changes

- **Refactored Architecture:** The `exptilt` and `el` engines share a
  unified structural design, ensuring consistent behavior for controls,
  standardization, and error handling.
- **Standardization:** Added `standardize = TRUE` argument to engines to
  improve numerical stability during optimization.
