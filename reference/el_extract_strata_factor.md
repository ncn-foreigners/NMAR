# Extract a strata factor from a survey.design object

Prefers strata already materialized in the `survey.design` object
(typically `design$strata`). When unavailable, attempts to reconstruct
strata from the original
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html) call.
When multiple stratification variables are supplied, their interaction
is used.

## Usage

``` r
el_extract_strata_factor(design)
```
