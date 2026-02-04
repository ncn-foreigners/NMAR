# Extract strata factor

Looks for strata already materialized in the `survey.design` object.
When unavailable, attempts to reconstruct strata from the original
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html) call.
When multiple stratification variables are supplied, their interaction
is used.

## Usage

``` r
el_extract_strata_factor(design)
```
