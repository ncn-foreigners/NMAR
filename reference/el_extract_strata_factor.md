# Extract a strata factor from a survey.design object

Uses the original svydesign() call stored in the object to recreate the
stratum labels as a single factor. When multiple stratification
variables are supplied, their interaction is used.

## Usage

``` r
el_extract_strata_factor(design)
```
