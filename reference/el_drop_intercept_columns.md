# Remove \`(Intercept)\` columns from a model matrix

Auxiliary constraints should never include an intercept, and the
missingness design injects its own intercept. We therefore strip the
automatically generated \`(Intercept)\` column before running downstream
validation.

## Usage

``` r
el_drop_intercept_columns(mat)
```
