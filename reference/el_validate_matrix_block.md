# Validate a design-matrix block and optionally return respondent rows

Applies the respondent mask, enforces the "no NA, no zero-variance"
policy, and returns the subset that downstream routines should use.

## Usage

``` r
el_validate_matrix_block(
  mat,
  mask,
  row_map,
  label,
  severity = c("error", "warn")
)
```
