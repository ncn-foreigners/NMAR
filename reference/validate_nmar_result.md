# Validate nmar_result

Ensures both the child class and the parent schema are satisfied. The
validator also back-fills defaults so downstream code can rely on the
presence of optional components without defensive checks.

## Usage

``` r
validate_nmar_result(x, class_name)
```

## Details

This helper is the single authority on the \`nmar_result\` schema. It
expects a list that already carries class `c(class_name, "nmar_result")`
and at least a primary estimate stored in `y_hat`. All other components
are optional. When they are `NULL` or missing, the validator supplies
safe defaults:

- Core scalars: `se` (numeric, default `NA_real_`), `estimate_name`
  (character, default `NA_character_`), `converged` (logical, default
  `NA`).

- `model`: list with `coefficients` and `vcov`, both defaulting to
  `NULL`.

- `weights_info`: list with `values` (default `NULL`) and
  `trimmed_fraction` (default `NA_real_`).

- `sample`: list with `n_total`, `n_respondents`, `is_survey`, and
  `design`, defaulted to missing/empty values.

- `inference`: list with `variance_method`, `df`, and `message`, all
  defaulted to missing values.

- `diagnostics`, `meta`, and `extra`: defaulted to empty lists, with
  `meta` carrying `engine_name`, `call`, and `formula` when unset.

Engine constructors should normally call
[`new_nmar_result()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/new_nmar_result.md)
rather than invoking this function directly.
[`new_nmar_result()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/new_nmar_result.md)
attaches classes and funnels all objects through
`validate_nmar_result()` so downstream S3 methods can assume a
consistent structure.
