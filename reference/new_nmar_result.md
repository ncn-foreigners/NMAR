# Construct Result Object (parent helper)

Builds an \`nmar_result\` list using the shared schema and validates it.
Engines must pass named fields; no legacy positional signature is
supported.

## Usage

``` r
new_nmar_result(...)
```

## Details

Engine-level constructors should call this helper with named arguments
rather than assembling result lists by hand. At minimum, engines should
supply `estimate` (numeric scalar) and `converged` (logical). All other
fields are optional:

- `estimate_name`: label for the primary estimand (defaults to
  `NA_character_` if omitted).

- `se`: standard error for the primary estimand (defaults to `NA_real_`
  when not available).

- `model`, `weights_info`, `sample`, `inference`, `diagnostics`, `meta`,
  `extra`: lists that may be partially specified or `NULL`;
  [`validate_nmar_result()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/validate_nmar_result.md)
  will back-fill missing subfields with safe defaults.

- `class`: engine-specific result subclass name, e.g.
  `"nmar_result_el"`; it is combined with the parent class
  `"nmar_result"`.

Calling `new_nmar_result()` ensures that every engine returns objects
that satisfy the shared schema and are immediately compatible with
parent S3 methods such as [`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html), `tidy()`,
`glance()`, and [`weights()`](https://rdrr.io/r/stats/weights.html).
