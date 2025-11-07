# Validate nmar_result structure

Ensures both the child class and the parent schema are satisfied. The
validator also back-fills defaults so downstream code can rely on the
presence of optional components without defensive checks.

## Usage

``` r
validate_nmar_result(x, class_name)
```
