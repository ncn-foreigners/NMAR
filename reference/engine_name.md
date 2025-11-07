# Canonical engine name

Returns a stable, machine-friendly identifier for an engine object. This
identifier is also used in \`nmar_result\$meta\$engine_name\` to keep a
consistent naming scheme between configurations and results.

## Usage

``` r
engine_name(x)
```

## Arguments

- x:

  An object inheriting from class \`nmar_engine\`.

## Value

A single character string, e.g. "empirical_likelihood".
