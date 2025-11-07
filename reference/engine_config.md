# Extract engine configuration

Returns the underlying configuration of an engine as a named list. This
is intended for programmatic inspection (e.g., parameter tuning,
logging). The returned object should be treated as read-only.

## Usage

``` r
engine_config(x)
```

## Arguments

- x:

  An object inheriting from class \`nmar_engine\`.

## Value

A named list of configuration fields.
