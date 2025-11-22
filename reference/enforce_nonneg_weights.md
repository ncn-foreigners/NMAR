# Enforce (near-)nonnegativity of weights

Softly enforces nonnegativity of a numeric weight vector. Large negative
values (beyond a tolerance) are treated as errors; small negative values
(for example, from numerical noise) are truncated to zero.

## Usage

``` r
enforce_nonneg_weights(weights, tol = 1e-08)
```

## Arguments

- weights:

  numeric vector of weights.

- tol:

  numeric tolerance below which negative values are treated as numerical
  noise and clipped to zero.

## Value

A list with components:

- `ok`:

  logical; `TRUE` if no clearly negative weights were found.

- `message`:

  character; diagnostic message when `ok` is `FALSE`, otherwise `NULL`.

- `weights`:

  numeric vector of adjusted weights (original if `ok` is `FALSE`,
  otherwise with small negatives clipped to zero).
