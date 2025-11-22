# Trim weights by capping and proportional redistribution

Applies a cap to a nonnegative weight vector and, when feasible,
redistributes excess mass across the remaining positive entries so that
the total sum is preserved. When the requested cap is too tight to
preserve the total mass, all positive entries are set to the cap and the
total sum decreases.

## Usage

``` r
trim_weights(weights, cap, tol = 1e-12, warn_tol = 1e-08)
```

## Arguments

- weights:

  numeric vector of weights.

- cap:

  positive numeric scalar; maximum allowed weight, or `Inf` to disable
  trimming.

- tol:

  numeric tolerance used when testing whether a rescaling step respects
  the cap.

- warn_tol:

  numeric tolerance used when testing whether the total sum has been
  preserved.

## Value

A list with components:

- `weights`:

  numeric vector of trimmed weights.

- `trimmed_fraction`:

  fraction of entries at or very close to the cap (within `tol`).

- `preserved_sum`:

  logical; `TRUE` if the total sum of weights is preserved to within
  `warn_tol`.

- `total_before`:

  numeric; sum of the original weights.

- `total_after`:

  numeric; sum of the trimmed weights.

## Details

Internally, a simple water-filling style algorithm is used on the
positive weights: the largest weights are successively saturated at the
cap and the remaining weights are rescaled by a common factor chosen to
maintain the total sum.
