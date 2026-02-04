# Weighted linear algebra

Compute X' diag(w) X efficiently. If w \>= 0, use SPD
crossprod(X\*sqrt(w)). Otherwise, fall back to X' (diag(w) X) via
crossprod(X, X\*w).

## Usage

``` r
shared_weighted_gram(X, w)
```
