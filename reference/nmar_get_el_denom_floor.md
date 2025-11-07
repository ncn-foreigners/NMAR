# EL denominator floor (global, consistent)

Returns the small positive floor \\\delta\\ used to guard the empirical
likelihood denominator \\D_i(\theta)\\ away from zero. This guard must
be applied consistently in the estimating equations, analytic Jacobian,
and post-solution weight construction. Advanced users can override via
\`options(nmar.el_denom_floor = 1e-8)\`.

## Usage

``` r
nmar_get_el_denom_floor()
```
