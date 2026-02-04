# NMAR numeric settings

NMAR numeric settings

## Usage

``` r
nmar_get_numeric_settings()
```

## Value

A named list with entries \`eta_cap\`, \`grad_eps\`, and \`grad_d\`.

## Details

Centralized access to numeric thresholds used across the package.

\- \`nmar.eta_cap\`: scalar \> 0. Caps the response-model linear
predictor to avoid extreme link values in Newton updates. Default 50. -
\`nmar.grad_eps\`: finite-difference step size epsilon for numeric
gradients of smooth functionals. Default 1e-6. - \`nmar.grad_d\`:
relative step adjustment for numeric gradients. Default 1e-3.
