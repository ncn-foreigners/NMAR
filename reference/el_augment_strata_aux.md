# Strata augmentation for survey designs

Augments the auxiliary design with strata dummies (dropping one level)
and appends stratum-share means when weights and N_pop are available.
Intended for survey workflows only.

## Usage

``` r
el_augment_strata_aux(
  aux_design_full,
  strata_factor,
  weights_full,
  N_pop,
  auxiliary_means
)
```
