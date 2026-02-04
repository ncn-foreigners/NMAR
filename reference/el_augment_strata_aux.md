# Strata augmentation for survey designs

Augments the auxiliary design with strata dummies and appends
stratum-share means when user-supplied `auxiliary_means` are present.

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
