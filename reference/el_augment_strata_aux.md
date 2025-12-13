# Strata augmentation for survey designs

Augments the auxiliary design with strata dummies (dropping one level)
and appends stratum-share means when user-supplied `auxiliary_means` are
present. This is the Wu-style strategy of adding stratum indicators to
the auxiliary calibration block in pseudo empirical likelihood for
surveys.

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
