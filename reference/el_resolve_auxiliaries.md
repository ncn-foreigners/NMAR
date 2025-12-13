# EL auxiliary design resolution and population means

Computes the respondent-side auxiliary matrix and the population means
vector used for centering \\X - \mu_x\\. When `auxiliary_means` is
supplied, only respondent rows are required to be fully observed; NA
values are permitted on nonrespondent rows. When `auxiliary_means` is
`NULL`, auxiliaries must be fully observed in the full data used to
estimate population means.

## Usage

``` r
el_resolve_auxiliaries(
  aux_design_full,
  respondent_mask,
  auxiliary_means,
  weights_full = NULL
)
```
