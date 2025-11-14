# Enforce respondents-only requirements for a given design

When all rows are respondents we no longer observe the population size
in the data. This helper ensures we were given \`n_total\` (communicated
here via the boolean flag) and, when auxiliary constraints are present,
verifies that population means were supplied as well.

## Usage

``` r
el_require_population_inputs(
  design,
  population_total_supplied,
  auxiliary_means,
  context_label
)
```
