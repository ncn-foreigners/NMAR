# Build the combined EL input specification

Constructs the parsed design matrices and augments the data with the
respondent indicator so iid and survey entry points can share the same
downstream workflow.

## Usage

``` r
el_build_input_spec(
  formula,
  data,
  weights_full = NULL,
  population_total = NULL,
  population_total_supplied = FALSE,
  is_survey = FALSE,
  design_object = NULL,
  auxiliary_means = NULL
)
```
