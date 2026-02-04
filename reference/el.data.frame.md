# Empirical likelihood for data frames

Empirical likelihood for data frames

## Usage

``` r
# S3 method for class 'data.frame'
el(
  data,
  formula,
  auxiliary_means = NULL,
  standardize = TRUE,
  trim_cap = Inf,
  control = list(),
  on_failure = c("return", "error"),
  variance_method = c("bootstrap", "none"),
  bootstrap_reps = 500,
  n_total = NULL,
  start = NULL,
  trace_level = 0,
  family = logit_family(),
  ...
)
```
