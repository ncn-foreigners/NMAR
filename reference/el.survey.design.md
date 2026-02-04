# Empirical likelihood estimator for survey designs

Empirical likelihood estimator for survey designs

## Usage

``` r
# S3 method for class 'survey.design'
el(
  data,
  formula,
  auxiliary_means = NULL,
  standardize = TRUE,
  strata_augmentation = TRUE,
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
