# Launch EL estimation once design matrices are parsed

Runs the shared respondent prep, auxiliary resolution, and solver call
so IID and survey methods stay aligned.

## Usage

``` r
el_run_core_analysis(
  call,
  formula,
  input_spec,
  variance_method,
  auxiliary_means,
  standardize,
  trim_cap,
  control,
  on_failure,
  family,
  bootstrap_reps,
  start,
  trace_level,
  extra_user_args = list()
)
```
