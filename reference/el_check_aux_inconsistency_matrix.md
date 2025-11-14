# Check auxiliary means consistency against respondents' sample support

Check auxiliary means consistency against respondents' sample support

## Usage

``` r
el_check_aux_inconsistency_matrix(
  aux_matrix_resp,
  provided_means = NULL,
  threshold = 8
)
```

## Arguments

- aux_matrix_resp:

  Respondent-side auxiliary design matrix.

- provided_means:

  Optional named numeric vector of auxiliary means aligned to the matrix
  columns.

- threshold:

  numeric, z-score threshold for flagging

## Value

list(max_z = numeric(1) or NA, cols = character())
