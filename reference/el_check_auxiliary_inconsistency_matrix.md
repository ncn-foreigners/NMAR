# Check auxiliary means consistency against respondents sample support.

Computes a simple z-score diagnostic comparing user-supplied auxiliary
means to the respondents' sample means.

## Usage

``` r
el_check_auxiliary_inconsistency_matrix(
  auxiliary_matrix_resp,
  provided_means = NULL
)
```

## Arguments

- auxiliary_matrix_resp:

  Respondent-side auxiliary design matrix.

- provided_means:

  Optional named numeric vector of auxiliary means aligned to the matrix
  columns.

## Value

list(max_z = numeric(1) or NA, cols = character())
