# Validate Data for NMAR Analysis

A robust function to validate a data frame or survey object before
performing NMAR analysis. The function checks for common errors like
missing variables, data type inconsistencies, and inappropriate variable
overlaps. It provides detailed, actionable error messages to facilitate
debugging.

## Usage

``` r
validate_data(
  data,
  outcome_variable,
  covariates_for_outcome,
  covariates_for_missingness = character(),
  allow_outcome_in_missingness = FALSE,
  allow_covariate_overlap = FALSE,
  allow_respondents_only = FALSE
)
```

## Arguments

- data:

  A data frame or a survey object.

- outcome_variable:

  A string specifying the outcome variable, which is expected to contain
  NA values.

- covariates_for_outcome:

  A character vector of covariates explaining the outcome.

- covariates_for_missingness:

  A character vector of covariates explaining missingness.

- allow_outcome_in_missingness:

  Logical; allow the outcome to also appear in the response-model
  covariates (default \`FALSE\`).

- allow_covariate_overlap:

  Logical; allow overlap between outcome and response covariate sets
  (default \`FALSE\`).

- allow_respondents_only:

  Logical; allow datasets with no missing outcome values
  (respondents-only). When TRUE, the caller is expected to provide any
  additional information required by the engine (e.g., total sample
  size).

## Value

Returns \`invisible(NULL)\` on success, stopping with a descriptive
error on failure.
