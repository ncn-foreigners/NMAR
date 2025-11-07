# Prepare common NMAR design components

Prepare common NMAR design components

## Usage

``` r
prepare_nmar_design(
  task,
  standardize = TRUE,
  auxiliary_means = NULL,
  include_response = TRUE,
  include_auxiliary = TRUE,
  data = task$data,
  design_weights = NULL
)
```

## Arguments

- task:

  An object created by \[new_nmar_task()\].

- standardize:

  Logical flag forwarded to engines.

- auxiliary_means:

  Optional named vector of auxiliary means.

- include_response:

  Logical; include response predictors.

- include_auxiliary:

  Logical; include auxiliary predictors.

- data:

  Optional override of the data-frame backing the task.

- design_weights:

  Optional numeric vector of design weights.

## Value

A list containing the trimmed data, predictor sets, weights, survey
design (if applicable), formula, and standardization settings.
