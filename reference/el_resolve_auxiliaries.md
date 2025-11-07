# Resolve auxiliary design and population means

Centralizes the construction of the auxiliary model matrix on
respondents and the corresponding population means used for centering
(X - mu_x).

## Usage

``` r
el_resolve_auxiliaries(
  full_data,
  respondent_data,
  aux_formula,
  auxiliary_means
)
```

## Arguments

- full_data:

  data.frame or survey.design with all sampled units.

- respondent_data:

  data.frame of respondents only.

- aux_formula:

  RHS-only formula for auxiliaries (no intercept) or NULL.

- auxiliary_means:

  optional named numeric vector of population means in the auxiliary
  model-matrix columns.

## Value

list(matrix, means, has_aux) where \`matrix\` is the respondent-side
auxiliary design on the unscaled space, \`means\` is a named numeric
vector aligned to its columns (or NULL), and \`has_aux\` is a logical
flag.

## Details

Rules: - If `aux_formula` is NULL or results in zero columns,
auxiliaries are disabled (returns an empty matrix and NULL means). - If
`auxiliary_means` is provided, it is matched/reordered to the respondent
design columns; unmatched names are dropped. If nothing matches,
auxiliaries are disabled. - Else, compute means from the full data: \*
survey.design: design-weighted column means using `weights(full_data)`
\* data.frame: simple column means Factor levels absent among
respondents are dropped via intersection of model-matrix column names
between full and respondent data.
