# Rebuild a partitioned formula y ~ aux \| response

Given a base formula whose RHS contains auxiliaries only (the normalized
task\$formula) and a character vector of response-only predictors,
construct a new formula that partitions the RHS as \`aux \| response\`
using language objects (no deparse/paste). If \`response_predictors\` is
empty, returns the base formula unchanged.

## Usage

``` r
nmar_rebuild_partitioned_formula(base_formula, response_predictors, env = NULL)
```

## Arguments

- base_formula:

  A two-sided formula (y ~ aux-only) whose environment is set.

- response_predictors:

  Character vector of response-only predictors.

- env:

  Optional formula environment override; defaults to base formula env.
