# Engine trait declarations

Public S3 generic returning a small list of declarative flags that
describe how an engine expects input validation to behave. This function
is also used internally by the input pipeline. Users can call it on an
engine object to introspect behaviour such as whether outcome variables
may appear in the missingness model or whether multiple outcomes are
supported.

## Usage

``` r
engine_traits(engine)
```

## Arguments

- engine:

  An object inheriting from class \`nmar_engine\`.

## Value

A named list of trait flags. The current fields are: -
\`allow_outcome_in_missingness\`: logical. -
\`allow_covariate_overlap\`: logical. - \`requires_single_outcome\`:
logical. Engines may add traits over time; callers should use \`\$\`
with care and rely on presence checks when needed.
