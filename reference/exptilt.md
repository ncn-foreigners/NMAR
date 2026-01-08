# Exponential tilting estimator

Generic for the exponential tilting (ET) estimator under NMAR. Methods
are provided for \`data.frame\` and \`survey.design\`.

## Usage

``` r
exptilt(data, ...)
```

## Arguments

- data:

  A \`data.frame\` or a \`survey.design\`.

- ...:

  Passed to class-specific methods.

## Value

An engine-specific NMAR result object (for example
`nmar_result_exptilt`).

## See also

\`exptilt.data.frame()\`, \`exptilt.survey.design()\`,
\`exptilt_engine()\`
