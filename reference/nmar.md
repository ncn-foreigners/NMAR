# Not Missing at Random (NMAR) Estimation

Provides a unified interface for Not Missing at Random (NMAR)
estimation. This function orchestrates the estimation process by
validating inputs and dispatching to the appropriate engine based on the
provided \`engine\` object. It ensures that all necessary data and model
specifications are correctly formatted before computation begins.

## Usage

``` r
nmar(formula, data, engine, trace_level = 0)
```

## Arguments

- formula:

  A two-sided formula of the form \`y_miss ~ aux1 + aux2 \| z1 + z2\`.
  The left-hand side is the outcome (with \`NA\` values indicating
  nonresponse). The right-hand side is split by \`\|\` into two parts: -
  left of \`\|\`: auxiliary variables (enter moment constraints); -
  right of \`\|\`: missingness (response) model predictors (enter the
  missingness model only). If \`\|\` is omitted, only auxiliary
  variables are used for both parsing and printing. The outcome variable
  is implicitly included in the response model.

- data:

  A data frame or \`survey.design\` containing the variables referenced
  by the formula.

- engine:

  An engine configuration object, typically created by an engine
  constructor function like \`exptilt()\`. This object defines the
  specific NMAR estimation method and its parameters. It must inherit
  from class \`nmar_engine\`.

- trace_level:

  Integer 0-3; controls verbosity level during estimation (default: 1):

  - 0: No output (silent mode)

  - 1: Major steps only (initialization, convergence, final results)

  - 2: Moderate detail (iteration summaries, key diagnostics)

  - 3: Full detail (all diagnostics, intermediate values)

## Value

An object containing the estimation results, whose structure will be
specific to the \`engine\` used. This might include estimated
parameters, convergence information, and other relevant output from the
chosen NMAR method.
