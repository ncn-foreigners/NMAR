# Riddles Simulation, Case 2: Exponential Mean

A simulated dataset of 500 observations based on Simulation Study I
(Model 1, Case 2) of Riddles, Kim, and Im (2016). The data features a
nonignorable nonresponse (NMAR) mechanism where the response probability
depends on the study variable \`y\`.

## Usage

``` r
riddles_case2
```

## Format

A data frame with 500 rows and 4 variables:

- x:

  Numeric. The auxiliary variable, x ~ Normal(0, 0.5).

- y:

  Numeric. The study variable with nonignorable nonresponse. \`y\`
  contains \`NA\`s for nonrespondents.

- y_true:

  Numeric. The complete, true value of \`y\` before missingness was
  introduced.

- delta:

  Integer. The response indicator (1 = responded, 0 = nonresponse).

## Source

Riddles, M. K., Kim, J. K., & Im, J. (2016). A
Propensity-Score-Adjustment Method for Nonignorable Nonresponse. Journal
of Survey Statistics and Methodology, 4(1), 1-31.

## Details

This dataset was generated using the following model parameters (n =
500):

- Density for x::

  x ~ Normal(mean = 0, variance = 0.5)

- Density for error::

  e ~ Normal(mean = 0, variance = 0.9)

- True Model (Case 2)::

  y_true = -2 + 0.5 \* exp(x) + e

- Response Model (NMAR)::

  logit(pi) = 0.8 - 0.2 \* y_true
