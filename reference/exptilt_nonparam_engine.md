# Nonparametric exponential tilting engine

Build a configuration for the nonparametric exponential-tilting EM
estimator used in NMAR problems with grouped outcomes and refusal
counts. Pass the resulting engine to \[nmar()\] together with the
appropriate formula and data.

## Usage

``` r
exptilt_nonparam_engine(refusal_col, max_iter = 100, tol_value = 1e-06)
```

## Arguments

- refusal_col:

  Column name in \`data\` containing refusal counts.

- max_iter:

  Maximum number of EM iterations.

- tol_value:

  Convergence tolerance for the EM updates.

## Value

A list of class \`c("nmar_engine_exptilt_nonparam", "nmar_engine")\`.
