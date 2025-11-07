# Extract Weights from NMAR Result

Extract Weights from NMAR Result

## Usage

``` r
# S3 method for class 'nmar_result'
weights(object, scale = c("probability", "population"), ...)
```

## Arguments

- object:

  An object of class `nmar_result`

- scale:

  Character: `"probability"` (default) or `"population"`.

  `"probability"`

  :   Returns p_i = w_tilde_i / sum_j w_tilde_j where sum_i p_i = 1
      (exact). This is the paper's canonical form (QLS 2002, Eq. 11).
      Use for computing means: y_bar = sum_i p_i \* y_i

  `"population"`

  :   Returns w_i = N_pop \* p_i where sum_i w_i = N_pop (exact). This
      follows survey package conventions. Use for computing totals:
      T_hat = sum_i w_i \* y_i = N_pop \* y_bar

- ...:

  Additional arguments (ignored)

## Value

Numeric vector of weights with length equal to number of respondents

## Details

The empirical likelihood estimator computes unnormalized masses
w_tilde_i = d_i / D_i that satisfy the constraint sum_i w_tilde_i =
sum_i d_i (without trimming). This method provides two standardized
representations:

**Mathematical guarantees** (hold even with trimming):

- `sum(weights(object, scale = "probability")) = 1` (within machine
  precision)

- `sum(weights(object, scale = "population")) = N_pop` (within machine
  precision)

- `weights(object, "population") = N_pop * weights(object, "probability")`
  (exact)

**Trimming effects**: When `trim_cap < Inf` and trimming is active, the
normalization identity sum_i w_tilde_i = sum_i d_i is violated. However,
this method still returns weights with correct sums by using the
formula: \$\$w_i = N_pop \* w_tilde_i / sum_j w_tilde_j\$\$

## References

Qin, J., Leung, D., & Shao, J. (2002). Estimation with survey data under
nonignorable nonresponse or informative sampling. *Journal of the
American Statistical Association*, 97(457), 193-200.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- nmar(y_miss ~ x, data = df, engine = el_engine())

# Probability weights (default): sum to 1
w_prob <- weights(res)
sum(w_prob) # Exactly 1.0

# Population weights: sum to N_pop
w_pop <- weights(res, scale = "population")
sum(w_pop) # Exactly nrow(df)

# Relationship (exact):
all.equal(w_pop, nrow(df) * w_prob)
} # }
```
