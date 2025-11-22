# Extract Weights from NMAR Result

Extract Weights from NMAR Result

## Usage

``` r
# S3 method for class 'nmar_result'
weights(object, scale = c("probability", "population"), ...)
```

## Arguments

- object:

  An object of class `nmar_result`.

- scale:

  Character: `"probability"` (default) or `"population"`.

  `"probability"`

  :   Returns \\p_i = \tilde w_i / \sum_j \tilde w_j\\ where \\\sum_i
      p_i = 1\\. This is the canonical form for computing means, for
      example \\\bar y = \sum_i p_i y_i\\.

  `"population"`

  :   Returns \\w_i = N\_\mathrm{pop} p_i\\ where \\\sum_i w_i =
      N\_\mathrm{pop}\\. This follows survey conventions and can be used
      for totals \\\hat T = \sum_i w_i y_i = N\_\mathrm{pop} \bar y\\.

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of weights with length equal to the number of
respondents.

## Details

By convention, NMAR engines that expose analysis weights store
unnormalized respondent masses \\\tilde w_i\\ on the analysis scale in
the `weights_info$values` component, and record the population size
\\N\_\mathrm{pop}\\ in `sample$n_total`. For the empirical likelihood
(EL) engine these masses are \\\tilde w_i = d_i / D_i(\theta)\\ as in
Qin, Leung, and Shao (2002); exponential tilting engines may store the
corresponding tilted masses.

When an engine follows this convention, this helper standardizes those
masses into probability and population weights with the following
properties (up to floating-point error, and even under trimming):

- `sum(weights(object, scale = "probability")) = 1`;

- `sum(weights(object, scale = "population")) = N_pop`;

- `weights(object, "population") = N_pop * weights(object, "probability")`.

Engines that use different internal weighting schemes can either map
their weights into this convention (by populating `weights_info$values`
and `sample$n_total` accordingly) or provide engine-specific methods
`weights.nmar_result_<method>()` if a different interpretation is
required.

**Trimming effects**: When `trim_cap < Inf` and trimming is active, the
stored unnormalized masses \\\tilde w_i\\ may no longer satisfy
identities such as \\\sum_i \tilde w_i = \sum_i d_i\\. This helper
always renormalizes the stored masses via \$\$w_i = N\_\mathrm{pop}
\tilde w_i / \sum_j \tilde w_j,\$\$ so that the returned probability-
and population-scale weights satisfy the stated sums regardless of
trimming or denominator floors used internally.

## References

Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
under nonignorable nonresponse or informative sampling. *Journal of the
American Statistical Association*, 97(457), 193-200.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- nmar(y_miss ~ x, data = df, engine = el_engine())

# Probability weights (default): sum to 1
w_prob <- weights(res)
sum(w_prob) # 1 (up to numerical precision)

# Population weights: sum to N_pop
w_pop <- weights(res, scale = "population")
sum(w_pop) # nrow(df) for IID data

# Relationship (exact up to floating-point error):
all.equal(w_pop, nrow(df) * w_prob)
} # }
```
