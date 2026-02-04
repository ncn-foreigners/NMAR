# Nonparametric exponential tilting engine for NMAR

Constructs a configuration for the nonparametric exponential tilting
estimator under nonignorable nonresponse. This engine implements the
"Fully Nonparametric Approach" from \*\*Appendix 2\*\* of Riddles et al.
(2016). The estimator uses an Expectation-Maximization (EM) algorithm to
directly estimate the nonresponse odds \\O(x_1, y)\\ for aggregated,
categorical data.

## Usage

``` r
exptilt_nonparam_engine(refusal_col = "", max_iter = 100, tol_value = 1e-06)
```

## Arguments

- refusal_col:

  character; the column name in `data` that contains the aggregated
  counts of non-respondents (refusals).

- max_iter:

  integer; the maximum number of iterations for the EM algorithm.

- tol_value:

  numeric; the convergence tolerance for the EM algorithm. The loop
  stops when the sum of absolute changes in the odds matrix is less than
  this value.

## Value

An engine object of class
`c("nmar_engine_exptilt_nonparam","nmar_engine")`. This is a
configuration list; it is not a fit. Pass it to
[nmar](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).

## Details

This engine is designed for cases where all variables (outcomes \$Y\$,
response predictors \$X_1\$, and instrumental variables \$X_2\$) are
categorical, and the input `data` is pre-aggregated into strata.

The method assumes an instrumental variable \\X_2\\ is available. The
response probability is assumed to depend on \\X_1\\ and \$Y\$, but
\*not\* on \\X_2\\.

The EM algorithm iteratively solves for the nonresponse odds: \$\$
O^{(t+1)}(x_1^\*, y^\*) = \frac{M\_{y^\*x_1^\*}^{(t)}}{N\_{y^\*x_1^\*}}
\$\$ where \\M\_{y^\*x_1^\*}^{(t)}\\ is the expected count of
non-respondents (calculated in the E-step) and \\N\_{y^\*x_1^\*}\\ is
the observed count of respondents for a given stratum \$(x_1, y)\$.

The final output from the
[`nmar`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md)
call is an object containing `data_to_return`, an aggregated data frame
where the original 'refusal' counts have been redistributed into the
outcome columns (e.g., 'Voted_A', 'Voted_B') as expected non-respondent
counts.

## References

Minsun Kim Riddles, Jae Kwang Kim, Jongho Im A
Propensity-score-adjustment Method for Nonignorable Nonresponse *Journal
of Survey Statistics and Methodology*, Volume 4, Issue 2, June 2016,
Pages 215–245. (See \*\*Appendix 2\*\* for this specific method).

## Examples

``` r
# Test data (Riddles 2016, Table 9)
voting_data_example <- data.frame(
  Gender = rep(c("Male", "Male", "Male", "Male", "Female", "Female", "Female", "Female"), 1),
  Age_group = c("20-29", "30-39", "40-49", ">=50", "20-29", "30-39", "40-49", "50+"),
  Voted_A = c(93, 104, 146, 560, 106, 129, 170, 501),
  Voted_B = c(115, 233, 295, 350, 159, 242, 262, 218),
  Other = c(4, 8, 5, 3, 8, 5, 5, 7),
  Refusal = c(28, 82, 49, 174, 62, 70, 69, 211),
  Total = c(240, 427, 495, 1087, 335, 446, 506, 937)
)

np_em_config <- exptilt_nonparam_engine(
  refusal_col = "Refusal",
  max_iter = 100,
  tol_value = 0.001
)

# Formula: Y1 + Y2 + ... ~ X1_vars | X2_vars
# Here, Y = Voted_A, Voted_B, Other
#      x1 = Gender (response model)
#      x2 = Age_group (instrumental variable)
em_formula <- Voted_A + Voted_B + Other ~ Gender | Age_group

# \donttest{
results_em_np <- nmar(
  formula = em_formula,
  data = voting_data_example,
  engine = np_em_config,
  trace_level = 0
)

# View the final adjusted counts
# (Original counts + expected non-respondent counts)
print(results_em_np$data_final)
#>    Voted_A  Voted_B    Other Gender Age_group
#> 1 108.4560 118.8595 12.68447   Male     20-29
#> 2 137.3696 248.0971 41.53330   Male     30-39
#> 3 172.4092 305.7756 16.81518   Male     40-49
#> 4 705.4611 368.3589 13.18002   Male      >=50
#> 5 149.2433 166.3526 19.40414 Female     20-29
#> 6 180.9255 253.0418 12.03269 Female     30-39
#> 7 224.0130 271.4359 10.55110 Female     40-49
#> 8 693.1421 227.4771 16.38086 Female       50+
# }
```
