# Exponential Tilting (Nonparametric)

``` r
library(NMAR)
```

``` r
# Test data (Riddles 2016)
voting_data_example <- data.frame(
  Gender = rep(c("Male", "Male", "Male", "Male", "Female", "Female", "Female", "Female"), 1),
  Age_group = c("20-29", "30-39", "40-49", ">=50", "20-29", "30-39", "40-49", "50+"),
  Voted_A = c(93, 104, 146, 560, 106, 129, 170, 501),
  Voted_B = c(115, 233, 295, 350, 159, 242, 262, 218),
  Other = c(4, 8, 5, 3, 8, 5, 5, 7),
  Refusal = c(28, 82, 49, 174, 62, 70, 69, 211),
  Total = c(240, 427, 495, 1087, 335, 446, 506, 937)
)
```

``` r
np_em_config <- exptilt_nonparam_engine(
  refusal_col = "Refusal",
  max_iter = 100,
  tol_value = 0.001
)

em_formula <- Voted_A + Voted_B + Other ~ Gender | Age_group

results_em_np <- nmar(formula = em_formula, data = voting_data_example, engine = np_em_config, trace_level = 0)
```

``` r
print(results_em_np$data_final)
#>   Gender Age_group  Voted_A  Voted_B    Other Total
#> 1   Male     20-29 108.4560 118.8595 12.68447   240
#> 2   Male     30-39 137.3696 248.0971 41.53330   427
#> 3   Male     40-49 172.4092 305.7756 16.81518   495
#> 4   Male      >=50 705.4611 368.3589 13.18002  1087
#> 5 Female     20-29 149.2433 166.3526 19.40414   335
#> 6 Female     30-39 180.9255 253.0418 12.03269   446
#> 7 Female     40-49 224.0130 271.4359 10.55110   506
#> 8 Female       50+ 693.1421 227.4771 16.38086   937
```
