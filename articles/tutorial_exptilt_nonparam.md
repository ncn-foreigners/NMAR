# Exponential Tilting (Nonparametric)

``` r
library(NMAR)
```

``` r
# Test data (Riddles 2016)
head(voting)
#>   Gender Age_group Voted_A Voted_B Other Refusal Total
#> 1   Male     20-29      93     115     4      28   240
#> 2   Male     30-39     104     233     8      82   427
#> 3   Male     40-49     146     295     5      49   495
#> 4   Male       50+     560     350     3     174  1087
#> 5 Female     20-29     106     159     8      62   335
#> 6 Female     30-39     129     242     5      70   446
```

``` r
np_em_config <- exptilt_nonparam_engine(
  refusal_col = "Refusal",
  max_iter = 100,
  tol_value = 0.001
)

em_formula <- Voted_A + Voted_B + Other ~ Gender | Age_group

results_em_np <- nmar(formula = em_formula, data = voting, engine = np_em_config, trace_level = 0)
```

``` r
print(results_em_np$data_final)
#>   Gender Age_group  Voted_A  Voted_B    Other Total
#> 1   Male     20-29 108.4560 118.8595 12.68447   240
#> 2   Male     30-39 137.3696 248.0971 41.53330   427
#> 3   Male     40-49 172.4092 305.7756 16.81518   495
#> 4   Male       50+ 705.4611 368.3589 13.18002  1087
#> 5 Female     20-29 149.2433 166.3526 19.40414   335
#> 6 Female     30-39 180.9255 253.0418 12.03269   446
#> 7 Female     40-49 224.0130 271.4359 10.55110   506
#> 8 Female       50+ 693.1421 227.4771 16.38086   937
```
