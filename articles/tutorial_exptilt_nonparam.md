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
# em_formula <- list(
#   outcome = ~ NULL, #TODO - This part is useless
#   covariates_outcome = ~ Gender,
#   covariates_missingness = ~ Age_group
# )
em_formula <- Voted_A + Voted_B + Other ~ Gender
```

``` r
# Exptilt Nonparam engine configuration
np_em_config <- exptilt_nonparam_engine(
# outcome_cols = c("Voted_A", "Voted_B", "Other"),
  refusal_col = "Refusal",
  max_iter = 100,
  tol_value = 1e-6
)
```

``` r
results_em_np <- nmar(formula = Voted_A + Voted_B + Other ~ Gender | Age_group, data = voting_data_example, engine = np_em_config)
```

``` r
print(results_em_np$processed_data)
#>   Gender Age_group Voted_A Voted_B Other Refusal Total m_est_Voted_A
#> 1   Male     20-29      93     115     4      28   240      15.44272
#> 2   Male     30-39     104     233     8      82   427      33.33023
#> 3   Male     40-49     146     295     5      49   495      26.33477
#> 4   Male      >=50     560     350     3     174  1087     145.21679
#> 5 Female     20-29     106     159     8      62   335      43.31161
#> 6 Female     30-39     129     242     5      70   446      51.92976
#> 7 Female     40-49     170     262     5      69   506      54.01178
#> 8 Female       50+     501     218     7     211   937     192.18917
#>   m_est_Voted_B m_est_Other
#> 1      3.929895    8.627389
#> 2     15.367505   33.302260
#> 3     10.950682   11.714545
#> 4     18.678383   10.104824
#> 5      7.426093   11.262294
#> 6     11.135424    6.934817
#> 7      9.514932    5.473284
#> 8      9.559000    9.251829
```
