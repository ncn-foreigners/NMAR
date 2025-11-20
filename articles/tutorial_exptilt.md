# Exponential Tilting

``` r
library(NMAR)
```

``` r
set.seed(1108)
head(riddles_case1)
#>            x          y     y_true delta
#> 1 0.45763402 -1.7736289 -1.7736289     1
#> 2 0.51906283         NA -0.4653327     0
#> 3 0.46298539         NA -0.5328213     0
#> 4 1.59425312         NA -0.7447774     0
#> 5 1.33728505  0.6900794  0.6900794     1
#> 6 0.07122291  0.7444041  0.7444041     1
```

``` r
y_original <- riddles_case1$y_true # to compare with true values
```

``` r
# Exptilt engine configuration
exptilt_config <- exptilt_engine(
  y_dens = 'normal',
  control = list(maxit = 10),
  stopping_threshold = 0.01,
  standardize = FALSE,
  family = 'logit',
  bootstrap_reps = 50,
  variance_method = 'bootstrap'
)

formula = y ~x
res <- nmar(formula = formula, data = riddles_case1, engine = exptilt_config, trace_level = 0)
```

``` r
print(res)
#> NMAR Result
#> ------------
#> y mean: -1.003197 (0.000002)
#> Converged: TRUE 
#> Variance method: bootstrap 
#> Estimator: exponential_tilting
```

``` r
coef(res)
#> (Intercept)           y 
#>   0.8629256  -0.1722635
```

``` r
cat('True Y mean:          ', sprintf('%.4f', mean(y_original)), '\n')
#> True Y mean:           -1.0105
est <- as.numeric(res$y_hat)
se <- res$se
cat('Est Y mean (NMAR):    ', sprintf('%.4f', est),
    '  3σ interval: (', sprintf('%.4f', est - 1.5 * se),
    ', ', sprintf('%.4f', est + 1.5 * se), 'σ=', sprintf('%.4f', se), ')\n')
#> Est Y mean (NMAR):     -1.0032   3σ interval: ( -1.0032 ,  -1.0032 σ= 0.0000 )
cat('Naive Y mean (MAR):   ', sprintf('%.4f', mean(riddles_case1$y, na.rm = T)), '\n')
#> Naive Y mean (MAR):    -1.0716
```

``` r
# if survey is installed
if (requireNamespace("survey", quietly = TRUE)) {

  library('survey')
  surv_test_weights <- runif(length(riddles_case1$y), 0.9, 1.1)
  surv_test_data <- svydesign(ids = ~1, data = riddles_case1, weights = ~surv_test_weights)

  exptilt_config <- exptilt_engine(
     standardize = FALSE,
     on_failure = "error",
     bootstrap_reps = 50,
     supress_warnings = FALSE,
     auxiliary_means = NULL,
     control = list(maxit = 15, method = "Newton"),
     family = "logit",
     y_dens = "normal",
     variance_method = "delta",
     stopping_threshold = 1e-5
)


  formula = y ~ x
  res <- nmar(formula = formula, data = riddles_case1, engine = exptilt_config)
  print(res)

}
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survival'
#> The following object is masked from 'package:future':
#> 
#>     cluster
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
#> Warning: Delta variance failed to evaluate; using bootstrap instead.
#> NMAR Result
#> ------------
#> y mean: -1.003204 (0.000000)
#> Converged: TRUE 
#> Variance method: delta 
#> Estimator: exponential_tilting
```
