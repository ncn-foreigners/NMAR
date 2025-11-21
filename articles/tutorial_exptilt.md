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
#> y mean: -1.003197 (0.077627)
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
#> Est Y mean (NMAR):     -1.0032   3σ interval: ( -1.1196 ,  -0.8868 σ= 0.0776 )
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
     stopping_threshold = 1
)


  formula = y ~ x
  res <- nmar(formula = formula, data = surv_test_data, engine = exptilt_config, trace_level = 3)
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
#> [STEP-1] ============================================================ 
#> [STEP-1]   EXPTILT ESTIMATION STARTED 
#> [STEP-1] ============================================================ 
#> [INFO] Running with trace_level = 3 
#> [INFO] Formula: y ~ x 
#> [INFO]  
#> [INFO] -- DATA SUMMARY -- 
#> [INFO]   Total observations:      500 
#> [INFO]   Respondents:             368 (26.4%) 
#> [INFO]   Non-respondents:         132 (73.6%) 
#> [INFO]  
#> [INFO] -- MODEL SPECIFICATION -- 
#> [INFO]   Outcome variable:        y 
#> [INFO]   Auxiliary variables:     x 
#> [INFO]   Missingness predictors:  (intercept only) 
#> [INFO]   Response model:          logit 
#> [INFO]   Outcome density:         normal 
#> [INFO]   Standardization:         disabled 
#> [INFO]  
#> [INFO] -- PARAMETER INITIALIZATION -- 
#> [INFO]   Number of parameters:    2 
#> [DETAIL-3]  
#> [INFO]   Initial values: 
#> [INFO]     (Intercept)               =  0.0481  [response model intercept] 
#> [INFO]     y                         =  0.0928  [outcome effect on response] 
#> [INFO]  
#> [INFO] -- CONDITIONAL DENSITY ESTIMATION -- 
#> [INFO]   Selected distribution:   normal 
#> [INFO]   Density parameters:      3 
#> [INFO]   Fitted density model summary: 
#> [INFO]  
#> 
#> Call:
#> glm(formula = formula, family = gaussian(), data = data)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -1.06281    0.04866  -21.84   <2e-16 ***
#> x            1.11333    0.06712   16.59   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.8711681)
#> 
#>     Null deviance: 558.57  on 367  degrees of freedom
#> Residual deviance: 318.85  on 366  degrees of freedom
#> AIC: 997.58
#> 
#> Number of Fisher Scoring iterations: 2
#> 
#> [INFO]   Computing density matrices... 
#> [INFO]  
#> [INFO] -- NONLINEAR SOLVER (nleqslv) -- 
#> [INFO]   Early stopping threshold: 1.0000 
#> [INFO]   Method:                   Newton 
#> [INFO]   Control parameters: 
#> [INFO]     maxit = 15 
#> [INFO]   Solving... 
#> [INFO]  
#> [INFO]   OK Converged 
#> [INFO]   Iterations:               3 
#> [INFO]   Termination code:         1 
#> [INFO]   Max |score|:              0.000000 
#> [INFO]   Final parameter estimates (scaled): 
#> [INFO]     (Intercept)          =  0.865942  [response intercept] 
#> [INFO]     y                    = -0.159226  [outcome -> response]
#> Warning: Delta variance may be unreliable with the current sample; using
#> bootstrap instead.
#> [INFO]  
#> [INFO] -- VARIANCE ESTIMATION (Bootstrap) -- 
#> [INFO]   Bootstrap replications:   50 
#> [INFO]   OK Bootstrap complete 
#> [INFO]   Standard error:           0.046478 
#> [INFO]  
#> [RESULT] ============================================================ 
#> [RESULT]   ESTIMATION COMPLETE 
#> [RESULT] ============================================================ 
#> [RESULT]   Mean estimate:            -1.004793 
#> [RESULT]   Standard error:           0.046478 
#> [RESULT]   95% CI:                   [-1.095891, -0.913696] 
#> [INFO]  
#> [INFO]   Response model coefficients: 
#> [INFO]     (Intercept)         :  0.865942  (Intercept) 
#> [INFO]     y                   : -0.159226  (Effect of y on response prob.) 
#> [RESULT] ============================================================ 
#> NMAR Result
#> ------------
#> y mean: -1.004793 (0.046478)
#> Converged: TRUE 
#> Variance method: delta 
#> Estimator: exponential_tilting
```
