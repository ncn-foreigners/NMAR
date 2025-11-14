# Exponential Tilting

``` r
library(NMAR)
```

``` r
# Function to generate testdata, basing on riddle(2016)
generate_test_data <- function(n_rows = 500, n_cols = 1, case = 1, x_var = 0.5, eps_var = 0.9, a = 0.8, b = -0.2) {
# Generate X variables - fixed to match comparison
  X <- as.data.frame(replicate(n_cols, rnorm(n_rows, 0, sqrt(x_var))))
  colnames(X) <- paste0("x", 1:n_cols)

# Generate Y - fixed coefficients to match comparison
  eps <- rnorm(n_rows, 0, sqrt(eps_var))
  if (case == 1) {
# Use fixed coefficient of 1 for all x variables to match: y = -1 + x1 + epsilon
    X$Y <- as.vector(-1 + as.matrix(X) %*% rep(1, n_cols) + eps)
  }
  else if (case == 2) {
    X$Y <- -2 + 0.5 * exp(as.matrix(X) %*% rep(1, n_cols)) + eps
  }
  else if (case == 3) {
    X$Y <- -1 + sin(2 * as.matrix(X) %*% rep(1, n_cols)) + eps
  }
  else if (case == 4) {
    X$Y <- -1 + 0.4 * as.matrix(X)^3 %*% rep(1, n_cols) + eps
  }

  Y_original <- X$Y

# Missingness mechanism - identical to comparison
  pi_obs <- 1 / (1 + exp(-(a + b * X$Y)))

# Create missing values
  mask <- runif(nrow(X)) > pi_obs
  mask[1] <- FALSE # Ensure at least one observation is not missing
  X$Y[mask] <- NA

  return(list(X = X, Y_original = Y_original))
}
```

``` r
set.seed(1108)
res_test_data <- generate_test_data(n_rows = 500, n_cols = 1, case = 1)
x <- res_test_data$X
y_original <- res_test_data$Y_original # to compare with true values
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

formula = Y ~x1 | Y
res <- nmar(formula = formula, data = x, engine = exptilt_config, trace_level = 3)
#> Warning: Outcome variable (Y) found in missingness predictors; Performance with / without Y on the right side is the same
#> [STEP-1] ============================================================ 
#> [STEP-1]   EXPTILT ESTIMATION STARTED 
#> [STEP-1] ============================================================ 
#> [INFO] Running with trace_level = 3 
#> [INFO] Formula: Y ~ x1 | Y 
#> [INFO]  
#> [INFO] -- DATA SUMMARY -- 
#> [INFO]   Total observations:      500 
#> [INFO]   Respondents:             368 (73.6%) 
#> [INFO]   Non-respondents:         132 (26.4%) 
#> [INFO]  
#> [INFO] -- MODEL SPECIFICATION -- 
#> [INFO]   Outcome variable:        Y 
#> [INFO]   Auxiliary variables:     x1 
#> [INFO]   Missingness predictors:  (intercept only) 
#> [INFO]   Response model:          logit 
#> [INFO]   Outcome density:         normal 
#> [INFO]   Standardization:         disabled 
#> [INFO]  
#> [INFO] -- PARAMETER INITIALIZATION -- 
#> [INFO]   Number of parameters:    2 
#> [DETAIL-3]  
#> [INFO]   Initial values: 
#> [INFO]     (Intercept)               = -0.0432  [response model intercept] 
#> [INFO]     Y                         =  0.0584  [outcome effect on response] 
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
#> x1           1.11333    0.06712   16.59   <2e-16 ***
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
#> [INFO]   Early stopping threshold: 0.0100 
#> [INFO]   Control parameters: 
#> [INFO]     maxit = 10 
#> [INFO]   Solving... 
#> [INFO]  
#> [INFO]   OK Converged 
#> [INFO]   Iterations:               5 
#> [INFO]   Termination code:         1 
#> [INFO]   Max |score|:              0.000000 
#> [INFO]   Final parameter estimates (scaled): 
#> [INFO]     (Intercept)          =  0.862931  [response intercept] 
#> [INFO]     Y                    = -0.172256  [outcome -> response] 
#> [INFO]  
#> [INFO] -- VARIANCE ESTIMATION (Bootstrap) -- 
#> [INFO]   Bootstrap replications:   50 
#> [INFO]   OK Bootstrap complete 
#> [INFO]   Standard error:           0.062760 
#> [INFO]  
#> [RESULT] ============================================================ 
#> [RESULT]   ESTIMATION COMPLETE 
#> [RESULT] ============================================================ 
#> [RESULT]   Mean estimate:            -1.003200 
#> [RESULT]   Standard error:           0.062760 
#> [RESULT]   95% CI:                   [-1.126209, -0.880191] 
#> [INFO]  
#> [INFO]   Response model coefficients: 
#> [INFO]     (Intercept)         :  0.862931  (Intercept) 
#> [INFO]     Y                   : -0.172256  (Effect of Y on response prob.) 
#> [RESULT] ============================================================
```

``` r
print(res)
#> NMAR Result
#> ------------
#> Y mean: -1.003200 (0.062760)
#> Converged: TRUE 
#> Variance method: bootstrap 
#> Estimator: exponential_tilting 
#> Sample size: 500 (respondents: 368)
```

``` r
coef(res)
#> (Intercept)           Y 
#>   0.8629305  -0.1722559
```

``` r
cat('True Y mean:          ', sprintf('%.4f', mean(y_original)), '\n')
#> True Y mean:           -1.0105
est <- as.numeric(res$y_hat)
se <- res$se
cat('Est Y mean (NMAR):    ', sprintf('%.4f', est),
    '  3σ interval: (', sprintf('%.4f', est - 1.5 * se),
    ', ', sprintf('%.4f', est + 1.5 * se), 'σ=', sprintf('%.4f', se), ')\n')
#> Est Y mean (NMAR):     -1.0032   3σ interval: ( -1.0973 ,  -0.9091 σ= 0.0628 )
cat('Naive Y mean (MAR):   ', sprintf('%.4f', mean(x[!is.na(x$Y), 'Y'])), '\n')
#> Naive Y mean (MAR):    -1.0716
```

``` r
# if survey is installed
if (requireNamespace("survey", quietly = TRUE)) {

  library('survey')
  surv_test_weights <- runif(nrow(x), 0.9, 1.1)
  surv_test_data <- svydesign(ids = ~1, data = x, weights = ~surv_test_weights)

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


  formula = Y ~ x1
  res <- nmar(formula = formula, data = x, engine = exptilt_config)
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
#> Y mean: -1.003204 (0.063102)
#> Converged: TRUE 
#> Variance method: delta 
#> Estimator: exponential_tilting 
#> Sample size: 500 (respondents: 368)
```
