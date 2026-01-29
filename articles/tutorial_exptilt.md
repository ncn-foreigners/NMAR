# Exponential Tilting

``` r
library(NMAR)
```

## Overview

This vignette demonstrates the Parametric Exponential Tilting estimator
for handling Not Missing at Random (NMAR). This method implements the
maximum likelihood approach described in Riddles, Kim, and Im (2016)

## Input data

This method is applicable when: - Data is individual-level with both the
outcome variable and auxiliary covariates. - The outcome variable
follows a known distribution (e.g., normal).

### Required Structure

- **Outcome Variable**: A continuous variable subject to nonresponse
  (e.g., `y`).
- **Covariates for outcome model**: Variables trying to model Outcome
  Variable (e.g., `x`).
- **Instrumental Variable**: (optional) A variable affecting response
  probability but not directly related to the outcome.

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

## Engine Configuration

The estimation is configured using `exptilt_engine`. Critical parameters
include the assumed outcome density for respondents and the variance
estimation method

The model is specified via a formula: Outcome ~ Response_Covariates \|
Instrument.

- **LHS**: The outcome variable (Y) containing missing values.
- **RHS** (Left of \|): Covariates included in the response propensity
  model .
- **RHS** (Right of \|): The instrumental variable. The method assumes
  this variable predicts the outcome but is conditionally independent of
  the response status .

``` r
# Exptilt engine configuration
exptilt_config <- exptilt_engine(
  y_dens = 'normal', # Assume f(y|x, delta=1) is Normal
  family = 'logit', # Logit link for response probability
  variance_method = 'bootstrap',
  bootstrap_reps = 5,
  control = list(maxit = 10), # Solver control parameters
  stopping_threshold = 0.01, # Convergence threshold
  standardize = FALSE
)

formula = y ~x
res <- nmar(formula = formula, data = riddles_case1, engine = exptilt_config, trace_level = 0)
```

## Results

The function returns an object containing the estimated coefficients for
the response model and the adjusted mean estimate for the outcome.

The adjusted mean is calculated using propensity score weighting, where
the weights are derived from the fitted nonignorable response model

``` r
print(res)
#> NMAR Result (Exponential tilting)
#> -------------------------------
#> y mean: -1.003197 (0.077755)
#> Converged: TRUE 
#> Variance method: bootstrap 
#> Estimator: exponential_tilting 
#> 
#> Exptilt diagnostics:
#>   Loss value: 0.000000 0.000000 
#>   Iterations: 5 
#>   Variance method: bootstrap 
#>   Bootstrap reps: 5 
#>   Stopping threshold: 0.010000
```

``` r
coef(res)
#> (Intercept)           y 
#>   0.8629256  -0.1722635
```

Comparing the results to the true population mean (known from the
synthetic generation) demonstrates the bias correction.

``` r
y_true_mean <- mean(riddles_case1$y_true)
est_mean <- as.numeric(res$y_hat)
se <- res$se

cat('Estimated (NMAR) Mean:', sprintf('%.4f', est_mean),
    ' (SE:', sprintf('%.4f', se), ')\n')
#> Estimated (NMAR) Mean: -1.0032  (SE: 0.0778 )
cat('True Y Mean:          ', sprintf('%.4f', y_true_mean), '\n')
#> True Y Mean:           -1.0105
cat('Naive (MCAR) Mean:     ', sprintf('%.4f', mean(riddles_case1$y, na.rm = TRUE)), '\n')
#> Naive (MCAR) Mean:      -1.0716
```

## Survey designs

The engine supports survey.design objects. When provided, the algorithm
incorporates sampling weights into the likelihood estimation and the
final propensity score weighting.

The weights affect both the estimation of the respondent density f(y∣x)
and the solution to the mean score equation, ensuring population-level
inference

``` r
# if survey is installed
if (requireNamespace("survey", quietly = TRUE)) {

  library('survey')
  surv_test_weights <- runif(length(riddles_case1$y), 0.9, 1.1)
  surv_test_data <- svydesign(ids = ~1, data = riddles_case1, weights = ~surv_test_weights)

  exptilt_config <- exptilt_engine(
     standardize = FALSE,
     on_failure = "error",
     bootstrap_reps = 5,
     supress_warnings = FALSE,
     control = list(maxit = 15, method = "Newton"),
     family = "logit",
     y_dens = "normal",
     variance_method = "bootstrap",
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
#> [INFO]     (Intercept)               = -0.0496  [response model intercept] 
#> [INFO]     y                         =  0.0673  [outcome effect on response] 
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
#> [INFO]     (Intercept)          =  0.853226  [response intercept] 
#> [INFO]     y                    = -0.179754  [outcome -> response] 
#> [INFO]  
#> [INFO] -- VARIANCE ESTIMATION (Bootstrap) -- 
#> [INFO]   Bootstrap replications:   5 
#> [INFO]   OK Bootstrap complete 
#> [INFO]   Standard error:           0.072753 
#> [INFO]  
#> [RESULT] ============================================================ 
#> [RESULT]   ESTIMATION COMPLETE 
#> [RESULT] ============================================================ 
#> [RESULT]   Mean estimate:            -0.998354 
#> [RESULT]   Standard error:           0.072753 
#> [RESULT]   95% CI:                   [-1.140950, -0.855758] 
#> [INFO]  
#> [INFO]   Response model coefficients: 
#> [INFO]     (Intercept)         :  0.853226  (Intercept) 
#> [INFO]     y                   : -0.179754  (Effect of y on response prob.) 
#> [RESULT] ============================================================ 
#> NMAR Result (Exponential tilting)
#> -------------------------------
#> y mean: -0.998354 (0.072753)
#> Converged: TRUE 
#> Variance method: bootstrap 
#> Estimator: exponential_tilting 
#> 
#> Exptilt diagnostics:
#>   Loss value: 0.000000 0.000000 
#>   Iterations: 3 
#>   Variance method: bootstrap 
#>   Bootstrap reps: 5 
#>   Stopping threshold: 1.000000
```
