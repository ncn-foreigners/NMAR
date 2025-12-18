# Exponential Tilting (Nonparametric)

``` r
library(NMAR)
```

## Overview

This vignette provides a practical guide to the Nonparametric
Exponential Tilting estimator for handling nonignorable nonresponse
(NMAR) in survey data. This method implements the “Fully Nonparametric
Approach” detailed in Appendix 2 of Riddles, Kim, and Im (2016)

The Nonparametric Exponential Tilting engine estimates Not Missing at
Random (NMAR) bias using aggregated categorical data (contingency
tables). Unlike other engines in this package, **it does not process
individual rows** (microdata) but summary counts per stratum.

## Input Data

This estimator is specifically designed for scenarios where: - **Data is
Aggregated**: The input consists of respondent and nonrespondent counts
per stratum, rather than microdata. - **Variables are Categorical**:
Both the outcome of interest and the auxiliary covariates are discrete

### Required Structure:

- **Stratification Columns**: Categorical variables defining the groups
  (e.g., `Gender`, `Age_group`)
- **Outcome Counts**: Columns representing counts of different outcomes
  (e.g., `Voted_A`, `Voted_B`, `Other`)
- **Refusal Count**: A column indicating the number of nonresponses
  (e.g., `Refusal`)

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

## Engine Configuration

The `exptilt_nonparam_engine` specifies the column containing
nonrespondent counts and sets convergence criteria for the EM algorithm.

The model is specified via a two-part formula: `Outcome_Counts` ~
`Response_Covariates` \| `Instrument`: - **LHS**: Sum of outcome
columns. - **RHS**: Covariates influencing response probability (left of
\|) and the instrumental variable (right of \|)

``` r
np_em_config <- exptilt_nonparam_engine(
  refusal_col = "Refusal",
  max_iter = 100, # Maximum EM iterations
  tol_value = 0.1 # Convergence tolerance
)

em_formula <- Voted_A + Voted_B + Other ~ Gender | Age_group

results_em_np <- nmar(formula = em_formula, data = voting, engine = np_em_config, trace_level = 0)
```

## Results

The data_final object returns the reconstructed population. The
algorithm redistributes those counts into the outcome columns based on
the estimated nonresponse odds. These adjusted counts allow for the
direct calculation of corrected population proportions

``` r
print(results_em_np$data_final)
#>    Voted_A  Voted_B     Other Gender Age_group
#> 1 105.7271 129.7191  4.553824   Male     20-29
#> 2 129.8415 287.1473 10.011131   Male     30-39
#> 3 162.7556 326.6639  5.580552   Male     40-49
#> 4 669.4366 413.9703  3.593145   Male       50+
#> 5 134.0864 191.1140  9.799682 Female     20-29
#> 6 157.4218 282.6429  5.935294 Female     30-39
#> 7 201.3671 298.8497  5.783269 Female     40-49
#> 8 658.0469 270.0901  8.862970 Female       50+
```

Beyond the adjusted data, the result object contains diagnostic
information and intermediate matrices used in the EM algorithm. These
can be inspected to assess convergence and model internals.

Key components include: - `$fit_stats`: Contains the number of
iterations performed and the convergence status (TRUE/FALSE). -
`$loss_value`: The final sum of absolute differences in the odds matrix,
indicating the precision of convergence relative to tol_value. -
`$p_hat_y_given_x_matrix`: The estimated conditional probabilities of
the outcome given covariates, calculated from the respondent data. -
`$n_y_x_matrix` and `$m_x_vec`: The internally used matrices for
weighted respondent and nonrespondent counts, respectively.

``` r
print(results_em_np$fit_stats)
#> $iterations
#> [1] 2
#> 
#> $converged
#> [1] TRUE
```

## Survey designs

The engine supports survey.design objects. When provided, the algorithm
extracts sampling weights and scales the observed counts prior to the EM
procedure, ensuring estimates reflect the target population.

``` r
# Example: Integration with the survey package
if (requireNamespace("survey", quietly = TRUE)) {
  library(survey)

# Simulate sampling weights
  set.seed(42)

  des <- svydesign(ids = ~1, weights = abs(rnorm(nrow(voting), mean = 1, sd = 0.05)), data = voting)

  fit_survey <- nmar(
    formula = em_formula,
    data = des,
    engine = np_em_config,
    trace_level = 0
  )

  print(fit_survey$data_final)
}
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
#>    Voted_A  Voted_B     Other Gender Age_group
#> 1 113.0123 138.5787  4.860519   Male     20-29
#> 2 126.2632 278.9687  9.711809   Male     30-39
#> 3 165.7664 332.5443  5.676708   Male     40-49
#> 4 690.8249 426.8706  3.700616   Male       50+
#> 5 136.7994 194.9707 10.001422 Female     20-29
#> 6 156.5905 281.1372  5.905697 Female     30-39
#> 7 216.5900 321.4294  6.222056 Female     40-49
#> 8 654.9408 268.7998  8.824640 Female       50+
```
