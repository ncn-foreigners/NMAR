# Nonparametric Exponential Tilting Theory

## Overview

This vignette explains the matrix-vectorized implementation of the
**fully nonparametric** exponential tilting (EXPTILT) estimator, as
described in **Appendix 2** of Riddles et al. (2016). This method is
designed for **fully categorical data**, where the outcomes ($Y$), the
response-model covariates ($X_{1}$), and the instrumental-variable
covariates ($X_{2}$) are all discrete.

Unlike the parametric method, which estimates the parameters $\phi$ of a
response probability function $\pi(x,y;\phi)$, the nonparametric
approach directly estimates the **nonresponse odds** for each stratum.
This is achieved using an **Expectation-Maximization (EM) algorithm** to
find the maximum likelihood estimates of these odds.

The implementation (`exptilt_nonparam.data.frame`) assumes the input
`data` is an **aggregated data frame**, where each row represents a
unique stratum $x^{*} = \left( x_{1},x_{2} \right)$ and contains the
*counts* of respondents for each outcome $y^{*}$ and the *total count*
of nonrespondents.

## Notation and Main Objects

The implementation maps directly to the notation in Appendix 2. Let
$x_{1}$ be the covariates for the response model and $x_{2}$ be the
nonresponse instrumental variable. A full stratum is
$x^{*} = \left( x_{1},x_{2} \right)$, and a response-model-only stratum
is $x_{1}$.

The algorithm is built from a set of fixed matrices (computed once) and
one matrix that is iteratively updated.

### Fixed Objects (Pre-computed)

1.  **Respondent Counts $N_{y^{*}x^{*}}$** (code: `n_y_x_matrix`)
    - **Source:** `data[, outcome_cols]`
    - **Dimensions:**
      $\left( N_{\text{strata}} \times K_{\text{outcomes}} \right)$,
      where $N_{\text{strata}}$ is the number of
      $\left( x_{1},x_{2} \right)$ rows and $K_{\text{outcomes}}$ is the
      number of $Y$ categories.
    - **Definition:** The observed (weighted) count of respondents for
      stratum $x^{*}$ and outcome $y^{*}$.
      $$N_{y^{*}x^{*}} = \sum\limits_{i \in A}d_{i}\delta_{i}I\left( y_{i} = y^{*},x_{i} = x^{*} \right)$$
      (Note: The implementation assumes $d_{i} = 1$ unless the input
      `data` is pre-weighted).
2.  **Nonrespondent Counts $M_{x^{*}}$** (code: `m_x_vec`)
    - **Source:** `data[, refusal_col]`
    - **Dimensions:** $\left( N_{\text{strata}} \times 1 \right)$
    - **Definition:** The observed (weighted) total count of
      nonrespondents for stratum $x^{*}$. $$M_{x^{*}} = m_{x^{*}}$$
3.  **Respondent Proportions ${\widehat{P}}_{y^{*}|x^{*}}$** (code:
    `p_hat_matrix`)
    - **Source:** `n_y_x_matrix / rowSums(n_y_x_matrix)`
    - **Dimensions:**
      $\left( N_{\text{strata}} \times K_{\text{outcomes}} \right)$
    - **Definition:** The conditional probability (proportion) of
      observing outcome $y^{*}$ given stratum $x^{*}$, *among
      respondents*. This is a fixed, observed quantity used in the
      E-Step.
      $${\widehat{p}}_{y^{*}|x^{*}} = \frac{N_{y^{*}x^{*}}}{\sum\limits_{y}N_{yx^{*}}} = pr\left( y = y^{*}|x_{1} = x_{1}^{*},x_{2} = x_{2}^{*},\delta = 1 \right)$$
4.  **Aggregated Respondent Counts $N_{y^{*}x_{1}^{*}}$** (code:
    `n_y_x1_matrix`)
    - **Source:** `aggregate(n_y_x_matrix ~ data$x1_key, ...)`
    - **Dimensions:**
      $\left( N_{x_{1}} \times K_{\text{outcomes}} \right)$, where
      $N_{x_{1}}$ is the number of unique $x_{1}$ strata.
    - **Definition:** The observed (weighted) count of respondents for
      outcome $y^{*}$ in stratum $x_{1}$, summed over the instrument
      $x_{2}$. This is the denominator of the M-Step.
      $$N_{y^{*}x_{1}^{*}} = \sum\limits_{x_{2}}N_{y^{*}{(x_{1}^{*},x_{2})}} = n_{y^{*}x_{1}^{*}}$$

### Iterative Objects (The EM Algorithm)

1.  **Odds Matrix $O^{(t)}\left( x_{1},y \right)$** (code:
    `odds_matrix`)
    - **Dimensions:**
      $\left( N_{x_{1}} \times K_{\text{outcomes}} \right)$
    - **Definition:** The **parameter** being estimated. It represents
      the odds of nonresponse for a given $x_{1}$ stratum and outcome
      $y$. This is updated in the M-Step.
    - **Initialization (Step 0):** $O^{(0)}\left( x_{1},y \right) = 1$
      for all $\left( x_{1},y \right)$.
2.  **Expected Nonrespondent Counts $M_{y^{*}x^{*}}^{(t)}$** (code:
    `m_y_x_matrix`)
    - **Dimensions:**
      $\left( N_{\text{strata}} \times K_{\text{outcomes}} \right)$
    - **Definition:** The *expected* count of nonrespondents for stratum
      $x^{*}$ and outcome $y^{*}$, given the current odds $O^{(t)}$.
      This is computed in the E-Step.
3.  **Aggregated Expected Nonrespondent Counts
    $M_{y^{*}x_{1}^{*}}^{(t)}$** (code: `m_y_x1_matrix`)
    - **Dimensions:**
      $\left( N_{x_{1}} \times K_{\text{outcomes}} \right)$
    - **Definition:** The *expected* count of nonrespondents for outcome
      $y^{*}$ in stratum $x_{1}$, summed over the instrument $x_{2}$.
      This is the numerator of the M-Step.

## The Expectation-Maximization (EM) Algorithm

The function `exptilt_nonparam.data.frame` is a direct implementation of
the EM algorithm in Appendix 2. The goal is to find the
$\text{argmax}O\left( x_{1},y \right)$ that maximizes the observed data
likelihood, which is solved by iterating two steps.

### Step 1.1 (E-Step): Compute Expected Nonrespondent Counts

The E-Step computes the expected breakdown of nonrespondents into
outcome categories. It answers: “Given our current `odds_matrix`
$O^{(t)}$, what is the expected count $M_{y^{*}x^{*}}^{(t)}$ for each
full stratum $x^{*}$?”

- **Formula:**
  $$M_{y^{*}x^{*}}^{(t)} = M_{x^{*}}\frac{{\widehat{p}}_{y^{*}|x^{*}}O^{(t)}\left( x_{1}^{*},y^{*} \right)}{\sum\limits_{y}{\widehat{p}}_{y|x^{*}}O^{(t)}\left( x_{1}^{*},y \right)}$$
- **Implementation (`E-STEP` in code):** This is vectorized. The
  denominator is computed via
  `rowSums(p_hat_matrix * odds_joined_matrix)`. The full calculation is:
  `m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator`

### Step 1.2 (M-Step): Update Odds Matrix

The M-Step updates the nonresponse odds $O\left( x_{1},y \right)$ using
the expected counts from the E-Step.

1.  **Aggregate Expected Counts:** First, the expected nonrespondent
    counts $M_{y^{*}x^{*}}^{(t)}$ are aggregated over the instrument
    $x_{2}$ to get the total expected nonrespondents
    $M_{y^{*}x_{1}^{*}}^{(t)}$ for each $\left( x_{1},y \right)$ cell.
    - **Formula:**
      $M_{y^{*}x_{1}^{*}}^{(t)} = \sum_{x_{2}}M_{y^{*}{(x_{1}^{*},x_{2})}}^{(t)}$
    - **Implementation:**
      `m_y_x1_matrix <- aggregate(m_y_x_matrix ~ data$x1_key, ...)`
2.  **Update Odds:** The new odds $O^{(t + 1)}$ is the simple ratio of
    the *total expected nonrespondents* to the *total observed
    respondents* for each $\left( x_{1},y \right)$ cell.
    - **Formula:**
      $$O^{(t + 1)}\left( x_{1}^{*},y^{*} \right) = \frac{M_{y^{*}x_{1}^{*}}^{(t)}}{N_{y^{*}x_{1}^{*}}}$$
    - **Implementation:** `odds_matrix <- m_y_x1_matrix / n_safe`

### Step 2: Convergence

Steps 1.1 and 1.2 are repeated until the change in the `odds_matrix`
(measured by the sum of absolute differences) falls below the tolerance
`tol`.

## Final Estimates and Survey Weights

### Survey Weights

The `exptilt_nonparam.data.frame` function assumes the input `data` is
**already aggregated**. If a `survey.design` object is provided to
`exptilt_nonparam.survey.design`, an adapter (not shown here) must first
be used to create this aggregated table from the microdata. In that
case, $N_{y^{*}x^{*}}$ and $M_{x^{*}}$ represent **weighted sums of
counts** ($\sum d_{i}\ldots$), not simple counts. The EM algorithm and
all derived matrices (`p_hat_matrix`, `odds_matrix`) are then correctly
computed based on these weighted inputs, following the logic of the
paper.

### Final Adjusted Counts

The final output `data_to_return` is the object of primary interest for
analysis. It is constructed by: 1. Calculating the final expected
nonrespondent counts $M_{y^{*}x^{*}}^{(\text{final})}$ using the
converged `odds_matrix`. 2. Adding these expected counts to the original
observed respondent counts $N_{y^{*}x^{*}}$.

- **Definition:**
  `data_to_return[y^*, x^*] = N_{y^*x^*} + M_{y^*x^*}^{(\text{final})}`
- **Implementation:**
  `data_to_return[, outcome_cols] <- n_y_x_matrix_ordered + m_y_x_matrix_ordered`

This final adjusted table represents the completed dataset, where the
“Refusal” counts have been redistributed across the outcome columns
according to the NMAR model.

### Final Proportion ${\widehat{\theta}}_{j}$

The final population proportion for an outcome $j$,
${\widehat{\theta}}_{j}$, is the total (weighted) adjusted count for
outcome $j$ divided by the total (weighted) population. This is
calculated *from* the `data_to_return` object.

- **Formula (Unweighted):**
  $${\widehat{\theta}}_{j} = \frac{\sum\limits_{x^{*}}\left( N_{jx^{*}} + M_{jx^{*}}^{(\text{final})} \right)}{\sum\limits_{y}\sum\limits_{x^{*}}\left( N_{yx^{*}} + M_{yx^{*}}^{(\text{final})} \right)}$$
- **Implementation:** This is precisely what the “Adjusted Proportions”
  calculation in your notebook’s Chunk 5 performs on the
  `data_to_return` object.

This differs from the paper’s Step 3 formula, which is an IPW (Inverse
Probability Weighting) estimator. However, both methods are
asymptotically equivalent, and the “add-and-sum” method (your
`data_to_return`) is a more direct and intuitive application of the EM
algorithm’s goal.
