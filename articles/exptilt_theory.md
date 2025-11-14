# Exponential Tilting Theory

## Overview

This vignette explains the matrix-vectorized implementation of the
exponential tilting (EXPTILT) estimator used in this package. The
exposition focuses on concepts and linear-algebraic operations used in
the implementation (functions such as
`generate_conditional_density_matrix`, `generate_C_matrix`,
`generate_Odds`, and the vectorized `s_function`), and mirrors the
statistical equations that underlie the algorithm (see Riddles et
al. (2016) for the theoretical background).

The method operates on a dataset split into respondents and
non-respondents. Let $$n = n_{0} + n_{1}$$ where $n_{1}$ is the number
of respondents (observed $Y$) and $n_{0}$ is the number of
non-respondents (missing $Y$). The algorithm always treats these two
groups separately: the conditional-density model is fit on respondents
and used to impute support for non-respondents; the missingness model is
fit using covariates that include the candidate $y$ values.

We emphasize two distinct and disjoint sets of covariates throughout,
and we require they have empty intersection:
$$\mathcal{X}_{\text{outcome}} \cap \mathcal{X}_{\text{missingness}} = \varnothing.$$

- covariates_for_outcome (denote $X_{\text{out}}$): variables used to
  model the conditional density
  $f_{1}\left( y \mid X_{\text{out}} \right)$ on respondents;
- covariates_for_missingness (denote $X_{\text{miss}}$): variables used
  in the response probability $\pi\left( X_{\text{miss}},y;\phi \right)$
  (this includes candidate $y$ values when evaluating the missing-data
  expectation).

Distinguishing these sets clearly in notation and code is crucial: the
conditional density model is fit using `covariates_for_outcome` on
respondents only, while the missingness model uses
`covariates_for_missingness` and (importantly) takes $y$ as an
additional predictor when forming $\pi( \cdot ;\phi)$.

## Notation and main objects

Let:

- $n$ be the total number of units; split them into respondents
  (observed $Y$) and non-respondents (missing $Y$). Denote indices with
  $i$ for non-respondents and $k$ (or $j$) for respondent outcomes used
  as support.
- $\{ y_{j}\}_{j = 1}^{n_{1}}$ be the vector of observed outcome values
  (respondents).
- $X^{\text{out}}$ denote the matrix of covariates_for_outcome for
  respondents.
- $X^{\text{un}}$ denote the matrix of covariates_for_outcome for
  non-respondents.
- $X^{\text{miss,obs}}$ and $X^{\text{miss,un}}$ denote the
  covariates_for_missingness for respondents and non-respondents,
  respectively (we shorten these to $X_{\text{miss}}^{\text{obs}}$ and
  $X_{\text{miss}}^{\text{un}}$ when needed).
- $\pi\left( x_{\text{miss}},y;\phi \right)$ the response probability
  (no explicit link notation here): for a given covariate row used in
  the missingness model and a candidate $y$, $\pi( \cdot ;\phi)$ denotes
  the probability of response under parameter $\phi$.

We build three central matrices:

1.  Conditional-density matrix $F$ (denoted in code as
    `f_matrix_nieobs`):

    - Size: $n_{0} \times n_{1}$ where $n_{0}$ is the number of
      non-respondents and $n_{1}$ is the number of distinct respondent
      $y$ values used as support.
    - Entries:
      $$F_{ij} = f_{1}\left( y_{j} \mid x_{i}^{\text{un}};\widehat{\gamma} \right)$$  
      where $\widehat{\gamma}$ is the (estimated) parameter vector of
      the conditional-density model fit on respondents. This corresponds
      to the empirical approximation in equation (12):
      $${\widehat{f}}_{1}\left( y_{j} \right) \propto \sum\limits_{k:\,\delta_{k} = 1}f_{1}\left( y_{j} \mid x_{k};\widehat{\gamma} \right).$$

2.  Column-normalizer vector $C$ (denoted in code as `C_matrix_nieobs`):

    - Size: $n_{1} \times 1$ (a column vector).
    - Entries: the column-sums of the conditional densities evaluated at
      respondent covariates:
      $$C_{j} = C\left( y_{j};\widehat{\gamma} \right) = \sum\limits_{k:\,\delta_{k} = 1}f_{1}\left( y_{j} \mid x_{k}^{\text{obs}};\widehat{\gamma} \right).$$
      Conceptually this is the denominator that appears when fractional
      weights are formed (see below).

3.  Odds matrix $O(\phi)$ (constructed by `generate_Odds`):

    - Size: $n_{0} \times n_{1}$.

- Entries (for non-respondent $i$ and candidate $y_{j}$):
  $$O_{ij}(\phi) = \frac{1 - \pi\left( x_{i}^{\text{miss,un}},y_{j};\phi \right)}{\pi\left( x_{i}^{\text{miss,un}},y_{j};\phi \right)}$$
  The implementation exploits the separability of the linear predictor
  in the parameters:
  $$\eta_{ij} = \beta_{0} + X_{i}^{\text{miss,un}}\beta_{\text{miss}} + \beta_{y}y_{j},$$
  and uses [`outer()`](https://rdrr.io/r/base/outer.html) to form the
  $n_{0} \times n_{1}$ matrix efficiently.

Using these objects we form a non-normalized weight matrix
$$U_{ij}(\phi) = O_{ij}(\phi) \cdot F_{ij}$$ and then a normalized
fractional-weight matrix
$$W_{ij}(\phi) = \frac{U_{ij}(\phi)}{C_{j}} = \frac{O_{ij}(\phi)\, f_{1}\left( y_{j} \mid x_{i}^{\text{un}};\widehat{\gamma} \right)}{\sum\limits_{k:\,\delta_{k} = 1}f_{1}\left( y_{j} \mid x_{k}^{\text{obs}};\widehat{\gamma} \right)}.$$

These $W_{ij}$ match the weights appearing in the theoretical mean score
approximation (equation (15) in the notes): they are the fractional
contribution of each imputed $(i,j)$ pair to the conditional expectation
for unit $i$.

## The vectorized score S_2 and matrix algebra

Recall the mean-score approximation that leads to the estimating
equation used to solve for $\phi$ (cf. equation (13) in the notes):
$$S_{2}\left( \phi;{\widehat{\phi}}^{(t)},\widehat{\gamma} \right) = \sum\limits_{r = 1}^{n}\lbrack\delta_{r}s\left( \phi;\delta_{r},x_{r}^{\text{miss}},x_{r}^{\text{out}},y_{r} \right) + \left( 1 - \delta_{r} \right){\widetilde{E}}_{0}\{ s\left( \phi;\delta,x_{r}^{\text{miss}},x_{r}^{\text{out}},Y \right) \mid x_{r}^{\text{out}};{\widehat{\phi}}^{(t)},\widehat{\gamma}\}\rbrack = 0.$$

It is convenient to decompose this as the sum of observed and unobserved
contributions:
$$S_{2}(\phi) = S_{\text{obs}}(\phi) + S_{\text{un}}(\phi),$$ where
$$S_{\text{obs}}(\phi) = \sum\limits_{r:\,\delta_{r} = 1}s\left( \phi;\delta = 1,x_{r}^{\text{miss}},x_{r}^{\text{out}},y_{r} \right)$$
is the score contribution from respondents, and the missing-unit
contribution is approximated by the discrete-support expectation
$$S_{\text{un}}(\phi) \approx \sum\limits_{i:\,\delta_{i} = 0}\sum\limits_{j = 1}^{n_{1}}W_{ij}(\phi)\, s\left( \phi;\delta = 0,x_{i}^{\text{miss}},x_{i}^{\text{out}},y_{j} \right).$$
The remainder of this section explains how the double sum in
$S_{\text{un}}(\phi)$ is computed via matrix operations.

The crucial observation for vectorization is that the inner conditional
expectation for a non-respondent $i$ is approximated by a weighted
finite-sum over the respondent support $\{ y_{j}\}$:
$${\widetilde{E}}_{0}\{ s\left( \phi;\delta,x_{i}^{\text{miss}},x_{i}^{\text{out}},Y \right) \mid x_{i}^{\text{out}}\} \approx \sum\limits_{j = 1}^{n_{1}}W_{ij}(\phi)\, s\left( \phi;\delta = 0,x_{i}^{\text{miss}},x_{i}^{\text{out}},y_{j} \right).$$

Stacking the non-respondent expectations across all non-respondents
gives a single matrix operation. Let

- For each pair $(i,j)$ evaluate the vector-valued score s at
  $\left( \delta = 0,x_{i}^{\text{miss}},x_{i}^{\text{out}},y_{j} \right)$.
  Collect these quickly by exploiting the algebraic factorization of the
  score (see the `s_function` implementation): many parameter-specific
  components are separable in $i$ and $j$ which allows creation of
  low-memory representations.

- Denote by $S_{ij}^{(0)}$ the p-dimensional score vector at $(i,j)$.
  Organize these so that for each parameter index $m$ we have an
  $n_{0} \times n_{1}$ matrix of values
  $\left\lbrack S_{\bullet \bullet}^{(0)} \right\rbrack_{m}$
  (parameter-wise maps). Then the non-respondent contribution to the
  overall score vector is
  $$\sum\limits_{i = 1}^{n_{0}}\sum\limits_{j = 1}^{n_{1}}W_{ij}(\phi)\, S_{ij}^{(0)} = \left\lbrack \,\text{vec}(W \circ \left\lbrack S^{(0)} \right\rbrack_{1})\,,\,\ldots\,,\,\text{vec}(W \circ \left\lbrack S^{(0)} \right\rbrack_{p})\, \right\rbrack^{\top}$$
  where $\circ$ denotes elementwise multiplication and `vec` followed by
  an appropriate collapse (row-sum or column-sum) implements the inner
  summation depending on the parameter’s factorization. In concrete
  computational terms:

- For parameter components that multiply per-row (i.e. depend only on
  $x_{i}$ times a factor that is a function of $y_{j}$) we compute
  elementwise products between $W$ and a factor matrix and then row-sum
  across $j$ to get an $n_{0} \times 1$ contribution per non-respondent,
  and then sum across $i$.

- For intercept-like or column-wise components, a column-sum followed by
  weighted multiplication suffices.

In the implementation this reduces to a sequence of dense matrix
operations and row/column sums rather than explicit loops over the
expanded index set of length $n_{0} \times n_{1}$. This yields large
speed and memory benefits for real datasets.

Concise vectorized S_2 recipe (conceptual):

1.  Build $F$ (size $n_{0} \times n_{1}$) via
    `generate_conditional_density_matrix(model)`.
2.  Build $C$ (size $n_{1}$) via `generate_C_matrix(model)` by summing
    the conditional densities over respondents.
3.  For a candidate $\phi$ compute $O(\phi)$ (size $n_{0} \times n_{1}$)
    via `generate_Odds(model,\phi)`.
4.  Form $W(\phi) = \frac{O(\phi) \circ F}{\mathbf{1}_{n_{0}}C^{\top}}$
    (i.e. divide each column of $O \circ F$ by the corresponding scalar
    $C_{j}$).
5.  Compute observed-score sum:
    $S_{\text{obs}}(\phi) = \sum_{r:\,\delta_{r} = 1}s\left( \phi;\delta = 1,x_{r}^{\text{miss}},x_{r}^{\text{out}},y_{r} \right)$
    (this is small: one score vector per respondent).
6.  Compute non-respondent expected-score: use $W(\phi)$ and
    parameter-wise factor matrices derived from $s( \cdot )$ to compute
    $$S_{\text{un}}(\phi) = \sum\limits_{i = 1}^{n_{0}}\sum\limits_{j = 1}^{n_{1}}W_{ij}(\phi)\, s\left( \phi;\delta = 0,x_{i}^{\text{miss}},x_{i}^{\text{out}},y_{j} \right)$$
    implemented via matrix multiplications and row/column sums.
7.  Return the total score
    $S_{2}(\phi) = S_{\text{obs}}(\phi) + S_{\text{un}}(\phi)$.

The root-finder then searches for $\phi$ such that $S_{2}(\phi) = 0$. In
the package this is implemented by repeatedly forming $O(\phi)$ for the
current candidate $\phi$, computing $W(\phi)$ and $S_{2}(\phi)$, and
letting `nleqslv` perform the iteration.

## Mean estimation and survey weights

After the response-model parameters $\widehat{\phi}$ are obtained the
package reports a mean estimate of the target outcome using an
inverse-probability reweighting form. Let respondents be indexed by
$j = 1,\ldots,n_{1}$; denote by
$\pi_{j} = \pi\left( x_{j}^{\text{miss}},y_{j};\widehat{\phi} \right)$
the fitted response probabilities evaluated at each respondent (here we
write $x_{j}^{\text{miss}}$ for the covariates used in the missingness
model evaluated at respondent $j$). With design weights $w_{j}$ (by
default $w_{j} \equiv 1$ for non-survey data) the point estimator
computed in the implementation is
$$\widehat{\mu} = \frac{\sum\limits_{j = 1}^{n_{1}}w_{j}\, y_{j}/\pi_{j}}{\sum\limits_{j = 1}^{n_{1}}w_{j}/\pi_{j}}.$$

Notes on survey weights and $f_{1}$ fitting:

- The conditional-density fit used to build $F$ can incorporate
  respondent sampling weights when provided (the implementation passes
  respondent weights to the density-fitting routine). Thus
  $f_{1}( \cdot \mid x;\gamma)$ may be a weighted fit when
  `design_weights` or a survey `design` is supplied.
- The mean-estimate (mean-est) also uses design weights $w_{j}$ in both
  numerator and denominator as shown above. In code these are provided
  via `model$design_weights` and default to 1 for non-survey use.

In short: survey weights enter two places — (i) the conditional-density
estimation (if requested) and (ii) the final mean calculation through
the weighted IPW-type ratio in (mean-est).

## Arguments passed to `exptilt` (summary)

The `exptilt` / `exptilt.data.frame` user-facing function accepts a
number of arguments that control model specification and computation;
the most important are:

- `data`: a `data.frame` with outcome and covariates.
- `formula`: a partitioned formula of the form
  `y ~ aux1 + aux2 | miss1 + miss2` where the left part (before `|`)
  lists `covariates_for_outcome` and the right part lists
  `covariates_for_missingness` (the package helper splits these
  automatically).
- `auxiliary_means`: (optional) population or target means for auxiliary
  variables used for scaling.
- `standardize`: logical, whether to standardize features before
  fitting.
- `prob_model_type`: character, either `"logit"` or `"probit"` to select
  the response probability family.
- `y_dens`: choice of conditional-density family; `"auto"`, `"normal"`,
  `"lognormal"`, or `"exponential"`.
- `variance_method`: `"delta"` or `"bootstrap"` for variance estimation.
- `bootstrap_reps`: number of bootstrap replications when bootstrap is
  used.
- `control`: list of control parameters forwarded to the nonlinear
  solver (`nleqslv`).
- `stopping_threshold`: numeric; the sup-norm threshold for early
  stopping of the score (see above).
- `on_failure`: behavior on failure (`"return"` or `"error"`).
- `supress_warnings`: logical to silence certain warnings.
- `design_weights`: optional vector of respondent design weights (or
  full-sample weights which are subset internally).
- `survey_design`: optional `survey.design` object; when provided some
  internal logic uses the survey path.
- `trace_level`: integer controlling verbosity.
- `sample_size`: integer for stratified subsampling used when data are
  large.
- `outcome_label`, `user_formula`: utility arguments used for
  presentation and bookkeeping.

These arguments appear in the `exptilt.data.frame` function signature
and control how the matrices $F$, $C$, and $O$ are built and how the
solver is run.

## Connection to the EM viewpoint

The EM-like update displayed in the theoretical notes (equation (14))
$$\left. {\widehat{\phi}}^{(t + 1)}\leftarrow{\text{solve}\mspace{6mu}}S_{2}\left( \phi \mid {\widehat{\phi}}^{(t)},\widehat{\gamma} \right) = 0 \right.$$
is exactly what the implementation achieves: for a fixed
conditional-density estimate $\widehat{\gamma}$ and a current
${\widehat{\phi}}^{(t)}$, the expectations for the missing units are
approximated by the discrete support on observed $y_{j}$ and the
resulting equation for $\phi$ is solved (via root-finding). The
heavy-lift is performed by the matrix calculus described earlier —
constructing $F$, $C$, $O$ and then computing $W$ and multiplying by the
parameter-wise score factors.

## Stopping criterion (maximum-norm)

Practical optimization requires a convergence criterion. The
implementation uses the maximum absolute component of the vector-valued
score as a stopping rule. Concretely, if the solver produces a score
vector $S_{2}(\phi)$ with
$$\max\limits_{m = 1,\ldots,p}\left| S_{2,m}(\phi) \right| < \varepsilon$$
for a user-specified `stopping_threshold = \varepsilon`, the algorithm
treats this as converged. In the code this is used as an early-exit
inside the target-function passed to the nonlinear solver: when the
score’s sup-norm is below the threshold, a zero vector is returned to
signal convergence and avoid further unnecessary computations.

This choice matches the intuition that the root-finder should stop when
the estimating equations are (componentwise) negligibly small.

## Practical tutorial: from raw data to matrix operations (conceptual steps)

1.  Fit the conditional density $f_{1}(y \mid x;\gamma)$ using
    respondents and `covariates_for_outcome`. This gives a function that
    evaluates $f_{1}(y,x)$ for any $y$ and any covariate row $x$ (used
    for both respondents and non-respondents).

2.  Evaluate the conditional density at the Cartesian product of
    non-respondent covariates and observed respondent $y$ values to form
    $F$ (done by `generate_conditional_density_matrix`). This is the
    empirical imputation support. Think of rows as target
    non-respondents and columns as candidate respondent outcomes.

3.  Evaluate the same conditional density at respondent covariates to
    form the column-normalizer $C$ (done by `generate_C_matrix`) — this
    is simply the column-sum of densities over respondents.

4.  For each trial value $\phi$ of the response-model parameters,
    compute the odds matrix $O(\phi)$ using the separable linear
    predictor and link function (implemented efficiently via
    [`outer()`](https://rdrr.io/r/base/outer.html) in code). Combine $O$
    with $F$ and normalize columns by $C$ to obtain $W(\phi)$.

5.  Use the vectorized `s_function` to obtain parameter-specific factor
    matrices for the non-respondent-imputed scores; multiply
    (elementwise) by $W(\phi)$ and reduce (row/column sums) to compute
    the non-respondent contribution to $S_{2}(\phi)$.

6.  Add the observed-respondent score and use a root-finder
    (e.g. `nleqslv`) to find $\phi$ with $S_{2}(\phi) = 0$. The solver
    may use the maximum-norm stopping threshold described above to exit
    early.

## Why the vectorization matters (practical remarks)

- Memory: the naive expansion to an explicit dataset of size
  $n_{0} \times n_{1}$ would store duplicated covariate rows and blow up
  memory. The implemention exploits separability (intercept +
  $X_{i}^{\text{miss}}\beta_{\text{miss}}$ + $\beta_{y}y_{j}$) and
  vectorized R primitives (`outer`, matrix multiplications, column/row
  sums) to avoid large temporary allocations.
- Speed: elementwise operations on dense matrices plus BLAS-accelerated
  matrix multiplication are much faster than interpreted loops in R for
  typical dataset sizes.
- Clarity: organizing logic as the three matrices $F$, $C$, and $O$,
  followed by elementwise combination and reductions, makes the
  relationship between the statistical approximation and the
  implementation transparent and easier to reason about.

## Closing notes and references

This vignette mapped the implementation functions to the math in the
theoretical notes and showed how the EM-like mean-score step reduces to
a small set of matrix operations. The implementation follows the ideas
described in Riddles et al. (2016) for exponential tilting under NMAR:
fit the conditional density on respondents, approximate the missing-data
expectation by a finite sum over observed $y$ values, and solve the
resulting estimating equations for the missingness-model parameters.

If you would like, I can produce a minimal worked numeric example (R
code) that builds the matrices and demonstrates the matrix
multiplications used here — implemented as small, readable functions for
teaching purposes — and place it in this vignette as an optional
appendix.
