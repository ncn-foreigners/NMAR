# Empirical Likelihood Theory

This document explains every mathematical object, equation, and
derivation behind the empirical likelihood (EL) estimator implemented in
the `nmar` package, and maps each concept to code. It covers both
data-frame (IID) and survey design use cases, handles arbitrary numbers
of response-model and auxiliary covariates, and supports both logit and
probit response families.

## Notation

### Units

- $i = 1,\ldots,n$ index respondents (those with observed $Y$)
- $R_{i} \in \{ 0,1\}$ is the response indicator; we work on observed
  subset $R_{i} = 1$

### Data

- **Outcome**: $Y_{i}$ (observed when $R_{i} = 1$; missing otherwise)
- **Response covariates**: row vector $Z_{i} \in {\mathbb{R}}^{K}$, from
  `model.matrix` of the response RHS
- **Auxiliary covariates**: row vector $X_{i} \in {\mathbb{R}}^{L}$
  (possibly $L = 0$), from auxiliary RHS (no intercept)
- **Population auxiliary means**: $\mu_{x} \in {\mathbb{R}}^{L}$, known;
  names match columns of $X$

Mapping to code: - `response_model_matrix` corresponds to $Z$ (includes
an intercept column) - `auxiliary_matrix` corresponds to $X$ (no
intercept); we center it in code as $X - \mu_{x}$

### Response Model (Family functions)

- **Linear predictor**: $\eta_{i} = Z_{i}\,\beta$
- **Response probability**:
  $w_{i} \equiv g\left( \eta_{i} \right) = {linkinv}\left( \eta_{i} \right)$
- **First derivative**:
  $\frac{dw}{d\eta}\left( \eta_{i} \right) = \mu_{\eta,i} = {mu.eta}\left( \eta_{i} \right)$
- **Second derivative**:
  $\frac{d^{2}w}{d\eta^{2}}\left( \eta_{i} \right) = {d2mu.deta2}\left( \eta_{i} \right)$
  Here `linkinv`, `mu.eta`, and `d2mu.deta2` refer to the chosen
  response family (logit or probit). We follow the paper’s $w_{i}$
  notation for the response probability and reserve $p_{i}^{\text{EL}}$
  for empirical-likelihood weights.

### Weight Re-parameterization

- $W \in (0,1)$ nuisance scalar; we parameterize via
  $z = \text{logit}(W)$ for stability and set $W = \text{plogis}(z)$
- $\lambda_{W} \in {\mathbb{R}}$ and $\lambda_{x} \in {\mathbb{R}}^{L}$
  are EL Lagrange multipliers for constraints; collected together in
  $\theta$

### EL Weights

- **Denominator**:
  $D_{i} = 1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{T}\lambda_{x}$
- **Base sampling weights**: $a_{i} = 1$ (IID) or $a_{i} =$ survey base
  weight for respondent $i$
- **EL weights for respondents**:
  $p_{i}^{\text{EL}} \propto a_{i}/D_{i}$ (proportionality normalized by
  totals below)

### Estimator

- $\widehat{Y} = \sum p_{i}^{\text{EL}}Y_{i}/\sum p_{i}^{\text{EL}}$

### Notation at a Glance

| Symbol                  | Meaning                                                                                   |
|-------------------------|-------------------------------------------------------------------------------------------|
| $i$                     | Respondent index (rows with observed $Y$)                                                 |
| $Y_{i}$                 | Outcome for unit $i$ (observed if $R_{i} = 1$)                                            |
| $Z_{i}$                 | Row of response design matrix (includes intercept)                                        |
| $X_{i}$                 | Row of auxiliary design (no intercept)                                                    |
| $\mu_{x}$               | Known population means of auxiliaries (vector)                                            |
| $\beta$                 | Response-model coefficients                                                               |
| $\eta_{i} = Z_{i}\beta$ | Linear predictor for response model                                                       |
| $w_{i}$                 | ${linkinv}\left( \eta_{i} \right)$ (logit: $plogis$; probit: $\Phi$)                      |
| $\mu_{\eta,i}$          | $\frac{dw_{i}}{d\eta_{i}}$                                                                |
| $\lambda_{W}$           | Multiplier for the $W$-constraint $\sum\left( w_{i} - W \right)/D_{i} = 0$                |
| $\lambda_{x}$           | Multipliers for auxiliary constraints $\sum\left( X_{i} - \mu_{x} \right)/D_{i} = 0$      |
| $D_{i}$                 | $1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{T}\lambda_{x}$ |
| $a_{i}$                 | Base weight (IID: 1; survey: design weight)                                               |
| $p_{i}^{EL}$            | Empirical-likelihood weight $\propto a_{i}/D_{i}$                                         |
| $\widehat{Y}$           | $\sum p_{i}^{EL}Y_{i}/\sum p_{i}^{EL}$                                                    |
| $F(\theta)$             | Estimating system (beta, W, and auxiliary equations)                                      |
| $A$                     | Jacobian $\partial F/\partial\theta$                                                      |

### Engines

- Family: “logit” (default) or “probit”; both use the log-likelihood
  score w.r.t. $\eta$:
  $s_{i} = \partial\log w_{i}/\partial\eta_{i} = \mu_{\eta,i}/w_{i}$
  (for respondents, $\delta_{i} = 1$)
- Scaling: optional standardization of design matrices and $\mu_{x}$ via
  nmar_scaling_recipe

## From Paper to Implementation: Core Ideas

The paper (Qin-Leung-Shao, JASA 2002) sets EL under nonignorable
response using:

- **Empirical likelihood weights** for respondents that satisfy:
  - Zero-sum residual:
    $\sum p_{i}^{\text{EL}}\left( w_{i} - W \right) = 0$
  - Auxiliary moments:
    $\sum p_{i}^{\text{EL}}\left( X_{i} - \mu_{x} \right) = 0$
- A **response model probability** $w_{i} = g\left( \eta_{i} \right)$,
  $\eta_{i} = Z_{i}\,\beta$

In our code, we adopt the same EL structure and estimating equations. We
extend it to arbitrary $Z$ and $X$, and to survey designs. For
uncertainty, we provide bootstrap variance (IID resampling and survey
replicate-weight bootstrap). Because $\widehat{Y}$ is a ratio-of-weights
estimator, any common normalization of
$p_{i}^{\text{EL}} \propto a_{i}/D_{i}$ cancels in $\widehat{Y}$; only
relative weights matter (the KKT multipliers $\lambda$ enforce the
constraints; normalization affects only a common scale that vanishes in
the ratio).

### Closed-form $\lambda_{W}$ (QLS Eq. 10, IID path)

QLS derive $\lambda_{W} = (N/n - 1)/(1 - W)$ when $a_{i} \equiv 1$ and
$n$ counts respondents. In our **IID (data-frame) path** we reuse this
relation, generalized to allow arbitrary base weights $a_{i}$:

- Let $n_{\text{resp\_weighted}} = \sum_{i:R_{i} = 1}a_{i}$ be the
  respondent-weighted total.
- Let $N_{\text{pop}}$ be the analysis-scale population total (either
  $n$ or a user-supplied `n_total`).

We set

$$\lambda_{W} = \frac{N_{\text{pop}}/n_{\text{resp\_weighted}} - 1}{1 - W},$$

which reduces to the original QLS formula when $a_{i} \equiv 1$. This
closed form is used **only** in the IID path to profile out
$\lambda_{W}$; for complex survey designs we instead treat $\lambda_{W}$
as a free parameter and solve for it jointly with
$\left( \beta,W,\lambda_{x} \right)$ via the design-weighted system
described in the survey extension section.

### Guarding and Numerical Stability

We solve the stacked system with a consistent guarding policy across
equations, Jacobian, and post-solution:

- Cap $\eta$:
  $\left. \eta\leftarrow\max\{\min\left( \eta,\,\eta_{\max} \right),\, - \eta_{\max}\} \right.$
- Compute $w = {linkinv}(\eta)$; clip $w$ to
  $\left\lbrack 10^{- 12},1 - 10^{- 12} \right\rbrack$ when used in
  ratios
- Guard denominators:
  $\left. D_{i}\leftarrow\max\{ D_{i},\,\delta\} \right.$ with a small
  $\delta > 0$
- In the Jacobian, multiply terms involving
  $\partial\left( 1/D_{i} \right)/\partial \cdot$ by
  $\mathbb{1}\{ D_{i}^{\text{raw}} > \delta\}$ so the analytic Jacobian
  matches the piecewise-smooth equations being solved

For the probit link,
$s_{i}(\eta) = \partial\log w/\partial\eta = \phi(\eta)/\Phi(\eta)$
(Mills ratio) is computed in the log domain for stability; its
derivative is $\frac{ds_{i}}{d\eta} = - \eta\, s_{i} - s_{i}^{2}$.

### Equation Crosswalk (QLS 2002 -\> This Vignette/Code)

- QLS (5): Discrete mass form for $p_{i}$ with two multipliers -\> Our
  $D_{i} = 1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{T}\lambda_{x}$
  and $p_{i}^{\text{EL}} \propto a_{i}/D_{i}$.
- QLS (7): $\sum\frac{x_{i} - \bar{x}}{1 + \cdots} = 0$ -\> Our
  auxiliary constraints
  $\sum a_{i}\left( X_{i} - \mu_{x} \right)/D_{i} = 0$.
- QLS (8): $\sum\frac{w_{i} - W}{1 + \cdots} = 0$ -\> Our $W$-equation
  $\sum a_{i}\left( w_{i} - W \right)/D_{i} = 0$.
- QLS (10): ${\widehat{\lambda}}_{2} = (N/n - 1)/(1 - W)$ -\> In the IID
  path we set
  $\lambda_{W} = \left( \left( N_{\text{pop}}/n_{\text{resp\_weighted}} \right) - 1 \right)/(1 - W)$;
  in the survey path $\lambda_{W}$ is solved from the additional linkage
  equation $g_{W}^{(1)}$.
- Estimator $\widehat{Y}$ in QLS -\> Our ratio
  $\widehat{Y} = \sum p_{i}^{EL}Y_{i}/\sum p_{i}^{EL}$ using
  $p_{i}^{EL} \propto a_{i}/D_{i}$.

### Likelihood and Profiling (sketch)

The paper’s semiparametric likelihood (their Eq. (2)) combines the
response mechanism $w_{i} = g\left( \eta_{i} \right)$ with the
nonparametric distribution $F$ of $(Y,X)$:

$$\mathcal{L}(\beta,W,F)\; \propto \;\prod\limits_{i = 1}^{n}w\left( Y_{i},X_{i};\beta \right)\, dF\left( Y_{i},X_{i} \right)\; \times \;(1 - W)^{N - n},$$

subject to (i) $\int dF = 1$, (ii) $\int X\, dF = \mu_{x}$, and (iii)
$\int w(Y,X;\beta)\, dF = W$. Discretizing $F$ at observed respondents
by assigning unknown masses $p_{i}$ and introducing multipliers
$\lambda$, the KKT conditions yield the familiar EL weight form with
denominator

$$D_{i}\; = \; 1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{\top}\lambda_{x},$$

and, with base weights $a_{i}$, the working weights are
$p_{i}^{\text{EL}} \propto a_{i}/D_{i}$.

Remark on conditioning: QLS’s Eq. (2) writes the first product as
$\prod_{i}\left\lbrack \, w\left( y_{i},x_{i};\beta \right)\, dF\left( y_{i},x_{i} \right)/W\, \right\rbrack$
so that it explicitly represents the likelihood of
$\left( Y_{i},X_{i} \right)$ conditional on $R_{i} = 1$. Multiplying by
the binomial term $W^{n}(1 - W)^{N - n}$ yields the same overall
likelihood as above because the $W^{- n}$ in the first factor cancels
the $W^{n}$ in the second. Both factorizations lead to the same
estimating equations and the same profiled log-likelihood form used
subsequently in QLS after introducing the multipliers.

### KKT and Denominator (details)

Introducing Lagrange multipliers
$\left( \lambda_{0},\lambda_{x},\lambda_{W} \right)$ for these
constraints and profiling the $p_{i}$’s gives the KKT stationarity
conditions (in a survey-weighted setting, with base weights $a_{i}$
acting as multiplicities in the empirical likelihood)

$$\frac{\partial}{\partial p_{i}}\lbrack\sum\limits_{j}a_{j}\,\log p_{j} - \lambda_{0}\left( \sum\limits_{j}p_{j} - 1 \right) - \lambda_{x}^{T}\sum\limits_{j}p_{j}\left( X_{j} - \mu_{x} \right) - \lambda_{W}\sum\limits_{j}p_{j}\left( w_{j} - W \right)\rbrack = 0,$$

which solve to

$$p_{i}\; \propto \;\frac{1}{\, 1 + \lambda_{x}^{T}\left( X_{i} - \mu_{x} \right) + \lambda_{W}\left( w_{i} - W \right)\,}\; \equiv \;\frac{1}{D_{i}}.$$

Normalizing to enforce $\sum p_{i} = 1$ yields
$p_{i} = \frac{D_{i}^{- 1}}{\sum_{j}D_{j}^{- 1}}$. In the presence of
base sampling weights $a_{i}$ (survey designs), the same derivation
gives the natural generalization

$$p_{i}^{\text{EL}}\; \propto \;\frac{a_{i}}{D_{i}}\quad\text{with}\quad D_{i} = 1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{T}\lambda_{x}.$$

This is exactly the working form used in our estimator. With survey base
weights, the working masses are proportional to $a_{i}/D_{i}$, so
$p_{i}^{\text{EL}} \propto a_{i}/D_{i}$. The EL weights
$p_{i}^{\text{EL}}$ are then used to build the mean estimator

$$\widehat{Y}\; = \;\frac{\sum\limits_{i}p_{i}^{\text{EL}}Y_{i}}{\sum\limits_{i}p_{i}^{\text{EL}}}.$$

The remaining unknowns $\left( \beta,W,\lambda_{x} \right)$ are
determined by the estimating equations below.

### Clarification: Relationship Between $W$ and $\lambda_{W}$

In the IID (data-frame) path, the EL multiplier for the probability
constraint is expressed as

$$\lambda_{W} = \frac{C}{1 - W},\quad{\text{with}\mspace{6mu}}C = \frac{N_{\text{pop}}}{n_{\text{resp\_weighted}}} - 1{\mspace{6mu}\text{and}\mspace{6mu}}W = \text{plogis}(z).$$

Intuition: In the EL KKT system, the constraint
$\sum p_{i}^{\text{EL}}\left( w_{i} - W \right) = 0$ sits alongside
normalization and (optionally) auxiliary constraints. Incorporating base
weights $a_{i}$ and the ratio between population and respondent totals
induces a scaling of the multiplier linked to the mass constraint.
Writing $\lambda_{W}$ in this scaled form keeps the parameter on a
numerically stable scale and lets the derivative structure (with respect
to $z$ via $W$) be handled cleanly. This is consistent with the EL
structure when the baseline mass is $n_{\text{resp\_weighted}}$ and the
“full population” target is $N_{\text{pop}}$, and it is exactly what the
IID code path uses to match the normalization implied by base weights.

Derivation sketch (KKT, IID case): The discretized semiparametric
likelihood (QLS, 2002) maximizes, over the unknown masses $\{ p_{i}\}$
at observed points and over $(\beta,W)$,

$$\ell\left( \beta,W,\lambda_{x},\lambda_{W} \right)\; = \;\sum\limits_{i = 1}^{n}\log w_{i}(\beta)\; + \;\left( N_{\text{pop}} - n_{\text{resp\_weighted}} \right)\log(1 - W)\; - \;\sum\limits_{i = 1}^{n}\log\!(1 + \left( X_{i} - \mu_{x} \right)^{\top}\lambda_{x} + \lambda_{W}\left( w_{i} - W \right)),$$

subject to the normalization and moment constraints that generate the EL
denominator. For the weighted-EL variant we work with unnormalized
respondent weights proportional to $a_{i}/D_{i}$; choosing the
conventional normalization

$$\sum\limits_{i = 1}^{n}\frac{a_{i}}{D_{i}}\, = \, n_{\text{resp\_weighted}} \equiv \sum\limits_{i = 1}^{n}a_{i}$$

recovers the same estimating system (and any common normalization
cancels in the ratio estimator
$\widehat{Y} = \sum p_{i}Y_{i}/\sum p_{i}$). Taking derivatives (KKT
conditions) and using that $\partial/\partial W$ of the second and third
terms produces opposing contributions, one obtains the system equivalent
to QLS (7)-(10). In particular, the first-order condition with respect
to the multiplier associated with the $W$-constraint yields, together
with the derivative with respect to $W$, the closed form

$$\lambda_{W}\; = \;\frac{\frac{N_{\text{pop}}}{n_{\text{resp\_weighted}}} - 1}{1 - W}\; = \;\frac{C}{1 - W},$$

which coincides with QLS (10) when $a_{i} \equiv 1$. This closed-form
relationship is used in the IID EL implementation to profile out
$\lambda_{W}$. In the **survey-design path**, by contrast, $\lambda_{W}$
is treated as an explicit unknown and the linkage between $W$ and
$\lambda_{W}$ is enforced through the additional equation

$$g_{W}^{(1)}\left( \beta,W,\lambda_{W},\lambda_{x} \right)\; = \;\frac{T_{0}}{1 - W} - \lambda_{W}\sum\limits_{i \in R}\frac{d_{i}}{D_{i}} = 0,$$

where $T_{0} = N_{pop} - \sum_{i \in R}d_{i}$ and $d_{i}$ are the design
weights. At the QLS simple-random-sampling limit (equal weights, no
auxiliaries) this system reduces to the same closed-form relation.

## Estimating Equations

**Unknown parameters**: $\beta \in {\mathbb{R}}^{K}$,
$z \in {\mathbb{R}}$ (for $W = \text{plogis}(z)$),
$\lambda_{x} \in {\mathbb{R}}^{L}$; define
$\theta = \left( \beta,z,\lambda_{x} \right)$.

Define $w_{i} = {linkinv}\left( \eta_{i} \right)$ and
$\mu_{\eta,i} = \frac{dw}{d\eta}\left( \eta_{i} \right)$ (denoted
`mu.eta(eta_i)` in code).

In the **IID (data-frame) path** all base weights are $a_{i} \equiv 1$,
so we can use the closed-form Qin-Leung-Shao (QLS) relation between $W$
and the EL multiplier for the response constraint. Writing
$$C = \frac{N_{\text{pop}}}{n_{\text{resp\_weighted}}} - 1,\qquad n_{\text{resp\_weighted}} = \sum\limits_{i}a_{i},\qquad N_{\text{pop}} = n,$$
QLS show that
$$\lambda_{W} = \frac{C}{1 - W} = \frac{N_{\text{pop}}/n_{\text{resp\_weighted}} - 1}{1 - W}.$$
Our IID implementation follows this and *profiles out* $\lambda_{W}$:
the unknowns for the Newton solver are
$\left( \beta,z,\lambda_{x} \right)$.

In the **survey path**, base weights are general design weights
$a_{i} = d_{i}$ and the corresponding QLS-style relation no longer has a
simple closed form. In that case we treat $\lambda_{W}$ as an additional
free parameter and include a separate equation linking $\lambda_{W}$ and
$W$ (see the “Survey extension” section below).

**Denominator**:
$D_{i} = 1 + \lambda_{W}\left( w_{i} - W \right) + \left( X_{i} - \mu_{x} \right)^{T}\lambda_{x}$,
with $D_{i} \geq \epsilon$ enforced numerically.

Define the score term $s_{i} = \mu_{\eta,i}/w_{i}$ (the unit-level
contribution to the log-likelihood score with respect to $\eta$). For
logit, $s_{i} = 1 - w_{i}$; for probit, $s_{i}$ behaves like
$\phi\left( \eta_{i} \right)/\Phi\left( \eta_{i} \right)$ when $w_{i}$
is bounded away from 0 via clipping (as implemented).

Intuition (why this score appears): for each respondent we observe
$R_{i} = 1$, so the Bernoulli log-likelihood contribution of the
response model is $\log w_{i}\left( \eta_{i} \right)$. Differentiating
w.r.t. the linear predictor gives
$$\frac{\partial}{\partial\eta_{i}}\log w_{i}\left( \eta_{i} \right)\, = \,\frac{1}{w_{i}}\,\frac{dw_{i}}{d\eta_{i}}\, = \,\frac{\mu_{\eta,i}}{w_{i}}\; \equiv \; s_{i}.$$

Thus $s_{i}$ measures the local sensitivity of the observed-response
likelihood to $\eta_{i}$. In the logit family,
$\mu_{\eta,i} = w_{i}\left( 1 - w_{i} \right)$ so
$s_{i} = 1 - w_{i}$-the familiar residual-like term; in the probit
family,
$s_{i} = \phi\left( \eta_{i} \right)/\Phi\left( \eta_{i} \right)$, the
(inverse) Mills ratio. The EL $\beta$-equations balance this likelihood
score against the EL penalty term $\lambda_{W}\,\mu_{\eta,i}/D_{i}$,
enforcing the calibration constraints while fitting the response model.

### The System of Estimating Equations $F(\theta) = 0$

**$\beta$-equations** ($K$ equations):
$$\sum a_{i}Z_{i}\left\lbrack s_{i} - \lambda_{W}\mu_{\eta,i}/D_{i} \right\rbrack = 0$$

**W-equation** (1 equation):
$$\sum a_{i}\left( w_{i} - W \right)/D_{i} = 0$$

**Auxiliary constraints** ($L$ equations):
$$\sum a_{i}\left( X_{i} - \mu_{x} \right)/D_{i} = 0$$

These are exactly how `el_build_equation_system` constructs the function
in code (`src_dev/engines/el/impl/equations.R`).

Intuition: the $\beta$-equations equate the score of the respondent
log-likelihood with the EL penalty term $\lambda_{W}\mu_{\eta,i}/D_{i}$;
the $W$-equation centers the modeled response probabilities around the
unconditional mean $W$ under the EL weights; the auxiliary equations
calibrate the centered auxiliaries to zero mean under the EL weights.

### Remarks

- For logit and probit, $s_{i}$ is the log-likelihood score
  $\partial\log w_{i}/\partial\eta_{i} = \mu_{\eta,i}/w_{i}$ (equals
  $1 - w_{i}$ for logit; behaves like $\phi/\Phi$ for probit when
  $w_{i}$ is clipped away from 0). This follows the paper’s MLE
  derivation; EL constraints supply the nonparametric part.

## Survey extension: design-weighted QLS system

The original QLS paper derives these equations under simple random
sampling, where each respondent has equal weight. In practice we often
work with complex survey designs and design weights
$d_{i} \approx 1/\pi_{i}$, where $\pi_{i}$ is the inclusion probability
for unit $i$. In our implementation we extend the QLS system using a
**design-weighted empirical likelihood**:

- For respondents $i \in R$ we use base weights $a_{i} = d_{i}$.
- We approximate the unknown distribution of $(Y,X)$ by a discrete
  measure
  $$F_{\text{EL}}(A) = \sum\limits_{i \in R}d_{i}p_{i}\,\mathbf{1}\{\left( y_{i},x_{i} \right) \in A\},$$
  where $p_{i} \geq 0$ and $\sum_{i}d_{i}p_{i} = 1$.
- Expectations under $F$ are represented by design-weighted sums
  $\sum_{i}d_{i}p_{i}( \cdot )$.

We impose the following constraints, which are the design-weighted
analogues of QLS (3):

- Normalization: $$\sum\limits_{i \in R}d_{i}p_{i} = 1.$$
- Response-rate constraint:
  $$\sum\limits_{i \in R}d_{i}p_{i}(w_{i}(\theta) - W) = 0.$$
- Auxiliary constraints (vector case):
  $$\sum\limits_{i \in R}d_{i}p_{i}\left( X_{i} - \mu_{x} \right) = 0.$$

Maximizing the design-weighted pseudo-likelihood under these constraints
yields EL weights of the same tilted form as in QLS:
$$p_{i}\left( \theta,W,\lambda_{x},\lambda_{W} \right)\; \propto \;\frac{1}{D_{i}},\qquad D_{i} = 1 + \lambda_{W}(w_{i}(\theta) - W) + \left( X_{i} - \mu_{x} \right)^{\top}\lambda_{x},$$
with the proportionality constant chosen such that
$\sum_{i}d_{i}p_{i} = 1$. In our implementation the **unnormalized EL
masses** are $$m_{i} = \frac{d_{i}}{D_{i}},$$ and the probability-scale
weights are $p_{i}^{EL} = m_{i}/\sum_{j}m_{j}$.

The corresponding design-weighted QLS estimating system in
$\left( \beta,W,\lambda_{W},\lambda_{x} \right)$ can be written as:

- Auxiliary block:
  $$g_{x}\left( \beta,W,\lambda_{W},\lambda_{x} \right) = \sum\limits_{i \in R}d_{i}\frac{X_{i} - \mu_{x}}{D_{i}} = 0.$$
- Response-rate constraint:
  $$g_{W}^{(2)}\left( \beta,W,\lambda_{W},\lambda_{x} \right) = \sum\limits_{i \in R}d_{i}\frac{w_{i}(\beta) - W}{D_{i}} = 0.$$
- Score equations for $\beta$:
  $$g_{\beta}\left( \beta,W,\lambda_{W},\lambda_{x} \right) = \sum\limits_{i \in R}d_{i}\left\lbrack \frac{\partial\log w_{i}(\beta)}{\partial\beta} - \lambda_{W}\frac{1}{D_{i}}\frac{\partial w_{i}(\beta)}{\partial\beta} \right\rbrack = 0.$$
- Linkage between $\lambda_{W}$ and the nonrespondent total:
  $$g_{W}^{(1)}\left( \beta,W,\lambda_{W},\lambda_{x} \right) = \frac{T_{0}}{1 - W} - \lambda_{W}\sum\limits_{i \in R}\frac{d_{i}}{D_{i}} = 0,$$
  where $T_{0} = N_{pop} - \sum_{i \in R}d_{i}$ on the analysis scale.

In code this system is implemented by
[`el_build_equation_system_survey()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_equation_system_survey.md)
in `src_dev/engines/el/impl/equations.R`. The parameter vector is
$$\theta_{\text{survey}} = \left( \beta,z,\lambda_{W},\lambda_{x} \right),\qquad z = \operatorname{logit}(W),$$
and the solver treats $\lambda_{W}$ as an explicit unknown. When all
design weights are equal and $N_{\text{pop}}$ and the respondent count
match the simple random sampling setup, this system reduces exactly to
the original QLS equations (6)-(10).

For survey designs we build an analytic Jacobian for this
design-weighted system whenever the response family supplies a second
derivative `d2mu.deta2` (logit and probit). The Jacobian structure
mirrors the IID case but with the expanded parameter vector
$\left( \beta,z,\lambda_{W},\lambda_{x} \right)$ and the additional
blocks for $g_{W}^{(1)}$ and $g_{W}^{(2)}$. When analytic derivatives
are not available, `nleqslv` falls back to numeric/Broyden Jacobians.

### Wu-style strata augmentation (survey designs)

For some stratified designs, especially when the NMAR mechanism varies
strongly across strata, it is important that the EL weights preserve the
**stratum composition** implied by the survey design. Following ideas
from Wu-style calibration, we augment the auxiliary vector with stratum
indicators when a `survey.design` object is provided:

- Recover a strata factor from the design (using its `strata=` call).
- Build dummy variables for strata (dropping one reference level).
- Compute stratum totals $N_{h}$ on the analysis scale from the design
  weights and convert to stratum shares $W_{h} = N_{h}/N_{pop}$.
- Append these stratum dummies to the auxiliary matrix $X$ and their
  targets $W_{h}$ to the auxiliary means.

The EL constraints then include additional terms of the form
$$\sum\limits_{i \in R}d_{i}p_{i}(\mathbf{1}\{ H_{i} = h\} - W_{h}) = 0$$
for each nonreference stratum $h$. This forces the EL weights to
reproduce the design-implied stratum shares while still adjusting within
strata for NMAR. In the implementation this augmentation is performed in
the survey entry point (`src_dev/engines/el/impl/survey.R`) before
auxiliary means are resolved, and the resulting augmented auxiliaries
flow through to
[`el_build_equation_system()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_equation_system.md)
or
[`el_build_equation_system_survey()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_equation_system_survey.md)
depending on the data type. The behavior is controlled by the logical
`strata_augmentation` argument of
[`el_engine()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md)
(default `TRUE`); it has an effect only when `data` is a `survey.design`
with defined strata.

Our survey EL implementation should be viewed as a design-weighted
analogue of QLS, informed by the pseudo empirical likelihood literature
(Chen and Sitter 1999; Wu 2005), rather than a verbatim implementation
of any single paper.

## Analytic Jacobian ($A$ Matrix, IID path)

For the IID (data-frame) path we differentiate $F(\theta) = 0$ with
respect to $\theta = \left( \beta,z,\lambda_{x} \right)$. Let:

- $\eta_{i} = Z_{i}\beta$,
  $w_{i} = \text{linkinv}\left( \eta_{i} \right)$,
  $\mu_{\eta,i} = \frac{dw}{d\eta}\left( \eta_{i} \right)$,
  $\mu''_{i} = \frac{d^{2}w}{d\eta^{2}}\left( \eta_{i} \right)$
- $W = \text{plogis}(z)$, $\frac{dW}{dz} = W(1 - W)$
- $\lambda_{W} = \frac{C}{1 - W}$, so
  $\frac{d\lambda_{W}}{dW} = \frac{C}{(1 - W)^{2}}$ and
  $\frac{d\lambda_{W}}{dz} = \frac{d\lambda_{W}}{dW} \cdot \frac{dW}{dz}$
- $X_{\text{centered},i} = X_{i} - \mu_{x}$

### Intermediate Derivatives

- $\left. s_{i} = \mu_{\eta,i}/w_{i}\Rightarrow\;\frac{ds_{i}}{d\eta_{i}} = \left( \mu\prime_{\eta,i}w_{i} - \mu_{\eta,i}^{2} \right)/w_{i}^{2} \right.$
  with
  $\mu\prime_{\eta,i} = \frac{d\mu_{\eta,i}}{d\eta_{i}} = \frac{d^{2}w}{d\eta_{i}^{2}} \equiv \mu''_{i}$
  (this is `d2mu.deta2(eta_i)` in code)
- $D_{i} = 1 + \lambda_{W}\left( w_{i} - W \right) + X_{\text{centered},i}^{T}\lambda_{x}$
  - $\frac{\partial D_{i}}{\partial\eta_{i}} = \lambda_{W}\mu_{\eta,i}$
  - $\frac{\partial D_{i}}{\partial z} = \frac{\partial\lambda_{W}}{\partial z} \cdot \left( w_{i} - W \right) - \lambda_{W} \cdot \frac{dW}{dz}$
  - $\frac{\partial D_{i}}{\partial\lambda_{x}} = X_{\text{centered},i}$

Define $\text{inv}_{i} = 1/D_{i}$ and the scalar term driving
$\beta$-equations:

$$T_{i} = s_{i} - \lambda_{W}\mu_{\eta,i}\text{inv}_{i},\quad s_{i} = \frac{\mu_{\eta,i}}{w_{i}}.$$

### Compute Its Derivatives

Using
$\,\mu\prime_{\eta,i} = d\mu_{\eta,i}/d\eta_{i} = {d2mu\, deta2}\left( \eta_{i} \right)$
and $\, dw_{i}/d\eta_{i} = \mu_{\eta,i}$,

$$\frac{\partial s_{i}}{\partial\eta_{i}} = \frac{\mu\prime_{\eta,i}w_{i} - \mu_{\eta,i}^{2}}{w_{i}^{2}}.$$

Also
$\,\frac{\partial\text{inv}_{i}}{\partial\eta_{i}} = - \text{inv}_{i}^{2} \cdot \frac{\partial D_{i}}{\partial\eta_{i}} = - \text{inv}_{i}^{2}\left( \lambda_{W}\mu_{\eta,i} \right)$.
Therefore

$$\frac{\partial T_{i}}{\partial\eta_{i}} = \frac{\mu\prime_{\eta,i}w_{i} - \mu_{\eta,i}^{2}}{w_{i}^{2}} - \lambda_{W}\mu\prime_{\eta,i}\text{inv}_{i} + \lambda_{W}^{2}\left( \mu_{\eta,i} \right)^{2}\text{inv}_{i}^{2}.$$

$$\frac{\partial T_{i}}{\partial z} = - \frac{\partial\lambda_{W}}{\partial z} \cdot \mu_{\eta,i}\text{inv}_{i} + \lambda_{W}\mu_{\eta,i}\text{inv}_{i}^{2} \cdot \frac{\partial D_{i}}{\partial z}$$

$$\frac{\partial T_{i}}{\partial\lambda_{x}} = \lambda_{W}\mu_{\eta,i}\text{inv}_{i}^{2} \cdot X_{\text{centered},i}$$

### Assemble Jacobian Blocks (with $a_{i}$ weights)

**$J_{\beta\beta}$ ($K \times K$)**:
$$J_{11} = \sum a_{i}Z_{i}^{T}\left\lbrack \frac{\partial T_{i}}{\partial\eta_{i}} \right\rbrack Z_{i}$$

**$J_{\beta z}$ ($K \times 1$)**:
$$J_{12} = \sum a_{i}Z_{i}^{T}\left\lbrack \frac{\partial T_{i}}{\partial z} \right\rbrack$$

**$J_{\beta\lambda}$ ($K \times L$)**:
$$J_{13} = \sum a_{i}Z_{i}^{T}\left\lbrack \frac{\partial T_{i}}{\partial\lambda_{x}} \right\rbrack$$

**$J_{z\beta}$ ($1 \times K$)**: derivative of W-equation w.r.t. $\beta$

Equation: $G_{W} = \sum a_{i}\left( w_{i} - W \right)\text{inv}_{i}$

$$\frac{\partial G_{W}}{\partial\eta_{i}} = a_{i}\left\lbrack \mu_{\eta,i}\text{inv}_{i} - \left( w_{i} - W \right)\text{inv}_{i}^{2}\left( \frac{\partial D_{i}}{\partial\eta_{i}} \right) \right\rbrack = a_{i}\left\lbrack \mu_{\eta,i}\text{inv}_{i} - \left( w_{i} - W \right)\text{inv}_{i}^{2}\left( \lambda_{W}\mu_{\eta,i} \right) \right\rbrack$$

Then: $J_{21} = \sum\frac{\partial G_{W}}{\partial\eta_{i}} \cdot Z_{i}$

**$J_{zz}$ ($1 \times 1$)**:
$$\frac{\partial G_{W}}{\partial z} = \sum a_{i}\left\lbrack - \frac{dW}{dz} \cdot \text{inv}_{i} - \left( w_{i} - W \right)\text{inv}_{i}^{2} \cdot \frac{\partial D_{i}}{\partial z} \right\rbrack$$

**$J_{z\lambda}$ ($1 \times L$)**:
$$\frac{\partial G_{W}}{\partial\lambda_{x}} = \sum a_{i}\left\lbrack - \left( w_{i} - W \right)\text{inv}_{i}^{2}X_{\text{centered},i} \right\rbrack$$

**$J_{\lambda\beta}$ ($L \times K$)**: constraints
$H(\lambda):\sum a_{i}\text{inv}_{i}X_{\text{centered},i} = 0$

$$\frac{\partial H}{\partial\eta_{i}} = - a_{i}\text{inv}_{i}^{2}\frac{\partial D_{i}}{\partial\eta_{i}}X_{\text{centered},i} = - a_{i}\text{inv}_{i}^{2}\left( \lambda_{W}\mu_{\eta,i} \right)X_{\text{centered},i}$$

Thus, component-wise
$J_{31} = \sum_{i}a_{i}\,\left( - \lambda_{W}\mu_{\eta,i}\,\text{inv}_{i}^{2} \right)\, X_{\text{centered},i}^{T}Z_{i}$.
In compact matrix form:

$$J_{31} = X_{\text{centered}}^{T}\operatorname{diag}\!( - a_{i}\,\lambda_{W}\,\mu_{\eta,i}\,\text{inv}_{i}^{2})Z.$$

**$J_{\lambda z}$ ($L \times 1$)**:
$$\frac{\partial H}{\partial z} = - \sum a_{i}\text{inv}_{i}^{2}\left( \frac{\partial D_{i}}{\partial z} \right)X_{\text{centered},i}$$

**$J_{\lambda\lambda}$ ($L \times L$)**:
$$\frac{\partial H}{\partial\lambda_{x}} = - X_{\text{centered}}^{T}\operatorname{diag}\left( a_{i}\,\text{inv}_{i}^{2} \right)X_{\text{centered}}.$$

These expressions match the unguarded analytic derivatives; in the code
(`src_dev/engines/el/impl/jacobian.R`), any terms involving derivatives
of $1/D_{i}$ are additionally multiplied by the active mask
$\mathbb{1}\{ D_{i}^{\text{raw}} > \delta\}$ to respect the denominator
floor used for numerical stability.

### Why Analytic A Helps

- Newton-Raphson (as used in our outer solve) linearizes $F(\theta)$
  near the current iterate:
  $F(\theta + \Delta) \approx F(\theta) + A(\theta)\,\Delta$. The update
  $\Delta$ solves $A\,\Delta = - F$, hence a high-quality $A$ is
  critical for fast, stable convergence.

## Variance Estimation

Analytic delta variance for the EL estimator has not yet been
implemented in the package. When requested via
`variance_method = "delta"`, the implementation returns `NA` with a
guidance message. We recommend `variance_method = "bootstrap"` for
standard errors in both IID and survey settings.

### Solving Strategy and Initialization

- Unknowns are $\theta = \left( \beta,z,\lambda_{x} \right)$ with
  $W = {plogis}(z)$. We solve the full stacked system $F(\theta) = 0$
  via Newton with the analytic Jacobian $A = \partial F/\partial\theta$
  using `nleqslv`.
- Globalization and scaling: we rely on `nleqslv`’s globalization
  (default `global = "qline"`, `xscalm = "auto"`) and enforce
  denominator positivity ($\min_{i}D_{i} \geq \varepsilon$) within
  equation evaluations. Optional standardization of design matrices
  improves conditioning.
- Initialization: by default $\beta$ starts at zeros in the scaled space
  (unless the user supplies `start$beta`), and $z$ is seeded at
  ${logit}\left( \text{observed response rate} \right)$. An internal
  last-chance Broyden retry may be used if Newton fails to converge;
  this is not a user-facing mode.

### Practical Identifiability and Diagnostics

The EL system balances the parametric response-model score against
calibration constraints. Identifiability can weaken in the following
situations:

- Weak or nearly collinear auxiliaries: if $X_{i} - \mu_{x}$ have little
  variation or are nearly collinear with the response score direction,
  the constraint block in $A = \partial F/\partial\theta$ becomes
  ill-conditioned.
- Inconsistent auxiliary means: if supplied $\mu_{x}$ are far from what
  the respondent sample can support (under the response model),
  denominators $D_{i}$ cluster near 0 and $\kappa(A)$ inflates.
- Heavy nonresponse or near-boundary $W$: when $W$ approaches 0 or 1,
  $\lambda_{W} = C/(1 - W)$ can spike and amplify sensitivity.

Diagnostics exposed by the implementation help assess these issues:

- `jacobian_condition_number` ($\kappa(A)$), `max_equation_residual`,
  denominator summaries (min, lower quantiles, median), weight
  concentration (max share, top-5 share, ESS), and the trimming
  fraction.

Mitigations include standardizing predictors, trimming extreme weights
(`trim_cap`), adding informative response-model predictors, and
preferring bootstrap variance when diagnostics indicate fragility.

## Survey Design Details

We extend QLS’s methodology to complex surveys in two complementary
ways:

- **Estimating equations with base weights:** All sums include the base
  weight $a_{i}$; set $a_{i}$ to the survey design weight for
  respondents. Totals $N_{\text{pop}}$ and
  $n_{\text{resp\_weighted}} = \sum a_{i}$ are computed from the design
  weights and used throughout the design-weighted system.

- **Nonrespondent total $T_{0}$ in the linkage equation:** In the
  survey-specific system we form
  $T_{0} = N_{\text{pop}} - n_{\text{resp\_weighted}}$ and enforce the
  linkage between $W$ and $\lambda_{W}$ through the equation
  $T_{0}/(1 - W) - \lambda_{W}\sum d_{i}/D_{i} = 0$ rather than using
  the closed-form
  $\lambda_{W} = \left( \left( N_{\text{pop}}/n_{\text{resp\_weighted}} \right) - 1 \right)/(1 - W)$.

- **Bootstrap variance via replicate weights:** For standard errors, we
  use bootstrap replicate-weight designs created with
  [`svrep::as_bootstrap_design`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html).
  For each replicate, the estimator is re-fit on a reconstructed design
  using that replicate’s weights, and
  [`survey::svrVar`](https://rdrr.io/pkg/survey/man/svrVar.html) is used
  to compute the variance of replicate estimates with appropriate
  scaling.

This matches the paper’s guidance to adapt the likelihood/estimating
framework to stratification or unequal-probability sampling while
relying on standard survey resampling for uncertainty. Analytic variance
has not been implemented yet.

Degrees-of-freedom: For confidence intervals, we use survey
degrees-of-freedom (t-quantiles) when a `survey.design` is supplied;
otherwise, we use normal quantiles.

## Scaling and Unscaling

### Scaling (optional; `standardize=TRUE`)

- **Compute a `nmar_scaling_recipe`**: for each column $j$ in $Z$ and
  $X$ (excluding intercept):
  - $\text{mean}_{j}$, $\text{sd}_{j}$; if $\text{sd}_{j} \approx 0$,
    set $\text{sd}_{j} = 1$ to avoid blow-ups.
- **Transform**:
  - $Z_{\text{scaled}}\lbrack,j\rbrack = \left( Z_{\text{un}}\lbrack,j\rbrack - \text{mean}_{j} \right)/\text{sd}_{j}$
  - $X_{\text{scaled}}\lbrack,j\rbrack = \left( X_{\text{un}}\lbrack,j\rbrack - \text{mean}_{j} \right)/\text{sd}_{j}$
  - $\mu_{x,\text{scaled}}\lbrack j\rbrack = \left( \mu_{x,\text{un}}\lbrack j\rbrack - \text{mean}_{j} \right)/\text{sd}_{j}$

### Unscaling $\beta$ and vcov

- **Construct linear map** $D$ of size $K \times K$:
  - For columns $j \neq$ intercept:
    $D\lbrack j,j\rbrack = 1/\text{sd}_{j}$
  - For intercept: adjust to absorb centering:
    $D\left\lbrack \text{intercept},j \right\rbrack = - \text{mean}_{j}/\text{sd}_{j}$
- **Transform**: $\beta_{\text{unscaled}} = D\beta_{\text{scaled}}$;
  $\text{vcov}_{\text{unscaled}} = D\,\text{vcov}_{\text{scaled}}\, D^{T}$

Code: centralized in `src_dev/shared/scaling.R`; engines call
[`validate_and_apply_nmar_scaling()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/validate_and_apply_nmar_scaling.md)
and
[`unscale_coefficients()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/unscale_coefficients.md).

## Bootstrap Variance

- **IID**:
  - Resample rows with replacement ($n$ to $n$), re-run estimator,
    compute $\text{var}$ of bootstrap $\widehat{Y}$s; warn if many
    failures; return $\sqrt{\text{var}}$.
- **Survey**:
  - Convert to bootstrap replicate-weight design via
    [`svrep::as_bootstrap_design`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.html).
  - For each replicate, re-construct a temporary design and run
    estimator; use
    [`survey::svrVar`](https://rdrr.io/pkg/survey/man/svrVar.html) to
    compute variance of replicate estimates (with scale/rscales).

Code mapping:

- Engine:
  `el_engine(..., family, standardize, trim_cap, variance_method, ...)`
  in `src_dev/engines/el/engine.R`
- Dispatch: `run_engine.nmar_engine_el(...)` in
  `src_dev/engines/el/run_engine.R` adapts the formula and forwards
  arguments to internal
  [`el()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el.md)
  methods.
  - [`el.data.frame()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_dataframe.md)
    /
    [`el.survey.design()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_survey.md)
    in `src_dev/engines/el/impl/dataframe.R` and
    `src_dev/engines/el/impl/survey.R` prepare inputs, call
    [`el_estimator_core()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_estimator_core.md),
    and wrap results.
- EL Core: `el_estimator_core(...)` in `src_dev/engines/el/impl/core.R`
  runs:
  - Construct $F(\theta)$ via
    [`el_build_equation_system()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_build_equation_system.md)
    (`src_dev/engines/el/impl/equations.R`).
  - Solve $F(\theta) = 0$ via `nleqslv` (Newton with analytic Jacobian
    when available, Broyden fallback).
  - Build EL weights, mean, and diagnostics.
- Jacobian: `el_build_jacobian(...)` in
  `src_dev/engines/el/impl/jacobian.R` returns analytic A whenever
  family supplies `d2mu.deta2` (logit, probit).
- Variance: Bootstrap variance is implemented in
  `src_dev/shared/bootstrap.R`.
- S3 result: `src_dev/engines/el/s3.R` provides print/tidy/glance and
  accessors.

### Practical Notes

- Denominator guard: $D_{i} \geq \varepsilon$ (default $10^{- 8}$)
  across all steps; diagnostics report extreme fractions.
- Eta cap option: you can adjust the $\eta$ cap via
  `options(nmar.eta_cap = 60)` (default is 50) to suit your data scale
  and link

### Algorithm

We solve the full stacked system $F(\theta) = 0$ with Newton using the
analytic Jacobian $A = \partial F/\partial\theta$ and globalization via
`nleqslv`. Denominator positivity ($\min_{i}D_{i} \geq \varepsilon$),
predictor standardization, and capped $\eta$ ensure numerical stability.
The estimating equations remain those of Qin, Leung and Shao (2002) for
the IID path, and a design-weighted analogue for the survey path.

``` text
Input: Z (response design), X (auxiliary design), mu_x (population means),
       a (base weights), family (logit/probit), trim_cap, tolerances.
Initialize: beta = 0 in scaled space (or user-supplied start),
            z = logit(observed response rate), lambda_x = 0
            (and lambda_W = 0 for survey designs).
Repeat until convergence of F(theta) = 0:
  1) Compute eta = Z beta, w = linkinv(eta), W = plogis(z).
     - IID (data.frame): set lambda_W = ((N_pop/n_resp_weighted) - 1)/(1 - W).
     - survey.design: use the current lambda_W component of theta.
  2) Evaluate full stacked equations using guarded denominators
     D_i = 1 + lambda_W (w_i - W) + (X_i - mu_x)^T lambda_x.
  3) Compute analytic Jacobian A = dF/dtheta (if available; else numeric/Broyden).
  4) Newton step: solve A * step = -F with globalization; enforce min D_i >= eps.
  5) Update theta <- theta + step.
Return: p_i \propto a_i / D_i and \hat{Y} = Sum p_i Y_i / Sum p_i.
```

## References

- Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data
  under nonignorable nonresponse or informative sampling. Journal of the
  American Statistical Association, 97(457), 193-200.
- Chen, J., and Sitter, R. R. (1999). A pseudo empirical likelihood
  approach for complex survey data. Biometrika, 86(2), 373-385.
- Wu, C. (2005). Algorithms and R codes for the pseudo empirical
  likelihood method in survey sampling. Canadian Journal of Statistics,
  33(3), 497-509.

## Appendix: EL Engine API Reference (User-Facing)

This appendix summarizes the key options of the EL engine (constructor:
[`el_engine()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_engine.md)),
their defaults, and recommended usage.

- **family** (default: “logit”)
  - Values: “logit”, “probit”, or a family object (list with `name`,
    `linkinv`, `mu.eta`, `d2mu.deta2`, `score_eta`).
  - Notes: We implement
    [`logit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/logit_family.md)
    and
    [`probit_family()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/probit_family.md).
    Both use the log-likelihood score
    `score_eta(eta, delta) = mu.eta(eta)/linkinv(eta)` (for
    respondents), i.e., $\partial\log p/\partial\eta$. This matches the
    paper’s semiparametric MLE equations and keeps the analytic Jacobian
    family-agnostic.
- **standardize** (default: TRUE)
  - Standardize $Z$/$X$ (and $\mu_{x}$) using a `nmar_scaling_recipe`
    for numerical stability. Coefficients and vcov are unscaled after
    solving.
- **trim_cap** (default: Inf)
  - Caps EL weights and redistributes mass. Improves robustness when
    extreme weights arise. Prefer `variance_method = "bootstrap"` when
    trimming is finite.
- **variance_method** (default: “none”)
  - “bootstrap”: IID resampling or survey replicate weights via `svrep`;
    preferred for SEs.
  - “none”: skip variance calculation.
- bootstrap_reps (default: 500)
  - Number of bootstrap replicates for `variance_method = "bootstrap"`.
    Increase for stability, decrease for speed.
- **control** (default: list())
  - Passed to `nleqslv` (e.g., `ftol`, `xtol`, `maxit`).

### Diagnostics (glance/print)

- `jacobian_condition_number` ($\kappa(A)$), `max_equation_residual`,
  denominator summaries (min, 1%/5%/median), weight concentration (max
  share/top-5/ESS), and trimming fraction. These help assess
  identification and numerical stability.
