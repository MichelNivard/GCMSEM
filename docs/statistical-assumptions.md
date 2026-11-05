# Statistical model, uncertainty, and interpretation

GCMSEM separates two estimators: an LD-aware estimator of a distribution of
joint SNP effects, and a structural equation model for that distribution. Good
software diagnostics do not establish the scientific assumptions of either
estimator. The validation scripts report the scenarios actually tested; they
do not certify coverage for every GWAS architecture.

## 1. The quantity being estimated

For genotype-standardized effects, let

```math
\widehat b_{ja}=\sum_k R_{jk}\beta_{ka}+\epsilon_{ja}.
```

Here `j` is a regression/tag SNP, `k` is a reference/source SNP, and `a` is a
trait. `beta` is the joint effect conditional on other variants; a GWAS reports
a noisy marginal association. GCMSEM estimates moments of the source effects,
including effects equal to zero. These moments are **not moments of people or
phenotypes**, and an SNP-effect variance is not automatically a heritability.

For a per-allele coefficient, the conversion is
`beta_sd = beta_allele * sqrt(2 * AF * (1-AF))`, with the same conversion for its
standard error. This uses reference-panel allele frequencies and an HWE
genotype-variance convention. The phenotype scale is retained: inverse-normal,
raw physical units, log odds, and linear-probability effects are not
interchangeable. There is no automatic liability-scale conversion. Ancestry,
build, alleles, phenotype units, variant inclusion, and sample ascertainment
must be compatible across the inputs.

The source universe and the tag universe are distinct. Restricting regression
rows to HM3 does not imply that only HM3 variants contribute to the association.
Source counts and mixture weights remain those of the declared reference
universe when tag SNPs are filtered or a jackknife group is deleted. Missing
source variants, long-range LD, a mismatched reference population, or a too
short LD window can change the target or bias the fit.

## 2. LD-power regressions

Source variants belong to disjoint annotation strata `c`. Within a stratum,
their effect vectors have common cumulants through order four. Effect vectors
at different source SNPs are independent under the working model. Thus

```math
L^{(q)}_{jc}=\sum_{k\in c} R_{jk}^{q},\qquad
\mathrm{cum}(b_{ja_1},\ldots,b_{ja_q})
=\sum_c L^{(q)}_{jc}\kappa_{c,a_1\ldots a_q}.
```

`L1` and `L3` retain LD signs. Replacing a cube by an absolute cube changes the
estimand. Odd-order moments also depend on the chosen effect-allele convention;
they are only interpretable with consistent signed LD and association alleles.

At each order a regression estimates the stratum coefficients. The design
must have sufficient independent variation: collinear annotation scores, too
few tag rows, or an unidentifiable noise intercept are genuine estimation
failures. Numerical regularization must not quietly turn such failures into
apparently identified estimates. Reported condition and rank diagnostics
should be reviewed.

The random-effects assumptions can fail when causal effects cluster within LD
blocks, depend on unmodeled annotations, or are correlated across nearby
variants. The block jackknife addresses dependence in the uncertainty
calculation; it does not repair a misspecified mean regression.

## 3. Sampling noise and sample overlap

The default error model is joint mean-zero Gaussian sampling error. At a SNP,

```math
\Omega_{j,ab}=s_{ja}s_{jb}I_{ab},
```

where `s` denotes the reported standard error on the analysis scale. After
fitting marginal means, a second-order regression includes both the LD scores
and `s_a*s_b` as predictors:

```math
E[(\widehat b_{ja}-\mu_{ja})(\widehat b_{jb}-\mu_{jb})]
=\sum_c L^{(2)}_{jc}\kappa_{c,ab}+s_{ja}s_{jb}I_{ab}.
```

This is an LDSC-style intercept model expressed on the effect scale. The
off-diagonal noise coefficients represent correlated estimation error,
including sample overlap; they are not genetic correlations or direct
estimates of the number of overlapping participants. Shared confounding can
also contribute. This model requires the relevant error contribution to have
the specified SE-product form and to be distinguishable from LD-dependent
signal. An intercept does not remove arbitrary population stratification,
selection bias, or nonlinear/binary-outcome approximation error.

Gaussian product corrections use the **whole** fitted noise matrix. For
centered estimates `x`, for example,

```math
H_{abc}=x_ax_bx_c-\Omega_{ab}x_c-\Omega_{ac}x_b-\Omega_{bc}x_a,
```

and the fourth-order correction subtracts the six covariance-times-product
terms and adds the three covariance-pairing terms. Repeated trait indices are
occurrences, so their multiplicities matter. The fourth-order LD regression
also removes the three products of fitted signal covariances to recover an
**excess cumulant**, rather than a raw fourth moment. No sample splitting is
needed.

Gaussian sampling errors do not require Gaussian true SNP effects. They do
exclude unknown nonzero third/fourth error cumulants. The identities are exact
when the required nuisance quantities are known; using estimated means,
covariances, and intercepts introduces finite-sample plug-in error. Every stage
must therefore be refitted inside every genomic deletion. A supplied fixed
noise matrix is treated as known; uncertainty from estimating it elsewhere is
not automatically propagated.

## 4. Pooling strata

Cumulants do not generally average across a mixture. GCMSEM first converts
stratum cumulants to raw moments, averages those raw moments using fixed source
variant proportions, and then converts back to cumulants. This retains
between-stratum differences in means and variances. A simple weighted average
of fourth cumulants would, for example, miss the excess kurtosis caused by a
mixture of Gaussian strata with different variances.

Even if SEM disturbances are independent **within** each MAF stratum, mixing
strata can make them dependent. A Gaussian disturbance within each stratum can
also become heavy-tailed after mixing. Consequently an independent-source SEM
for pooled moments is an additional assumption. Inspect conditional results
or fit stratum-specific models when this distinction is consequential.

## 5. Reference uncertainty and genomic deletions

Reference genotypes are finite samples. The conventional second-order
LD-score correction reduces finite-reference bias in squared correlations; it
is not a general unbiased estimator for all reference-panel designs. Third
and fourth LD powers are plug-in estimates, and their reference-sampling bias
is not removed. Unless an outer reference resampling analysis is supplied,
all reported uncertainty conditions on the fitted LD scores, annotation
assignment, and reference variant counts.

Genomic deletion groups should contain whole dependent LD blocks. Approximately
equal tag-SNP counts are useful, but no grouping should split an LD block just
to achieve equal lengths. Few groups, extremely uneven groups, strong
dependence between groups, or a dominant genomic region limit jackknife
reliability. Tag-row count is a grouping and weighting convention, not an
estimate of the number of independent genetic observations.

For unequal deletion fractions `f_b`, the weighted delete-group construction
uses

```math
\tau_b=\widehat m/f_b-(1/f_b-1)\widehat m_{-b},\qquad
\widehat m_J=\sum_b f_b\tau_b,
```

```math
\widehat V=\frac1B\sum_b
\frac{(\tau_b-\widehat m_J)(\tau_b-\widehat m_J)^T}{1/f_b-1}.
```

For equal groups this becomes the usual delete-block jackknife covariance.
The covariance formula is formally associated with the jackknife-corrected
estimate; reporting the full-data point estimate with that uncertainty is a
first-order convention. It is not an exact guarantee for heterogeneous LD
blocks or nonlinear fourth-order estimates. See
[Patterson's weighted-jackknife note](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf)
and Busing, Meijer and van der Leeden (1999), *Delete-m jackknife for unequal m*,
*Statistics and Computing* 9, 3–8.

The number of unique cumulants through order four is
`sum(choose(p+q-1,q), q=1,...,4)`. The full sampling covariance therefore grows
quickly, even when moment estimation is streamed. With `B` deletion groups,
its rank is at most `B-1`. Selecting a scientifically relevant subset of
moments can be necessary. A large SNP count does not remove this rank bound.

## 6. SEM meaning and inference

The R package propagates disturbance cumulants through a linear structural
model. For structural coefficient matrix `B`,

```math
z=(I-B)^{-1}e=T e,\qquad
K^{(q)}_z=T^{\otimes q}K^{(q)}_e.
```

Observed and latent node distributions are disturbance distributions. Marking
a disturbance Gaussian fixes its third cumulants and excess fourth cumulants
to zero. It does **not** force its observed outcome to be Gaussian after
receiving non-Gaussian upstream influences. Conversely, omitting third/fourth
observed moments from the fitting objective does not impose Gaussianity.
Fourth cumulant zero means ordinary univariate kurtosis three after
standardizing a nondegenerate variable.

DWLS chooses parameters using inverse marginal sampling variances. Robust
parameter uncertainty retains the full sampling covariance:

```math
W=\mathrm{diag}(1/\widehat V_{ii}),\quad
A=J^TWJ,\quad
\mathrm{Var}(\widehat\theta)
=A^{-1}J^TW\widehat VWJA^{-1}.
```

There is no additional division by GWAS sample size. A singular sampling
covariance does not itself prevent this sandwich computation because DWLS
does not invert it. These SEs propagate correlated moment uncertainty under
the regular, correctly specified SEM approximation. They are not generally
robust to an arbitrarily misspecified SEM; in particular, the usual
`J' W J` bread is not a general misspecification-robust Hessian. A
rank-deficient parameter Jacobian or non-convergence can prevent a usable
standard error. A local Gaussian quadratic-form approximation supplies an SE for Q and a
moment-matched scaled-chi-square diagnostic reference. It accounts for fitting
parameters by projecting off the model Jacobian, while treating V and that
tangent space as fixed. Its p-value is exploratory; calibrated higher-order
model-fit inference and a WLSMV correction are not established. See the
[derivation and limitations](../tutorials/R-results.md#how-the-approximate-se-is-calculated).

Exact structured derivatives remove finite-difference error; multiple starts
help detect numerical minima. Neither establishes identification. A zero
factor loading can leave its higher moments unidentified. Equal loadings,
fixed paths, factor scales, independent sources, and reciprocal stability are
substantive restrictions. Local derivative rank does not establish a unique
global model, a causal direction, or validity of a selected model's standard
errors after data-dependent selection. Reciprocal source relabeling can give
different paths with the same observed moments.

The default, `feasibility = "check"`, fits the **unpenalized DWLS objective**.
It does not constrain source variances to be nonnegative and does not impose
positive-semidefinite source moment matrices. Numerically converged candidate
solutions are compared by DWLS discrepancy without preferring an admissible
source distribution. A negative source variance, violated Pearson inequality,
or indefinite fourth-order source moment matrix is retained in the results
and prominently flagged. Such classical or higher-order Heywood cases are
not silently dropped, clipped, or replaced by a worse-fitting admissible
solution. Across simulation replicates or studies, filtering estimates by
admissibility can itself alter the distribution of reported estimates.

There are two distinct diagnostic questions. First, is a solution numerically
regular enough to calculate the local sandwich? Second, do its fitted source
moments satisfy necessary conditions for a probability distribution? A source
Heywood case can have a full-rank Jacobian and a well-defined **formal**
sandwich SE. That SE is retained, while `source_admissible` and
`inference_valid` make the problem explicit. Retaining a number does not make
it a validated probability-model confidence interval. Numerical failure,
non-convergence, rank deficiency, or an imposed active constraint can still
prevent regular SE calculation; those failures must not be hidden by
fabricating an SE. Point estimates and their diagnostics remain available.

An explicitly requested `feasibility = "penalty"` changes the estimation
procedure by adding a finite quadratic penalty for negative eigenvalues of
source moment matrices. If `M` is the source raw-moment matrix indexed by
monomials of degree zero, one, and two, this optional penalty is

```math
P(\theta)=\sum_i\{\min(\lambda_i(M(\theta)),0)\}^2,
\qquad Q_{\mathrm{pen}}=Q_{\mathrm{DWLS}}+\rho P.
```

Source units are normalized using fixed scales determined at initialization;
this scaling does not change positive semidefiniteness. The code differentiates
the matrix entries and squared distance to the PSD cone exactly. When this
optional procedure is selected, the penalty weight starts at `1e6` by default
and increases by factors of `1e3`, up to `1e12`, when independent, normalized
feasibility checks fail. With more than 12 source nodes, the penalty and
diagnostics use the source covariance and separate univariate moment
matrices; they do not test the full joint fourth-order matrix. This limitation
is reported in the source diagnostics.

This is **approximate numerical feasibility enforcement**, not an exact
semidefinite-constrained solver. A finite penalty can leave a violation, and a
large penalty can make optimization ill-conditioned. An optimizer success
code alone is insufficient: gradient, rank, and imposed-boundary checks
remain necessary. The reported `objective` remains the unpenalized DWLS
discrepancy; the penalized objective, penalty weight, and residual penalty are
stored separately. Penalized and unconstrained runs should be reported as
different estimators, including any change in bias or admissibility frequency.
The optional `stable` setting filters reciprocal solutions; it is not a
stability-constrained optimizer and need not find every stable solution.

Positive semidefinite covariance and fourth-order moment matrices, together
with univariate Pearson inequalities, are necessary conditions. They are not
a proof that an arbitrary multivariate truncated moment sequence has a
probability distribution. An imposed active moment boundary requires
nonregular inference even when it satisfies the numerical feasibility checks.
In the unconstrained fit, source-boundary proximity is a distributional warning
and is reported separately from numerical identification of the algebraic
moment model.

## 7. Appropriate claims from a run

A successful run establishes that the stated inputs were processed under the
recorded settings and that the reported numerical checks passed. Stronger
claims require evidence: ancestry and allele matching, suitable LD coverage,
an adequate effect architecture, sufficiently independent deletion groups,
and simulation or external validation of the intended scientific use. In
particular, a lower moment discrepancy alone is not proof of causation.
