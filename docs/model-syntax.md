# GCMSEM model syntax

GCMSEM models the second, third and fourth **cumulants of joint SNP effects**. A marginal fourth cumulant is excess kurtosis times variance squared, not the fourth central moment. Model equations refer to genomic effects, not directly to an individual-level phenotype regression.

## Define a model

```r
library(GCMSEM)
model <- gcmsem_model('
  f =~ a*x + a*y
  f ~~ 1*f
  y ~ b*x
', observed = c('x', 'y'),
   gaussian_observed = FALSE, gaussian_latent = FALSE)
```

This is `x = a f + u_x`, `y = b x + a f + u_y`. The shared label `a` imposes equal **direct** loadings in the model's input units. The total loading of `f` on `y` also includes the path through `x`. Only the `x -> y` path is present. The latent-source variance is fixed to one and all three source marginal third/fourth cumulants are free. The model has ten free parameters for twelve observed cumulants. Parameter counting is not an identification test; inspect the fitted Jacobian rank.

Defaults are explicit:

- Every observed node has a free residual variance. Every latent node has variance fixed at one, providing a scale; override this with an explicit `~~` term if you provide another scale constraint.
- Source covariances and mixed higher source cumulants default to zero.
- `gaussian_observed = TRUE` and `gaussian_latent = TRUE` fix marginal source third/fourth cumulants to zero. Set either to `FALSE` to free those cumulants for that class of nodes, or pass a character vector naming only the nodes to keep Gaussian.
- Explicit higher-cumulant syntax overrides the corresponding Gaussian default. Thus the flags specify defaults, not an immutable Gaussian constraint.
- Observed higher moments can be nonzero even when that node's own disturbance is Gaussian, because non-Gaussian upstream sources propagate through paths.
- There is no automatic equality of loadings, residual covariances, causal direction selection, or zero-confounder model selection.

## Paths, factors and covariances

| Syntax | Meaning |
|---|---|
| `y ~ b*x` | Path from x to y, labelled b |
| `f =~ x + y + z` | Latent factor with three freely estimated loadings |
| `x ~~ v*x` | Source/residual variance, labelled v |
| `x ~~ c*y` | Source/residual covariance |

In these statements, numbers before `*` fix a coefficient and labels before `*` name free parameters. Reusing a label ties coefficients, including across statements.

## Higher source cumulants: explicit tuples

Write each variable occurrence inside `k3(...)` or `k4(...)`; the tuple length specifies the cumulant order. Put a parameter label inside `free("label")` on the right of `=`. This separates variable names from parameter labels, even when they share the same spelling.

| Syntax | Meaning |
|---|---|
| `k3(x, x, y)` | Free source cumulant $\kappa_3(x,x,y)$ |
| `k3(x, x, y) = free()` | Same anonymous free parameter |
| `k3(x, x, y) = free("s")` | Same source cumulant, labelled s |
| `k3(x, y, y) = free("t")` | Source cumulant $\kappa_3(x,y,y)$, labelled t |
| `k4(x, x, y, z) = free("k")` | Source cumulant $\kappa_4(x,x,y,z)$, labelled k |
| `k4(x, x, x, x) = 0` | Fix x's marginal fourth cumulant to zero |
| `k3(f, f, f) = free("s_f")` | Free f's marginal third cumulant |
| `k4(f, f, f, f) = free("k_f")` | Free f's marginal fourth cumulant |

The repeated names are the indices: `k3(x, x, y)` contains two copies of x, while `k3(x, y, y)` contains two copies of y. A cumulant is symmetric, so `k3(y, x, x)` denotes the same quantity as `k3(x, x, y)`. Tuple ordering does not specify a causal direction.

The right side is optional. If present, it must be `free()`, `free("label")`, or a finite numeric constant, such as `0`, `-0.02` or `1e-6`. Single-quoted labels are also accepted. Labels and variable names must begin with a letter and contain letters, digits, underscores or periods. Reusing a label imposes equality on the original parameter scale, including across cumulant statements and path or covariance statements. For example, `k3(x, x, y) = free("x")` unambiguously uses x as both a variable name and a parameter label. Choose equality constraints with the parameters' units in mind.

These are declarations **inside the model string**, not R function calls. The parser does not evaluate R code, expressions or function calls in the model. Use an outer single-quoted R string when the labels use double quotes:

```r
model <- gcmsem_model('
  f =~ a*x + a*y
  y ~ b*x
  k3(f, f, f) = free("skew_f")
  k4(f, f, f, f) = free("fourth_f")
  k3(x, x, y) = free("xxy")
  k3(x, y, y) = free("xyy")
', observed = c('x', 'y'))
```

Here the Gaussian defaults still apply to undeclared marginal source cumulants. The four explicit declarations override their respective defaults; freeing one third cumulant does not automatically free a fourth cumulant. Mixed higher source cumulants remain zero unless declared. To fit only selected **observed** orders, use the separate `orders` argument to `gcmsem()` below.

The parser accepts `#` comments and newline or semicolon statement separators. `+` separates terms in path, factor and covariance statements; a higher-cumulant statement declares one tuple. GCMSEM deliberately does not implement lavaan's thresholds, nonlinear constraints, defined parameters, intercept models, groups or full parser grammar. Unsupported syntax errors explicitly.

### Migrate earlier higher-moment syntax

The directional-looking aliases `~--` and `~~-` are rejected with a replacement in the error message. Their historical meaning is preserved in that replacement; it is never silently reversed. The older `~~~` and `~~~~` operators remain temporarily readable for existing saved models and emit a deprecation warning with the explicit tuple expansion. Use tuples for new or revised models.

| Earlier statement | Replacement |
|---|---|
| `x ~-- y` | `k3(x, x, y)` |
| `x ~~- y` | `k3(x, y, y)` |
| `x ~~~ s*x*y` | `k3(x, x, y) = free("s")` |
| `x ~~~~ k*x*y*y` | `k4(x, x, y, y) = free("k")` |
| `f ~~~~ 0*f*f*f` | `k4(f, f, f, f) = 0` |

## Variables and model scope

All variables must be in `observed`, declared by `=~`, or passed in the `latent` argument. Latent regressions and reciprocal paths are allowed. A nonsingular `I-B` is required. `stable=TRUE` at fitting selects only available converged solutions whose feedback matrix has spectral radius less than one; this is an extra assumption and does not prove causal validity.

A multi-indicator factor needs no special two-trait mode:

```r
cfa <- gcmsem_model('
  F =~ x1 + x2 + x3 + x4
  x1 ~~ r*x2
', observed = c('x1','x2','x3','x4'), gaussian_latent = FALSE)
```

## Read and fit

```r
bundle <- read_moments('results/traits')
fit <- gcmsem(model, bundle, orders = 2:4, n_starts = 10,
              seed = 42, log_file = 'results/model.log')
summary(fit)
parameter_estimates(fit)
vcov(fit)
fit$starts
```

`orders=2` fits covariance information only; it does **not** change a source distribution or remove higher parameters already declared in the model. Higher parameters that then have no information cause a rank warning and withheld standard errors. Alternatively supply `keys=c('x|x','x|y','y|y','x|x|y')`. Full sampling covariance is subset using keys. The count of symmetric observed cumulants is `choose(N+1,2) + choose(N+2,3) + choose(N+3,4)`; no restriction to two traits is built in.

`gcmsem` uses DWLS for the estimates and the full supplied sampling covariance for sandwich uncertainty. The latter is already the covariance of the estimated moments: there is no additional division by GWAS sample size. A singular sampling covariance is accepted if it is positive semidefinite and every used diagonal is positive. Its rank and the sandwich rank are reported.

A compiled C++ tensor contraction and an independent R reference implementation are included; `implied_moments(..., backend="R")` selects the reference. The optimizer uses exact structured first derivatives of `(I-B)^-1` and tensor products. Parameter scales and starts are derived from observed SNP-effect variances, so small genomic units do not require artificial variance floors. Shared-label equalities remain on the original input scale. Estimates and standard errors are returned in that original model scale. User-specified fixed values are never rescaled to a different scientific model.

## Diagnostics and limits

The default fitter runs deterministic multistart `nlminb`, reports each start and selects the lowest raw DWLS discrepancy among converged solutions. Source admissibility never filters or replaces the default estimate. Negative free residual variances are allowed, as in unconstrained SEM; fixed values remain fixed. A requested `stable=TRUE` branch restriction is an explicit additional user assumption. It checks nonnegative source variances, source covariance PSD, marginal Pearson bounds and (up to twelve non-deterministic source nodes) PSD of the joint raw moment matrix for monomials of degree at most two. Known deterministic zero sources are removed after verifying that all incident cumulants are zero. These are **necessary**, not sufficient, conditions for a joint distribution.

The default `feasibility="check"` diagnoses source moment conditions without imposing them. A negative source variance emits a **CLASSICAL HEYWOOD CASE** warning. Impossible higher-order source moments emit a **HIGHER-ORDER HEYWOOD CASE** warning naming the failed Pearson bound or joint moment-matrix check. Such estimates can arise through sampling variation; they are retained, not clipped, dropped or replaced by a constrained solution.

All point estimates remain available. When convergence and the Jacobian permit a sandwich calculation, formal standard errors also remain available in `se`. `se_status="formal_only_flagged"` distinguishes these from `regular_local` results. The separate flags `source_admissible`, `formal_se_available` and `inference_valid` prevent mathematical computability from being confused with regular interpretable inference. Conditional Wald p-values are left unavailable for flagged fits. Only numerical failures (nonconvergence, large gradient, rank deficiency or invalid sandwich variance) make formal SEs unavailable; the reason is stated.

For an explicitly requested sensitivity analysis, `feasibility="penalty"` adds a squared negative-eigenvalue penalty with exact derivatives and increasing strength until independently normalized necessary source checks pass, or the maximum penalty is reached. This is a **different estimator**, labelled `DWLS with source-feasibility penalty`. The raw DWLS `objective` is separate from `penalized_objective` and `constraint_penalty`. Penalty enforcement does not guarantee exact feasibility or a global optimum, and it is never applied silently. For more than twelve source nodes the full joint moment-matrix check/penalty is skipped explicitly; covariance and marginal checks remain.

Boundaries and inadmissibility flag ordinary local inference without suppressing available formal estimates. Local rank does not prove global identification: reciprocal and overcomplete source assignments can produce different parameters with identical observed moments.

`Q` is a descriptive DWLS discrepancy, not a calibrated chi-square statistic. The package does not report WLSMV corrections or causal-direction p-values. Parameter p-values, when provided, are conditional Wald calculations for a regular fixed model and solution. Model selection, reference-LD uncertainty, and misspecification can require additional calibration beyond these SEs.
