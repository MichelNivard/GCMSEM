# Work with GCMSEM results in R

Run the [small worked example](Quickstart.md) first. It creates `moments`, `scaled`, `model` and `fit` in your R session. These examples use those objects.

`estimate_moments()` and `gcmsem()` return R objects and create no result files by default. Typing an object prints a compact overview. `summary(fit)` displays its diagnostics and parameter table; assigning `s <- summary(fit)` keeps a structured summary without printing. Use `as.data.frame()` to extract the moment or parameter table.

| R expression | Contents |
|---|---|
| `moments$moments` | Data frame of cumulant keys, orders, estimates and jackknife SEs |
| `moments$S`, `moments$S_SE` | SNP-effect covariance and SE matrices for all traits |
| `moments$R`, `moments$R_SE` | Genetic correlations and delta-method SEs using the full sampling covariance |
| `moments$V` | Full joint sampling covariance of the cumulant estimates |
| `moments$jackknife` | Matrix of genomic deletion estimates |
| `moments$metadata`, `moments$log` | Provenance, QC, settings and the Python estimation log |
| `parameter_estimates(fit)` | Numeric parameter estimates, formal SEs and inference-status labels |
| `coef(fit)`, `vcov(fit)` | Named coefficient vector and parameter sampling covariance |
| `fit$diagnostics`, `fit$log` | Convergence, identification, source-moment checks and fitting log |

`summary(moments)` shows the covariance and correlation matrices with SEs. Console previews show at most six traits; `print(moments, max.traits = Inf)` shows all traits. Every matrix in the object is complete.

Save the entire object with ordinary R serialization, or request tabular files explicitly:

```r
saveRDS(fit, "toy/fit.rds")
restored <- readRDS("toy/fit.rds")

export_gcmsem(moments, "toy/exported-moments")
export_gcmsem(fit, "toy/exported-fit")
reloaded <- read_moments("toy/exported-moments")

# This option still returns the fitted R object.
fit_with_files <- gcmsem(
  model, scaled, n_starts = 8, seed = 20261103,
  write.to.tsv = TRUE, out = "toy/another-fit"
)
```

Exports refuse existing files unless `overwrite = TRUE`. `estimate_moments()` accepts the same `write.to.tsv` and `out` arguments. `log_file` explicitly writes a fit log, while `verbose = TRUE` displays progress in the console. Logs remain available in the returned objects.

The Python CLI remains available for terminal workflows and writes versioned moment bundles; `read_moments(prefix)` loads them into the same R object format. [The estimator guide](../docs/moment-estimator.md) describes each bundle file. Multiple MAF strata produce conditional objects in `moments$conditional`. Pooling uses raw-moment mixture identities. Source counts remain fixed when tag rows are filtered or deleted. Jackknife groups preserve whole LD blocks and account for unequal deletion sizes.

## Assess model fit before interpreting paths

A fitted SEM predicts each moment you supplied. The **residual** is the estimated moment minus the prediction. GCMSEM's DWLS discrepancy, `Q`, adds the squared residuals after dividing each by its own sampling variance. A smaller Q means the model reproduces these moments more closely. Even a correct model has nonzero Q because the input moments are estimated with noise.

```r
fit_measures(fit)
residual_table <- data.frame(
  moment = names(fit$observed),
  observed = unname(fit$observed),
  predicted = unname(fitted(fit)),
  residual_in_SEs = unname(fit$residuals / sqrt(diag(fit$V)))
)
residual_table[order(-abs(residual_table$residual_in_SEs)), ]
```

`fit_measures()` returns an ordinary one-row data frame:

| Column | Interpretation |
|---|---|
| `Q` | The minimized DWLS discrepancy. |
| `Q_se` | A rough sampling SD for Q, using the observed residuals as an estimate of population misfit and including quadratic sampling noise. |
| `Q_se_delta` | The first-order contribution alone; it approaches zero when the model fits the moments exactly. |
| `Q_null_mean`, `Q_null_sd` | Approximate mean and SD of Q if the specified moment model is correct. They account for correlations among moments and for parameters estimated from those moments. |
| `Q_null_95` | Approximate 95th percentile of that null reference. |
| `Q_p_approx` | An exploratory upper-tail probability from a scaled chi-square reference matching the null mean and variance. |
| `Q_reference_scale`, `Q_reference_df` | The scale and effective degrees of freedom used by the approximation. |
| `residual_df` | Number of fitted moments minus free parameters; it need not equal the effective reference degrees of freedom. |
| `Q_status` | Whether the calculation is a local approximation, formally flagged, or unavailable. |

Start with the convergence and Jacobian-rank diagnostics. Then inspect Q relative to its approximate null distribution and look at the residual table. A large Q and several large residuals indicate that the specified model fails to explain the estimated moments. Do not interpret the paths of a badly fitting model as supported effects. Keep that fit and its diagnostics available so readers can see why the model was set aside.

Compare Q across candidate models only when they use the same input moments, scaling and weights. Fitting more parameters can reduce Q, so read the parameter counts and model restrictions too. The forward and reverse models in the [causal tutorial](Causal-examples.md) are not nested; subtracting their Q values does not give an ordinary chi-square difference test. Even for nested models, this DWLS implementation needs a dedicated difference-test calculation. The two Q statistics also share sampling error, so their separate SEs cannot be combined as if independent.

### How the approximate SE is calculated

Let $m$ contain estimated moments, $f(\widehat\theta)$ their fitted values and $D=\mathrm{diag}(1/\sqrt{V_{ii}})$. Define standardized residuals $e=D[m-f(\widehat\theta)]$, moment-error correlation $C=DVD$, and the standardized model Jacobian $J=D(\partial f/\partial\theta)$. The projection $P=I-JJ^+$ removes changes that can be absorbed by refitting the parameters. Here $J^+$ is the Moore–Penrose inverse; the implementation uses an SVD.

```math
C_r=PCP,\qquad Q=e^\top e,
```

```math
E_0(Q)\approx\mathrm{tr}(C_r),\qquad
\mathrm{Var}_0(Q)\approx2\mathrm{tr}(C_r^2),
```

```math
\widehat{\mathrm{SE}}(Q)
\approx\sqrt{2\mathrm{tr}(C_r^2)+4e^\top C_r e}.
```

These are Gaussian quadratic-form identities applied to a locally linearized fitted model. The last equation substitutes observed residuals for unknown population misfit. Because those residuals also contain noise, it can overestimate the sampling SD near a well-fitting model. `Q_null_sd` describes the correctly specified case directly. Neither quantity justifies a symmetric normal confidence interval around a nonnegative Q statistic.

The reference matches the first two null moments to $a\chi^2_\nu$, with $a=\mathrm{tr}(C_r^2)/\mathrm{tr}(C_r)$ and $\nu=\mathrm{tr}(C_r)^2/\mathrm{tr}(C_r^2)$. When standardized moment errors are independent, it reduces to the familiar chi-square reference with moments-minus-parameters degrees of freedom. Correlated moments generally produce a weighted sum of chi-squares, which this reference only approximates. See [Yuan and Bentler's quadratic-form approximations](https://doi.org/10.1348/000711009X449771).

The calculation holds `V`, the weights, and the fitted tangent space fixed. Nonlinear curvature, uncertainty in `V`, LD-reference uncertainty, source-boundary behavior and selection among models can affect its accuracy. **The package has not established calibrated Q p-values for higher-order GWAS moments.** Use this as a diagnostic reference, inspect residuals and validate consequential decisions by simulation. The existing parameter SE and higher-order coverage limitations remain relevant.

Nonconvergence, deficient Jacobian rank and penalized fitting produce `NA` Q uncertainty with a reason. Algebraically regular unconstrained Heywood fits retain formal Q diagnostics with a flag, just as they retain formal parameter SEs. Passing a fit check does not repair impossible source moments or establish a causal interpretation.
