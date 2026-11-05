# Specify a structural equation model in R

A structural equation model describes which variables influence one another, which share a latent source, and which source distributions supply the variation. Begin with the [small worked example](Quickstart.md) to obtain a moment object. Then define the structure you want to examine.

## Define paths and fit the model

The quickstart created `scaled`, a moment object for traits `A`, `B` and `C`. This specification says A influences B, and B influences C:

```r
model <- gcmsem_model('
  B ~ ab*A
  C ~ bc*B
', observed = c("A", "B", "C"), gaussian_observed = FALSE)
fit <- gcmsem(model, scaled, n_starts = 8, seed = 20261103)
```

The outcome is left of `~`; its predictor is right of it. `ab` and `bc` name parameters to estimate. The `observed` vector must match the trait names in the moment object. The model adds free trait-specific disturbance variances; undeclared cross-source terms are zero. Inspect `fit_measures(fit)` and `fit$diagnostics` before interpreting paths.

## Choose the source moments

`gaussian_observed = FALSE` frees the marginal third and fourth cumulants of the **observed variables' disturbances**. `gaussian_latent` controls latent disturbances separately. These switches do not force marginal observed variables to be Gaussian when they receive non-Gaussian inputs. `orders = 2:3` in `gcmsem()` selects which observed moments enter the fitting objective; it is a different operation.

Paths, loadings and covariances use familiar SEM notation. Higher cumulants use an explicit tuple of variables, with parameter labels kept outside that tuple:

| Syntax | Meaning |
|---|---|
| `y ~ b*x` | Path from x to y, labeled b |
| `f =~ a*x + a*y` | Factor with equal direct loadings a |
| `f ~~ 1*f` | Fixed factor disturbance variance |
| `x ~~ y` | Disturbance covariance |
| `k3(x, x, y) = free("s")` | Third cumulant of x,x,y disturbances, labeled s |
| `k3(x, y, y) = free("t")` | Third cumulant of x,y,y disturbances, labeled t |
| `k4(x, x, y, y) = free()` | Free fourth cumulant of x,x,y,y disturbances |
| `k3(f, f, f) = free("skew_f")` | Free factor third cumulant |
| `k4(f, f, f, f) = 0` | Fix the factor fourth cumulant to its Gaussian value |

Every argument inside `k3(...)` or `k4(...)` is a variable. The quoted name inside `free("...")` is a parameter label, even when it has the same spelling as a variable. For example, `k3(x, x, y) = free("y")` contains two x indices and one y index; its parameter is labeled y. Reusing a label imposes equality. `free()` or a bare `k3(x, x, y)` creates an unlabeled free parameter; a numeric right-hand side fixes the cumulant.

Repeated names specify repeated tensor indices; fourth-order parameters are **excess cumulants**, so Gaussian means zero. Tuple syntax also represents cross-cumulants involving three or four different traits. Higher cross-source cumulants not declared in syntax are zero. The old directional punctuation aliases have been retired rather than assigned a new meaning. See the [model syntax guide](../docs/model-syntax.md) for defaults and migration guidance. This is a documented subset and extension of lavaan-like syntax, not a drop-in replacement for every lavaan feature. Means are estimated by Python but are not currently fitted as an R mean structure.

Fit the specification with `fit <- gcmsem(model, moments)`. Follow [Working with R results](R-results.md) to inspect model fit before interpreting paths, or the [causal examples](Causal-examples.md) for complete specifications with path diagrams.
