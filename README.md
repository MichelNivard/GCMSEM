# GCMSEM — Genomic Comoment SEM

GCMSEM estimates covariance, co-skewness and fourth cumulants of joint SNP effects from GWAS summary statistics, then fits structural equation models to those quantities. It contains a Python pipeline and an independent R package. Neither requires MCMSEM or lavaan.

**Version 0.1.0 is an initial research release.** The software includes input validation, reproducible tests, logs and numerical diagnostics. Its higher-order inference is still being validated; passing software tests is not a certification of interval coverage or causal identification. See the [executed validation](docs/validation.md) and [statistical assumptions](docs/statistical-assumptions.md).

## What is included

* Signed LD-power score construction from PLINK reference genotypes or signed LD blocks, disjoint MAF strata, and separate reference/source and regression/tag sets.
* Allele harmonization, explicit effect-scale conversion, estimated cross-trait sampling-error covariance, and Gaussian sampling-noise correction without sample splitting.
* Means and cumulants through order four for arbitrary trait counts, selected higher-order tuples, trait covariance/correlation matrices, and their joint genomic jackknife sampling covariance.
* Standalone R model syntax for paths, factors, covariance, third and fourth source cumulants, fixed parameters and equality labels.
* Unconstrained DWLS fitting with exact derivatives, multiple starts, full-covariance sandwich SEs, prominent Heywood warnings, and separate checks for identification, source-moment admissibility and convergence. Source-moment feasibility penalties are optional.

## Install

From the repository root:

```bash
python -m pip install .
R CMD INSTALL .
```

Python requires NumPy and SciPy. R uses its standard libraries and a small C++ kernel; a working C++ compiler is required for source installation. `jsonlite` is optional for richer bundle metadata parsing. GWASBrewer is used only by the separate simulation validation.

## Run a complete small example

```bash
python examples/make_toy_inputs.py --out toy
gcmsem score --manifest toy/ld_manifest.tsv --out toy/scores.npz
gcmsem estimate --sumstats toy/sumstats.tsv --scores toy/scores.npz --noise-intercept toy/noise.tsv --out toy/result
Rscript examples/fit_bundle.R toy/result
```

The toy input uses three traits, unequal synthetic LD block sizes, non-Gaussian effects and correlated sampling error. This quickstart supplies the known generating noise covariance; omit `--noise-intercept` to estimate it from LD regressions, as in an ordinary GWAS analysis. It requires no downloads. The [GWASBrewer validation](validation/) uses the real HapMap3 chromosome 19 LD pattern, actual varying block lengths, and both estimated and supplied overlap correction. The toy example is not a substitute for that validation.

For real GWAS, use an appropriate signed LD reference and an explicit sumstats manifest:

| trait | path | effect_scale |
|---|---|---|
| Height | height.tsv.gz | allele |
| Weight | weight.tsv.gz | allele |
| BMI | bmi.tsv.gz | allele |

Save this as a tab-separated file. Each GWAS table requires `SNP`, `A1`, `A2`, `BETA`, `SE`; optional `EAF` means frequency of `A1`. The manifest scale is `allele` or `genotype_standardized`. The phenotype's original units are preserved. Do not supply odds ratios as `BETA`; use the underlying log-odds coefficient where appropriate. Binary GWAS are not automatically transformed to liability scale.

See [score construction](docs/scores.md) for reference formats and finite-reference corrections, and the [moment estimator guide](docs/moment-estimator.md) for GWAS inputs, equations, selection of higher moments, and output files. The harmonizer uses a common complete set of tags across traits; it does not combine pair-specific SNP subsets into a joint covariance. It keeps zero-effect SNPs and applies no association-significance filter.

## Specify a model in R

```r
library(GCMSEM)
data <- read_moments("toy/result")

model <- gcmsem_model('
  B ~ ab*A
  C ~ bc*B
', observed = c("A", "B", "C"), gaussian_observed = FALSE)

fit <- gcmsem(model, data, n_starts = 8, log_file = "toy/sem.log")
summary(fit)
parameter_estimates(fit)
```

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

Repeated names specify repeated tensor indices; fourth-order parameters are **excess cumulants**, so Gaussian means zero. Tuple syntax also represents cross-cumulants involving three or four different traits. Higher cross-source cumulants not declared in syntax are zero. The old directional punctuation aliases have been retired rather than assigned a new meaning. See the [model syntax guide](docs/model-syntax.md) for defaults and migration guidance. This is a documented subset and extension of lavaan-like syntax, not a drop-in replacement for every lavaan feature. Means are estimated by Python but are not currently fitted as an R mean structure.

## Output and uncertainty

Python writes a labeled, versioned bundle:

| File suffix | Contents |
|---|---|
| `.moments.tsv` | Canonical moment keys, estimates and jackknife SEs |
| `.sampling_cov.tsv` | Full joint sampling covariance, with matching row/column keys |
| `.jackknife.tsv` | Genomic deletion estimates |
| `.covariance.tsv` | SNP-effect covariance matrix across all traits |
| `.correlation.tsv` | SNP-effect correlation matrix across all traits |
| `.metadata.json` | Versions, settings, inputs, hashes, QC and statistical assumptions |
| `.log` | Timestamped commands, stage progress, intercepts and diagnostics |

Multiple MAF strata additionally produce conditional bundles. Pooling uses raw-moment mixture identities, not a simple average of fourth cumulants. Source counts remain fixed when tag rows are filtered or deleted. Jackknife groups preserve whole LD blocks and account for unequal deletion sizes.

R minimizes a diagonal-weighted discrepancy. Its sandwich SEs use the entire supplied sampling covariance, without dividing by sample size again. `Q` is reported as an **uncalibrated discrepancy**, not a WLSMV or chi-square test.

**The default retains Heywood cases.** Free disturbance variances and higher cumulants are not clipped to the admissible region, and the multistart search does not prefer a worse fit merely because its source moments are admissible. Negative variances or impossible higher source moments generate prominent warnings and diagnostic flags; point estimates and computable formal sandwich SEs remain available. This preserves the estimator's sampling variability across replications. Inadmissibility is distinct from nonconvergence or nonidentification, for which an ordinary numerical SE may not be available. Report all statuses when aggregating studies; silently discarding flagged fits can introduce selection bias.

`feasibility="penalty"` is an explicitly requested alternative estimator. It can change sampling bias and the distribution of estimates, particularly near a boundary. Neither formal SEs nor an admissible solution resolve model-selection uncertainty, calibrated interval coverage, or equivalent causal structures.

## Scaling and computation

The count of all symmetric moments through order four is $\binom{N+4}{4}-1$. Python computes higher-order responses one tuple at a time, supports selecting tuples, and checks expected memory requirements. The full joint sampling covariance still needs quadratic storage in the number of selected moments. Its rank is at most $B-1$ for $B$ deletion groups; DWLS does not require inverting that matrix.

For the path matrix $A$, the R engine uses

$$
T=(I-A)^{-1},\qquad \mathrm{d}T=T\,(\mathrm{d}A)\,T.
$$

It then applies product rules to sparse source cumulants. Fitting uses these exact gradients with `nlminb`. This is an exact specialized derivative engine, not a general torch/autodiff graph. Native tensor contraction accelerates the repeated calculation, while the R implementation provides an independent reference path. The [RAM-style engine explanation](docs/engine.md) describes how this generalizes to user-defined models. Large-scale performance and calibration should be assessed for the requested trait/moment set.

## Development and provenance

```bash
python -m pip install '.[test]'
python -m pytest tests/python -q
R CMD build . --no-build-vignettes --no-manual
R CMD check GCMSEM_0.1.0.tar.gz --no-manual --no-vignettes
```

GitHub Actions runs Python tests and R package checks. See [contributing](CONTRIBUTING.md), [references](docs/references.md), [statistical assumptions](docs/statistical-assumptions.md), and [validation](docs/validation.md). All empirical reference and GWAS datasets retain their own licenses and must be obtained separately.
