# GCMSEM — Genomic Comoment SEM

GCMSEM uses GWAS summary statistics to study how traits share genetic effects and to fit structural equation models (SEMs) containing directed paths, latent common causes, or both. Users specify their own models and choose which higher moments to include.

A genetic covariance describes whether SNP effects on two traits tend to vary together. Higher moments add information about the **shape of their joint distribution**: co-skewness describes joint asymmetry, and fourth-order co-moments describe how large effects occur together. Two causal models can imply the same covariance while predicting different higher moments. Those differences can help distinguish models when their structural and distributional assumptions hold.

This idea has a history. Non-Gaussian causal methods such as [LiNGAM (Shimizu et al., 2006)](https://jmlr.org/papers/v7/shimizu06a.html) use distributional information beyond covariance to identify linear causal structure under independent disturbances and no hidden confounding. In genetics, [Luke O'Connor and Alkes Price's LCV model (2018)](https://doi.org/10.1038/s41588-018-0255-0) uses mixed fourth moments of SNP effects to distinguish genetic correlation from partial genetic causality. Our earlier [MCMSEM work (Tamimy et al., 2022)](https://doi.org/10.31235/osf.io/ynam2) models covariance, co-skewness and co-kurtosis to examine causal direction, reciprocal paths and latent confounding under specified assumptions.

GCMSEM brings higher-moment SEM to GWAS summary statistics: it estimates moments of the underlying SNP-effect distribution using LD information, then lets the user specify an SEM for those moments. Its directed paths can represent causal hypotheses; their interpretation depends on the model's assumptions and identification. The [causal tutorial](https://github.com/MichelNivard/GCMSEM/wiki/Causal-examples) walks through those decisions.

## Two packages, one analysis

| Component | What it does | How you use it |
|---|---|---|
| **Python package** | Constructs signed LD-power scores, aligns GWAS summary statistics, estimates SNP-effect moments and their uncertainty, and assembles covariance and correlation matrices for any number of traits. | Run the `gcmsem` command in a terminal, or call it from R through `estimate_moments()`. |
| **R package** | Defines SEMs, propagates source moments through paths and factors, fits parameters, and returns model diagnostics, uncertainty and readable summaries. | Use `gcmsem_model()`, `gcmsem()` and ordinary R objects. Export tables only when requested. |

The workflow resembles [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM): estimate genetic moments and their sampling covariance, then fit a model in R. The use of higher moments follows [MCMSEM](https://github.com/MichelNivard/MCMSEM). GCMSEM has its own parser and fitting engine; its familiar paths and factors are extended with explicit `k3(...)` and `k4(...)` source-cumulant specifications.

## From summary statistics to a fitted model

### 1. Obtain LD-power scores

Nearby SNPs are correlated through linkage disequilibrium (LD), so a GWAS association reflects effects at several variants. GCMSEM uses sums of signed LD correlations raised to powers one through four to estimate moments of the underlying SNP-effect distribution. It can compute these scores from PLINK reference genotypes or supplied signed LD blocks using `gcmsem score`; see [score construction](docs/scores.md).

For the worked examples, use the **precomputed European scores at HapMap3 tag SNPs** described on [Data and scores](https://github.com/MichelNivard/GCMSEM/wiki/Data-and-scores). They cover 1,074,062 regression SNPs, with 1,719,440 reference source variants across five minor-allele-frequency strata. The source variants are the variants whose effects contribute to the associations; regression SNPs are the rows used to estimate those effects' moments. Download the precomputed scores from [Zenodo, DOI 10.5281/zenodo.23977299](https://zenodo.org/records/23977299). The [download page](https://github.com/MichelNivard/GCMSEM/wiki/Home#data-release-status) gives checksum-verified commands for the scores and GWAS files. It contains derived scores and metadata, with no reference genotypes or original HapMap3 dataset.

### 2. Clean and align the GWAS files

Call R's [`estimate_moments()`](man/estimate_moments.Rd) with one input file per trait and the score file. Its Python [`harmonize_sumstats()`](python/gcmsem/io.py) routine:

1. Checks the required SNP, allele, effect and standard-error columns, and rejects duplicate SNP identifiers.
2. Matches SNPs to the score tags; removes invalid numbers, nonpositive SEs, invalid alleles and mismatches. Optional sample-size and imputation-quality thresholds are recorded in the log.
3. Aligns effect alleles to the reference, changing the sign of an effect when the alleles are reversed. It checks allele frequencies when supplied. Ambiguous A/T and C/G SNPs are dropped by default; retaining them requires frequencies that uniquely resolve their orientation.
4. Converts declared per-allele effects and SEs to standardized-genotype units using reference allele frequencies. The phenotype units remain unchanged.
5. Retains the common valid SNP set across all traits and records per-trait exclusions and the final count. Zero-effect SNPs are retained; there is no association-significance filter.

The [causal tutorial](https://github.com/MichelNivard/GCMSEM/wiki/Causal-examples) shows the input columns and original-download preparation commands. The [15-trait tutorial](https://github.com/MichelNivard/GCMSEM/wiki/Fifteen-trait-correlations) builds a manifest and runs this procedure for a full correlation matrix.

### 3. Estimate moments and their uncertainty

`estimate_moments()` estimates means, covariances, third cumulants and fourth cumulants, up to the requested order. A cumulant separates the information at its order from products of lower-order moments; Gaussian distributions have zero third and fourth cumulants. Set `max_order = 2` for genetic correlations or `max_order = 4` for the higher-moment SEM examples.

GWAS estimates contain sampling error, which can be correlated across traits sharing participants. The engine estimates that error covariance from LD regressions and corrects products of GWAS estimates under a Gaussian sampling-error model. It repeats the estimation while leaving out genomic groups to obtain a joint sampling covariance, `V`. You receive the moments, covariance matrix `S`, correlation matrix `R`, their SEs, `V`, and timestamped logs in one R object. See [the estimator guide](docs/moment-estimator.md) for the equations and assumptions.

### 4. Fit and assess an SEM

In R, `gcmsem_model()` defines paths, factors and source moments; `gcmsem()` estimates their parameters by **diagonally weighted least squares (DWLS)**. It minimizes the sum of squared differences between estimated and model-predicted moments, each divided by that moment's sampling variance. Exact derivatives of the predicted moments make the repeated optimizer calculations efficient and avoid finite-difference approximations. Parameter SEs use the full `V`, including correlations between estimated moments.

The fit reports the number of moments and parameters, the rank and condition number of the model Jacobian (how moments change with parameters), the remaining gradient, optimizer convergence and results from multiple starts. Too few independent moment changes reveal local nonidentification; a large condition number indicates weak numerical separation of parameters. These checks assess the selected solution. They do not establish global uniqueness or causal direction.

Source distributions receive separate checks: negative variances, whether the covariance matrix is positive semidefinite, and the relationship between skewness and kurtosis. For standardized skewness $s$ and excess kurtosis $k$, every distribution with finite fourth moment must satisfy $k \geq s^2-2$. GCMSEM also checks the joint fourth-order moment matrix for models with at most 12 source nodes. Passing these necessary checks is useful evidence, but does not prove a complete multivariate distribution exists.

**Classical and higher-order Heywood cases are retained and prominently warned about.** The default does not constrain or discard them: sampling noise can push otherwise meaningful estimates outside the admissible region. An optional feasibility penalty changes the estimator and must be requested explicitly. [Working with R results](https://github.com/MichelNivard/GCMSEM/wiki/R-results) explains these flags, the Q fit discrepancy, its approximate SE and the exploratory model-fit reference.

## Effective polygenicity from one or more GWAS

Use [`estimate_polygenicity()`](man/estimate_polygenicity.Rd) with GWAS files and the
GCMSEM score archive to estimate effective variant counts directly in R. Each
trait uses its own eligible SNPs. The returned object contains counts and SEs,
source cumulants, intermediate quantities, joint sampling covariance, all genomic
deletions, influential-block diagnostics and logs. `summary()` explains the
results; `plot()` displays denominator sensitivity. The [polygenicity tutorial](https://github.com/MichelNivard/GCMSEM/wiki/Polygenicity)
starts with a standalone analysis and documents every argument, then validates
the method against S-LD4M with the author's reference inputs.

## Install

Use Python 3.10 or newer, R 4.1 or newer and a C++ compiler. In a terminal:

```bash
git clone https://github.com/MichelNivard/GCMSEM.git
cd GCMSEM
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Start R with that Python environment activated, or pass its executable to `estimate_moments(python = ...)`. The Python package installs NumPy and SciPy; R uses `jsonlite` for communication with it. R-only fitting from existing moments needs no Python process. Try the [small complete example](https://github.com/MichelNivard/GCMSEM/wiki/Quickstart), which creates its own inputs.

## Wiki guide

| Page | What you will learn |
|---|---|
| [Home / installation and downloads](https://github.com/MichelNivard/GCMSEM/wiki) | Install both packages, obtain the external datasets and choose a tutorial. |
| [Small complete example](https://github.com/MichelNivard/GCMSEM/wiki/Quickstart) | Generate three synthetic GWAS files and run an entire analysis in R. |
| [Model specification](https://github.com/MichelNivard/GCMSEM/wiki/Model-specification) | Define directed paths, factors, equality constraints, and Gaussian or free higher source moments. |
| [Working with R results](https://github.com/MichelNivard/GCMSEM/wiki/R-results) | Inspect objects, assess Q and its approximate uncertainty, interpret warnings and save results. |
| [Height/weight and BMI/SBP](https://github.com/MichelNivard/GCMSEM/wiki/Causal-examples) | Download GWAS files, estimate moments and compare illustrated causal models before interpreting paths. |
| [Fifteen-trait correlations](https://github.com/MichelNivard/GCMSEM/wiki/Fifteen-trait-correlations) | Compute genetic correlations from summary statistics, run standard LDSC and compare 105 trait pairs in an R scatterplot. |
| [Kurtosis and polygenicity](https://github.com/MichelNivard/GCMSEM/wiki/Polygenicity) | Run `estimate_polygenicity()` for height, use depression to diagnose weak precision, then study the 15-trait S-LD4M validation comparison. |
| [Data and scores](https://github.com/MichelNivard/GCMSEM/wiki/Data-and-scores) | Understand the LD reference, score contents, source datasets and checksum-verified downloads. |
| [Standard LDSC environment](https://github.com/MichelNivard/GCMSEM/wiki/Standard-LDSC-environment) | Install the original LDSC software and its standard reference and weights for the comparison. |
| [Execution report](https://github.com/MichelNivard/GCMSEM/wiki/Execution-report) | Inspect the actual analyses, numerical checks and remaining validation limitations. |
| [Rights and attribution](https://github.com/MichelNivard/GCMSEM/wiki/DATA_RIGHTS) | Find the terms and attribution attached to each external dataset. |

## Validation and development

GCMSEM is new software. The [executed validation](docs/validation.md) includes known-model recovery, analytical-derivative checks and GWASBrewer simulations with realistic LD and unequal blocks. Higher-order interval coverage still needs improvement; the [statistical assumptions](docs/statistical-assumptions.md) explain the conditions behind moment estimation, SEM identification and uncertainty. R package 0.4.0 and Python package 0.2.0 have separate version numbers because they provide different parts of the workflow.

For the implementation and its RAM-style derivative structure, see [the engine guide](docs/engine.md). See [references](docs/references.md) and [contributing](CONTRIBUTING.md) for provenance and development checks.
