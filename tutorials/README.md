# GCMSEM tutorials and downloads

GCMSEM estimates genetic moments from GWAS summary statistics and fits user-defined structural equation models in R. Start with installation below, then choose a worked example. Each tutorial describes the files, shows the commands and prints output from an executed run. Numerical data are hosted externally; code, documentation and illustrative figures are kept in this repository.

| Page | Purpose |
|---|---|
| [Small complete example](Quickstart.md) | Generate three synthetic GWAS tables, estimate their moments and fit a model without downloading empirical data. |
| [Model specification](Model-specification.md) | Define paths, factors, constraints and Gaussian or free higher source moments. |
| [Working with R results](R-results.md) | Inspect objects, check fit with Q and its approximate uncertainty, interpret warnings and save results. |
| [Height/weight and BMI/SBP](Causal-examples.md) | Obtain real GWAS summary statistics, estimate moments, draw six models and assess fit before interpreting paths. |
| [Fifteen-trait correlations](Fifteen-trait-correlations.md) | Compute a full genetic-correlation matrix from GWAS files, run standard LDSC, and compare 105 pairs in R. |
| [Kurtosis and polygenicity](Polygenicity.md) | Estimate effective polygenicity for 15 traits, compare original S-LD4M with GCMSEM, and inspect the scatterplot and uncertainty. |
| [Data and scores](Data-and-scores.md) | Understand the precomputed HapMap3-tag scores, input datasets and download checksums. |
| [Standard LDSC environment](Standard-LDSC-environment.md) | Install the original LDSC software and standard reference/weight files. |
| [Execution report](Execution-report.md) | Inspect measured results, checks and validation limits. |
| [Rights and attribution](DATA_RIGHTS.md) | Find the reuse terms and attribution for each external dataset. |

## Install

Use Python 3.10 or newer, R 4.1 or newer and a C++ compiler. Run these commands in a terminal:

```bash
git clone https://github.com/MichelNivard/GCMSEM.git
cd GCMSEM
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

This installs R package version 0.3.0 and Python engine version 0.1.0. `estimate_moments()` uses the Python engine and requires `jsonlite`. Keep the environment activated when starting R. If starting R elsewhere, provide the environment's executable through the function's `python` argument.

The [small worked example](Quickstart.md) generates three synthetic GWAS tables and scores locally. Run it immediately after installation to try `estimate_moments()` → `standardize_moments()` → `gcmsem()` without downloading the empirical datasets.

The standard LDSC comparison uses a separate Python 2.7 environment for original `bulik/ldsc`. Follow [Standard-LDSC-environment.md](Standard-LDSC-environment.md), then the [15-trait tutorial](Fifteen-trait-correlations.md). GCMSEM uses the Python 3 environment above.

## Data release status

The numerical data are publicly available at [Zenodo, DOI 10.5281/zenodo.23977299](https://zenodo.org/records/23977299). All four archives were verified against their expected sizes and checksums. The checked-in [manifest](data_manifest.json) pins their exact bytes and the published record number. No Zenodo login or API token is needed to download them.

From the repository root, download and extract all four archives:

```bash
python tutorials/download_data.py --out tutorial-data
```

For the causal tutorial, request just the scores and four GWAS tables:

```bash
python tutorials/download_data.py --asset scores --asset causal --out tutorial-data
```

For the correlation tutorial, use `--asset scores --asset panel`. Add `--asset results` only if you want the archived numerical results. The downloader checks each archive's size and SHA-256 before extraction. If you download the ZIP files manually from Zenodo, extract them into the same `tutorial-data/` directory, preserving the paths inside.

## What the data release contains

The score archive contains derived signed LD-power scores for 1,074,062 HapMap3 tag SNPs. Its European 503-person reference uses GRCh37 coordinates, five disjoint MAF strata and 1 Mb rolling LD windows. The 1,719,440 source variants have a different role from the regression tags: filtering a tag does not change the source counts. No reference genotypes or original HapMap3 dataset are redistributed.

The causal archive supplies height, weight, raw BMI and raw SBP GWAS subsets. The panel archive supplies the 15 GWAS tables for the correlation analysis and the historical prealigned array. The optional results archive supplies the original numerical bundles and the complete standard LDSC comparison. The updated tutorials recompute GCMSEM analyses and fit diagnostics from those inputs.

[The public Wiki](https://github.com/MichelNivard/GCMSEM/wiki) publishes these pages for browsing. Their Markdown sources remain versioned here. Software validation and limitations are described in the [execution report](Execution-report.md), including the remaining uncertainty about higher-order interval coverage.
