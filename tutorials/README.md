# GCMSEM empirical tutorials

These worked examples use real GWAS summary statistics and signed LD-power scores at HapMap3 regression SNPs. All numerical inputs and results are distributed separately from GitHub. No reference genotypes, PLINK BED/BIM/FAM files, or original HapMap3 dataset are included in the tutorial archives.

| Tutorial | What it runs |
|---|---|
| [Height/weight and BMI/SBP](Causal-examples.md) | Summary statistics → cumulants and joint sampling covariance → standalone R SEM, including retained Heywood cases |
| [Fifteen-trait correlations](Fifteen-trait-correlations.md) | All 105 pairs with GCMSEM and untouched stock LDSC using its published standard reference and weights |
| [Data and scores](Data-and-scores.md) | Reference definition, download integrity, source reconstruction and reuse limits |
| [Rights and attribution](DATA_RIGHTS.md) | Original GWAS policies and reference attribution |

## Data release status

The numerical release is prepared for Zenodo, and account authorization is complete. **The upload is incomplete; the record is not yet published.** Zenodo has reserved DOI `10.5281/zenodo.23977299` for the draft. This is not yet a working public data release. The checked-in [manifest](data_manifest.json) records exact archive sizes and SHA-256 checksums. Once publication is verified, its `zenodo_record_id` will be set to the actual published record.

After publication, run from the repository root:

```bash
python tutorials/download_data.py --out tutorial-data
```

If using a supplied, already-downloaded archive before deposition, extract it into `tutorial-data/`, preserving the paths in the archive. The examples then run locally without Zenodo credentials. A replacement host must serve the same pinned bytes; updates require a new data manifest.

## Install

Use Python 3.12 or newer, R and a C++ compiler. A separate Python environment avoids changing an existing analysis setup:

```bash
git clone https://github.com/MichelNivard/GCMSEM.git
cd GCMSEM
python -m venv .venv
source .venv/bin/activate
python -m pip install .
R CMD INSTALL .
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
```

The standard LDSC comparison uses an isolated Python 2.7 environment for untouched original `bulik/ldsc`. Follow [Standard-LDSC-environment.md](Standard-LDSC-environment.md), then the [15-trait tutorial](Fifteen-trait-correlations.md). GCMSEM continues to run under Python 3; do not install old LDSC dependencies into the GCMSEM environment.

## What is reproduced

The scores use a European 503-person reference, GRCh37 coordinates, five disjoint MAF strata and a 1 Mb rolling LD window. The 1,074,062 regression tags and 1,719,440 source variants have different roles: filtering a GWAS tag does not change source counts. Scores and effects use the same minor-allele orientation; signed odd powers are retained. GCMSEM's HM3 SNP selection is retained for both methods' inputs; standard LDSC then applies its ordinary munging and pair-specific QC.

The [execution report](Execution-report.md) records the completed real-data runs, all six SEM fits, all 105 standard LDSC pairs and software validation.

The SEM examples use the normal file-based harmonizer. The correlation tutorial uses a documented pre-harmonized 15-trait array to reproduce the exact common panel. The cache lacks original GWAS EAF, so the ordinary harmonizer drops three palindromic tags rather than inventing EAF from reference frequencies. These two routes are explicitly recorded; they must not be mixed silently.

The higher-order analyses remain research examples. Formal sandwich SEs are reported when computable, including for inadmissible source moments. Convergence, identification, admissibility and inference flags accompany every fit. The earlier simulation validation showed limited fourth-order interval coverage; these empirical tutorials do not establish causal identification or resolve that calibration limitation.

## Wiki availability

The repository is public and its [GitHub Wiki](https://github.com/MichelNivard/GCMSEM/wiki) is enabled. These linked Markdown pages provide the versioned tutorial source. The Wiki is readable by everyone; editing is restricted to repository collaborators.
