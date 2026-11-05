# Genetic correlations for 15 traits in R

GCMSEM is new software, and genetic correlation provides a useful first comparison with an established method. In this tutorial you will compute a 15 × 15 genetic-correlation matrix from GWAS summary statistics, run standard LDSC on the same input traits, read its results into R, and compare all 105 distinct pairs in a scatterplot. The goal is to check empirical agreement between the procedures and learn the complete analysis workflow.

The comparison uses **original, unmodified bulik/ldsc 1.0.1 with its published reference and weighting scores**. Agreement here checks the second-order estimates. It does not establish the validity of higher-order causal models. Both methods retain their own ordinary filtering, weighting and uncertainty calculations.

## 1. Install the software

Run these commands in a terminal. They download the source code, install the Python engine, and install the R package:

```bash
git clone https://github.com/MichelNivard/GCMSEM.git
cd GCMSEM
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Use Python 3.10 or newer for GCMSEM. Start R from this `GCMSEM` directory, or open it as your RStudio working directory. Terminal commands below also assume this directory.

## 2. Get the scores and summary statistics

The inputs are available from [Zenodo, DOI 10.5281/zenodo.23977299](https://zenodo.org/records/23977299). Download and verify the two input archives with:

```bash
python tutorials/download_data.py \
  --asset scores --asset panel --out tutorial-data
```

If you already have `gcmsem-scores-v1.zip` and `gcmsem-panel-v1.zip`, put them in the working directory and extract them:

```bash
mkdir -p tutorial-data
unzip gcmsem-scores-v1.zip -d tutorial-data
unzip gcmsem-panel-v1.zip -d tutorial-data
```

The numerical data live in these external archives. The repository contains the analysis scripts and the tutorial figure.

| File inside `tutorial-data` | Contents and use |
|---|---|
| `scores/eur_hm3.scores.npz` | GCMSEM's precomputed signed LD scores, SNP alleles, positions, MAF strata, and source-variant counts. The tag SNPs use the HM3 selection. |
| `sumstats/*.tsv.gz` | One compressed, allele-aligned GWAS table per trait, with `SNP`, `A1`, `A2`, `BETA`, `SE`, `N`, `CHR`, `BP`, and `REF_AF`. |
| `panel.npz` | The exact prealigned SNP-by-trait arrays used for the archived 15-trait GCMSEM run. |

`BETA` and `SE` in these prepared GWAS files are on the genotype-standardized scale; `A1` is the effect allele and `A2` the other allele. `REF_AF` is the reference-panel frequency of `A1`. `N` is the SNP-specific sample size. The score archive contains derived scores and metadata; reference genotypes are not distributed.

The 15 input files are:

| R trait name | GWAS filename in `sumstats/` | Phenotype |
|---|---|---|
| `DrinksPerWeek` | `GSCAN_DrinksPerWeek.tsv.gz` | GSCAN drinks per week, public release excluding 23andMe |
| `VitaminD` | `30890_raw.tsv.gz` | Vitamin D |
| `SBP` | `4080_irnt.tsv.gz` | Systolic blood pressure |
| `Insomnia` | `1200.tsv.gz` | Insomnia |
| `Depression` | `20544_11.tsv.gz` | Diagnosed depression |
| `LDL` | `30780_raw.tsv.gz` | Direct LDL |
| `FEV1` | `3063_irnt.tsv.gz` | Forced expiratory volume in one second |
| `HeartAttack` | `6150_1.tsv.gz` | Heart attack |
| `Angina` | `6150_2.tsv.gz` | Angina |
| `Height` | `50_irnt.tsv.gz` | Standing height |
| `Weight` | `21002_irnt.tsv.gz` | Weight |
| `BMI` | `21001_irnt.tsv.gz` | Body mass index |
| `DBP` | `4079_irnt.tsv.gz` | Diastolic blood pressure |
| `FVC` | `3062_irnt.tsv.gz` | Forced vital capacity |
| `Waist` | `48_irnt.tsv.gz` | Waist circumference |

The numbered phenotypes come from the Neale Lab UK Biobank release. `_irnt` denotes an inverse-rank-normalized phenotype. This tutorial uses the IRNT BMI and SBP files; the causal tutorial uses their raw-phenotype files. See [Data and scores](Data-and-scores.md) for the source releases and preparation details.

## 3. Compute genetic correlations in R

Start R with the GCMSEM Python environment activated. Load the R package and the plotting helper:

```r
library(GCMSEM)
options(digits = 6)
source("tutorials/correlations.R")
data_dir <- "tutorial-data"
```

Create a manifest with one row per trait. The names below become the row and column names of your matrices; the codes identify the downloaded files.

```r
codes <- c(
  DrinksPerWeek = "GSCAN_DrinksPerWeek", VitaminD = "30890_raw",
  SBP = "4080_irnt", Insomnia = "1200", Depression = "20544_11",
  LDL = "30780_raw", FEV1 = "3063_irnt", HeartAttack = "6150_1",
  Angina = "6150_2", Height = "50_irnt", Weight = "21002_irnt",
  BMI = "21001_irnt", DBP = "4079_irnt", FVC = "3062_irnt", Waist = "48_irnt"
)
manifest <- data.frame(
  trait = names(codes),
  path = file.path(data_dir, "sumstats", paste0(unname(codes), ".tsv.gz")),
  effect_scale = "genotype_standardized"
)
stopifnot(all(file.exists(manifest$path)))

python <- file.path(getwd(), ".venv", "bin", "python")
gcmsem <- estimate_moments(
  sumstats = manifest,
  scores = file.path(data_dir, "scores/eur_hm3.scores.npz"),
  max_order = 2,
  n_blocks = 200,
  python = python,
  verbose = TRUE
)
gcmsem
```

This call performs the analysis. It reads and harmonizes every GWAS file, converts effects to the declared scale, estimates SNP-effect means and covariances using the supplied LD scores, and refits the calculation in 200 genomic deletions to estimate uncertainty. `max_order = 2` requests the moments needed for correlations. The result is an R object; no permanent output files are created unless requested.

The normal harmonizer drops unresolved palindromic SNPs. `REF_AF` is a reference frequency and is not substituted for an observed GWAS frequency. In this 15-trait intersection the final common panel has 1,012,060 SNPs; the potentially affected tags were already absent from the intersection. The R reference run took about 14 minutes. The full reference source counts are preserved when tag SNPs are filtered. The reference run prints:

```text
GCMSEM moments
Traits: DrinksPerWeek, VitaminD, SBP, Insomnia, Depression, LDL, FEV1, HeartAttack, Angina, Height, Weight, BMI, DBP, FVC, Waist
SNPs: 1,012,060
Cumulants: 135 (order 1: 15, order 2: 120)
Sampling covariance V: 135 x 135   Jackknife groups: 200

Genetic correlations ($R)
              DrinksPerWeek VitaminD     SBP Insomnia Depression     LDL
DrinksPerWeek        1.0000   0.0870 -0.0615   0.1024     0.0579 -0.0724
VitaminD             0.0870   1.0000  0.0483   0.0113    -0.1821 -0.1110
SBP                 -0.0615   0.0483  1.0000   0.0435    -0.0592 -0.1594
Insomnia             0.1024   0.0113  0.0435   1.0000     0.3999 -0.0156
Depression           0.0579  -0.1821 -0.0592   0.3999     1.0000  0.0463
LDL                 -0.0724  -0.1110 -0.1594  -0.0156     0.0463  1.0000
Showing 6 of 15 traits; the full matrix is in $R.
```

Extract complete matrices using `$`:

```r
gcmsem_cov <- gcmsem$S       # SNP-effect covariances: 15 by 15
gcmsem_rg <- gcmsem$R        # Correlations: 15 by 15
gcmsem_rg_se <- gcmsem$R_SE  # Delta-method standard errors
gcmsem_v <- gcmsem$V         # Joint sampling covariance of all 135 cumulants

anthropometric <- c("Height", "Weight", "BMI")
round(gcmsem_rg[anthropometric, anthropometric], 3)
```

```text
       Height Weight    BMI
Height  1.000  0.470 -0.181
Weight  0.470  1.000  0.783
BMI    -0.181  0.783  1.000
```

The correlation between traits is their estimated covariance divided by the product of their estimated SNP-effect SDs. Its SE propagates uncertainty from those three quantities using the full `V`. `summary(gcmsem)` shows matrices with SEs; `print(gcmsem, max.traits = Inf)` displays every trait. The log remains in `gcmsem$log`, and the deletion estimates are in `gcmsem$jackknife`.

Save the completed R analysis before switching back to the terminal:

```r
saveRDS(gcmsem, "gcmsem-correlations.rds")
```

## 4. Run standard LDSC

Set up original Python 2.7 LDSC using [Standard LDSC environment](Standard-LDSC-environment.md). That page defines the shell variable `gcmsem_ldsc_python`. It is a separate environment from the Python 3 environment used by GCMSEM.

In the terminal, fetch the pinned original LDSC checkout and its standard reference files:

```bash
python tutorials/setup_standard_ldsc.py \
  --out external/standard-ldsc --python2 "$gcmsem_ldsc_python"
python tutorials/download_standard_ldsc.py --out external/standard-reference
```

The software is bulik/ldsc commit `2fdeeb3b44379408794154993dbd6101b8946b7e`. `external/standard-reference/w_hm3.snplist` lists SNPs and alleles for the standard `--merge-alleles` step. `eur_w_ld_chr/` contains the published European LD scores, regression weights and source counts. The same published prefix is used for `--ref-ld-chr` and `--w-ld-chr`, following the [official LDSC tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

Prepare the ordinary LDSC input columns from the same downloaded GWAS files:

```bash
python tutorials/prepare_standard_ldsc_inputs.py \
  --sumstats tutorial-data/sumstats --out tutorial-data/standard-inputs
```

Each output contains `SNP A1 A2 BETA P N`. The adapter computes a two-sided normal `P` from `BETA / SE`. Stock `munge_sumstats.py` then performs its usual filtering and allele checks. GCMSEM's existing HM3 SNP selection is retained in each trait input, as intended for this comparison. LDSC takes its ordinary pair-specific intersections.

Run munging and all 105 regressions:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tutorials/run_standard_ldsc.py \
  --inputs tutorial-data/standard-inputs \
  --hm3 external/standard-reference/w_hm3.snplist \
  --ref-prefix external/standard-reference/eur_w_ld_chr \
  --ldsc-dir external/standard-ldsc \
  --python2 "$gcmsem_ldsc_python" \
  --out tutorial-results/standard_ldsc
```

This runner invokes original `munge_sumstats.py` and `ldsc.py`, retaining their statistical defaults. Its only I/O tuning is a 500,000-row munging batch size, independently checked to leave the decompressed munged results identical. Every command and native LDSC log is retained under the output directory.

For example, after munging, height/weight uses this ordinary LDSC command:

```bash
"$gcmsem_ldsc_python" external/standard-ldsc/ldsc.py \
  --rg tutorial-results/standard_ldsc/munged/50_irnt.sumstats.gz,tutorial-results/standard_ldsc/munged/21002_irnt.sumstats.gz \
  --ref-ld-chr external/standard-reference/eur_w_ld_chr/ \
  --w-ld-chr external/standard-reference/eur_w_ld_chr/ \
  --out tutorial-results/standard_ldsc/height_weight_example
```

The main runner creates:

| Output | Contents |
|---|---|
| `ldsc_standard_rg.tsv` | The 15 × 15 genetic-correlation matrix collected from stock LDSC summaries. |
| `ldsc_standard_rg_se.tsv` | The corresponding standard errors. |
| `all_pairs.tsv` | One row per pair, including correlations, SEs and SNP counts. |
| `munged/` | The formatted summary statistics and munging logs. |
| `fits/` | Individual native LDSC regression outputs and logs. |
| `run.log`, `configuration.json` | Commands, settings and source/reference checksums. |

The reference run completes all 105 pairs, using 1,012,052–1,057,882 SNPs per pair. The matrices retain stock LDSC's printed four-decimal-place precision.

## 5. Read the LDSC correlations into R

Return to R. If you started a new session, reload the GCMSEM object with `gcmsem <- readRDS("gcmsem-correlations.rds")`, load `GCMSEM`, source `tutorials/correlations.R`, and set `gcmsem_rg <- gcmsem$R`.

Read LDSC's two tables as named numeric matrices:

```r
ldsc_dir <- "tutorial-results/standard_ldsc"
ldsc_rg <- as.matrix(read.delim(
  file.path(ldsc_dir, "ldsc_standard_rg.tsv"),
  row.names = 1, check.names = FALSE
))
ldsc_rg_se <- as.matrix(read.delim(
  file.path(ldsc_dir, "ldsc_standard_rg_se.tsv"),
  row.names = 1, check.names = FALSE
))

traits <- rownames(gcmsem_rg)
stopifnot(setequal(traits, rownames(ldsc_rg)),
          setequal(traits, colnames(ldsc_rg)),
          setequal(traits, rownames(ldsc_rg_se)),
          setequal(traits, colnames(ldsc_rg_se)))
ldsc_rg <- ldsc_rg[traits, traits]
ldsc_rg_se <- ldsc_rg_se[traits, traits]

anthropometric <- c("Height", "Weight", "BMI")
round(ldsc_rg[anthropometric, anthropometric], 3)
```

```text
       Height Weight    BMI
Height  1.000  0.428 -0.169
Weight  0.428  1.000  0.818
BMI    -0.169  0.818  1.000
```

`row.names = 1` uses the first column as row names; `check.names = FALSE` preserves the trait names. Name-based indexing makes both methods use the same order. Each matrix contains every pair twice, so the comparison will take only one triangle and omit the diagonal.

## 6. Compare all 105 pairs and draw the scatterplot

The functions in [`correlations.R`](https://github.com/MichelNivard/GCMSEM/blob/main/tutorials/correlations.R) create a comparison object and draw it using base R graphics:

```r
comparison <- compare_correlations(gcmsem, ldsc_rg, ldsc_rg_se)
comparison
```

```text
GCMSEM and standard LDSC correlation comparison
15 traits; 105 distinct trait pairs

Pearson correlation            0.982855
RMSE                           0.0512594
Mean GCMSEM - LDSC             -0.0036361
Largest absolute difference    0.155582
```

`comparison$pairs` is a data frame with one row per pair. `difference` is GCMSEM minus LDSC:

```r
head(comparison$pairs, 3)
```

```text
         trait1   trait2     gcmsem    ldsc  difference ldsc_se
1 DrinksPerWeek VitaminD  0.0869961  0.0539  0.03309615  0.0197
2 DrinksPerWeek      SBP -0.0614626 -0.0342 -0.02726258  0.0200
3 DrinksPerWeek Insomnia  0.1023968  0.0937  0.00869675  0.0248
```

Draw the annotated plot in your R graphics window:

```r
plot(comparison)
```

![Scatterplot of the 105 trait-pair correlations: standard LDSC on the horizontal axis, GCMSEM on the vertical axis, with a dashed equality line.](https://raw.githubusercontent.com/MichelNivard/GCMSEM/main/tutorials/figures/correlation-scatterplot.png)

Each point is one trait pair. The dashed line marks equal estimates. The orange point labels Heart Attack/Angina, the largest difference in this analysis. Pearson correlation measures how similarly the estimates vary across pairs; RMSE measures their typical separation on the correlation scale. Pairs share traits and the methods share GWAS inputs, so these are descriptive agreement statistics.

To make your own simple version, the complete plotting code is:

```r
pairs <- comparison$pairs
limits <- range(c(pairs$ldsc, pairs$gcmsem))
plot(pairs$ldsc, pairs$gcmsem,
     xlim = limits, ylim = limits, asp = 1,
     pch = 19, col = "steelblue",
     xlab = "Standard LDSC genetic correlation",
     ylab = "GCMSEM SNP-effect correlation")
abline(a = 0, b = 1, lty = 2, col = "grey40")
```

To save the annotated figure, explicitly open a graphics file:

```r
png("correlation-scatterplot.png", width = 1600, height = 1500, res = 220)
plot(comparison)
dev.off()
```

Your estimates remain in `gcmsem`, `ldsc_rg`, `ldsc_rg_se`, and `comparison`. You can save those R objects for a later session:

```r
saveRDS(list(gcmsem = gcmsem, comparison = comparison), "correlation-analysis.rds")
```

TSV export is optional:

```r
comparison <- compare_correlations(
  gcmsem, ldsc_rg, ldsc_rg_se,
  write.to.tsv = TRUE, out = "tutorial-results/correlation-comparison"
)
```

This writes `.pairs.tsv` and `.statistics.tsv`, and still returns `comparison`. The [executed console transcript](https://github.com/MichelNivard/GCMSEM/blob/main/tutorials/output/correlation-comparison.txt) records the recomputed values used in this figure.

## 7. Optional archived results and exact historical reproduction

The primary workflow above estimates GCMSEM correlations from GWAS files and runs standard LDSC. To inspect the previously distributed numerical results, download the optional results archive:

```bash
python tutorials/download_data.py --asset results --out tutorial-data
```

`results/standard_ldsc_verified/` contains the completed stock LDSC run used as the reference here. To compare with those results while keeping your newly computed GCMSEM object, set `ldsc_dir <- "tutorial-data/results/standard_ldsc_verified"` and repeat sections 5–6.

`results/correlations/gcmsem.*` contains the historical GCMSEM moment bundle. That archived analysis used a prealigned panel of 1,012,060 SNPs. The normal GWAS-file route above produced the same common intersection and reproduced the correlations in the executed check. You can load that specific historical result with `read_moments()` or recompute it from the archived `panel.npz`:

```bash
python tutorials/run_correlations.py \
  --panel tutorial-data/panel.npz \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --out tutorial-results/historical-panel --stage gcmsem
```

```r
historical <- read_moments("tutorial-results/historical-panel/gcmsem")
```

That separate reproduction preserves the original allele alignment. The main tutorial's console output and scatterplot use the normal harmonizer and freshly estimated 1,012,060-SNP analysis. Both analyses retain the score source universe of 1,719,440 variants across five MAF strata.

GCMSEM's correlation describes centered SNP effects pooled over those strata. Standard LDSC uses its own published reference, weighting, effect assumptions and jackknife. Strong agreement is useful empirical evidence, while differences should be examined using each procedure's assumptions and uncertainty.
