# Height/weight and BMI/blood-pressure SEMs in R

This tutorial starts with SNP-level GWAS summary statistics, estimates their genetic cumulants, and fits an SEM in R. You will first fit a height-to-weight model, inspect its parameters and diagnostics, and then repeat the analysis for BMI and systolic blood pressure (SBP). The complete R scripts also fit the reverse direction and a reciprocal sensitivity model.

The examples use Neale Lab UK Biobank round-2 GWAS results:

| R trait name | Source phenotype | Units |
|---|---|---|
| `Height` | `50_irnt` | Inverse-normal-transformed height |
| `Weight` | `21002_irnt` | Inverse-normal-transformed weight |
| `BMI` | `21001_raw` | kg/m² |
| `SBP` | `4080_raw` | mmHg |

## 1. Install the software and obtain the files

Follow the [installation commands](README.md#install), then run the terminal commands below from the `GCMSEM` repository directory. Keep the Python environment activated when starting R so that `estimate_moments()` can find the installed Python engine.

**Data availability:** the numerical archives are prepared, but public deposition is still pending. The [data release page](README.md#data-release-status) records the current status. If you have the supplied archives, extract them now:

```bash
mkdir -p tutorial-data
unzip gcmsem-scores-v1.zip -d tutorial-data
unzip gcmsem-causal-v1.zip -d tutorial-data
```

Once the data release is published, this command downloads, verifies and extracts the same two archives:

```bash
python tutorials/download_data.py --asset scores --asset causal --out tutorial-data
```

Here are the files used in this tutorial:

| File under `tutorial-data/` | Contents |
|---|---|
| `scores/eur_hm3.scores.npz` | Precomputed signed LD-power scores and reference metadata for 1,074,062 HapMap3 regression SNPs; orders 1–4, five MAF strata and 1 Mb LD windows |
| `sumstats/50_irnt.tsv.gz` | Prepared height summary statistics |
| `sumstats/21002_irnt.tsv.gz` | Prepared weight summary statistics |
| `sumstats/21001_raw.tsv.gz` | Prepared BMI summary statistics |
| `sumstats/4080_raw.tsv.gz` | Prepared SBP summary statistics |
| `provenance/` | Original source URLs, phenotype definitions, checksums and preparation records |

The compressed summary-statistic files are ordinary tab-separated text. The `.npz` score file is read directly by the estimator. The [scores page](Data-and-scores.md) explains how these scores were constructed.

If you want to start from the original GWAS downloads, the [final section](#9-start-from-the-original-neale-lab-downloads) gives the exact `wget` and preparation commands.

## 2. Inspect the height file in R

Start R in the repository directory:

```r
library(GCMSEM)

height_head <- read.delim(
  gzfile("tutorial-data/sumstats/50_irnt.tsv.gz"),
  nrows = 3
)
print(height_head, digits = 4, row.names = FALSE)
```

The first three records are:

```text
        SNP A1 A2       BETA       SE      N CHR     BP  REF_AF
  rs3131969  A  G -0.0009732 0.001167 360388   1 754182 0.12823
 rs12562034  A  G  0.0011891 0.001097 360388   1 768448 0.09245
  rs4040617  G  A -0.0009679 0.001128 360388   1 779322 0.11829
```

Each row describes one SNP:

| Column | Meaning in the prepared files |
|---|---|
| `SNP` | Variant rs identifier, used to match the score file |
| `A1`, `A2` | Effect allele and other allele; effects have been aligned to the score reference |
| `BETA`, `SE` | Association estimate and standard error per SD of genotype dosage |
| `N` | GWAS sample size for that SNP |
| `CHR`, `BP` | Chromosome and GRCh37 base-pair coordinate |
| `REF_AF` | Reference-panel frequency of `A1`; the original GWAS effect-allele frequency is unavailable in these prepared files |

The phenotype units remain those in the first table. For example, standardizing genotype dosage leaves the BMI outcome in kg/m².

## 3. Estimate the height/weight moments

A manifest is a small R data frame: one row per trait, with its name, input file and effect scale. The trait names become the variable names in your SEM.

```r
hw_files <- data.frame(
  trait = c("Height", "Weight"),
  path = c("tutorial-data/sumstats/50_irnt.tsv.gz",
           "tutorial-data/sumstats/21002_irnt.tsv.gz"),
  effect_scale = "genotype_standardized"
)

hw_moments <- estimate_moments(
  hw_files,
  scores = "tutorial-data/scores/eur_hm3.scores.npz",
  max_order = 4,
  n_blocks = 200
)

hw_moments
```

The R reference run took about three minutes. It matches alleles, estimates cumulants through order four, and calculates their joint sampling covariance using chromosome-respecting jackknife groups. The returned R object prints:

```text
GCMSEM moments
Traits: Height, Weight
SNPs: 1,074,059
Cumulants: 14 (order 1: 2, order 2: 3, order 3: 4, order 4: 5)
Sampling covariance V: 14 x 14   Jackknife groups: 200

Genetic correlations ($R)
       Height Weight
Height  1.000  0.472
Weight  0.472  1.000
```

The object keeps the quantities you need for analysis:

```r
hw_moments$S                 # SNP-effect covariance matrix
hw_moments$R                 # Genetic correlation matrix
hw_moments$R_SE              # Correlation SEs from the full sampling covariance
head(as.data.frame(hw_moments))  # Individual cumulants and their SEs
V <- vcov(hw_moments)        # Joint covariance of the 14 estimated cumulants
```

The 14 cumulants comprise two means, three second-order terms, four third-order terms and five fourth-order terms. This SEM uses the 12 centered terms of orders 2–4. Order-four terms are excess cumulants: they equal zero for Gaussian sources.

## 4. Specify and fit the height-to-weight model

First put both SNP-effect distributions on a convenient scale:

```r
hw_scaled <- standardize_moments(hw_moments)
hw_scaled$scaling$scale
```

`standardize_moments()` divides each trait's effects by its estimated SNP-effect SD and transforms the entire sampling covariance consistently. These scale constants are held fixed during fitting and SE calculations. In this example the common source has equal direct loadings **on this SNP-effect-SD scale**; choosing another scale changes that model restriction.

Now define the model as an R string:

```r
height_to_weight <- '
  Weight ~ forward*Height

  F =~ 1*Height + 1*Weight
  F ~~ shared_variance*F
  k3(F, F, F) = 0
  k4(F, F, F, F) = 0

  Height ~~ variance_x*Height
  Weight ~~ variance_y*Weight
  k3(Height, Height, Height) = free("skew_x")
  k3(Weight, Weight, Weight) = free("skew_y")
  k4(Height, Height, Height, Height) = free("excess_x")
  k4(Weight, Weight, Weight, Weight) = free("excess_y")
'

hw_model <- gcmsem_model(
  height_to_weight,
  observed = c("Height", "Weight"),
  gaussian_observed = FALSE,
  gaussian_latent = TRUE
)
```

Read the syntax line by line:

- `Weight ~ forward*Height` estimates a directed path named `forward`.
- `F =~ 1*Height + 1*Weight` defines a common source with both direct loadings fixed to one.
- `F ~~ shared_variance*F` estimates its variance freely, including negative values.
- `k3(F, F, F) = 0` and `k4(F, F, F, F) = 0` give this source Gaussian third- and fourth-order cumulants.
- `Height ~~ variance_x*Height` and its weight counterpart estimate trait-specific source variances.
- `free("skew_x")` and `free("excess_x")` name free third- and fourth-order source cumulants. Every variable in a cumulant is written explicitly inside the tuple.

Fit the model with 32 starting points and retain the result in `hw_fit`:

```r
hw_fit <- gcmsem(
  hw_model, hw_scaled,
  orders = 2:4,
  n_starts = 32,
  seed = 20261105,
  feasibility = "check",
  control = list(iter.max = 5000, eval.max = 10000, rel.tol = 1e-11)
)

summary(hw_fit)
```

The reference data produce a negative common-source variance. GCMSEM emits a prominent Heywood warning and retains the estimate and its computable formal SE. The actual summary is:

```text
GCMSEM DWLS fit
Observed traits: Height, Weight
Moments: 12   Free parameters: 8   Jacobian rank: 8/8
Converged: yes   Uncalibrated Q: 1.797
Formal sandwich SEs: available   Inference checks: flagged
*** HEYWOOD CASE: source moment checks failed ***
Negative source variance: F
Estimates and available formal SEs retained.
Flags: source moment feasibility checks failed

Parameter estimates
                               term           label estimate formal SE           se_status      z  p
                   Height ~~ Height      variance_x    1.183    0.1095 formal_only_flagged   10.8 NA
         k3(Height, Height, Height)          skew_x    3.887     2.777 formal_only_flagged    1.4 NA
 k4(Height, Height, Height, Height)        excess_x    735.6     171.2 formal_only_flagged  4.296 NA
                   Weight ~~ Weight      variance_y   0.9935   0.09489 formal_only_flagged  10.47 NA
         k3(Weight, Weight, Weight)          skew_y    4.599     2.446 formal_only_flagged   1.88 NA
 k4(Weight, Weight, Weight, Weight)        excess_y    406.8     185.5 formal_only_flagged  2.193 NA
                             F ~~ F shared_variance  -0.1829    0.0438 formal_only_flagged -4.175 NA
                    Weight ~ Height         forward   0.6549   0.04192 formal_only_flagged  15.62 NA
```

The fit converged with a full-rank parameter Jacobian. Its source-moment checks failed, as the summary states. The `forward` path describes the specified model for genetic effect distributions; a causal interpretation requires the assumed source structure and homogeneous effects. `Q` is an uncalibrated discrepancy. The formal SEs are conditional on this model, the fixed scales and the selected solution; see the [execution report](Execution-report.md) for the higher-order coverage validation.

## 5. Work with the fitted R object

`summary()` returns an inspectable summary object, and the usual R extractors return vectors, matrices and data frames:

```r
hw_summary <- summary(hw_fit)
hw_summary$diagnostics

parameters <- parameter_estimates(hw_fit)
parameters[parameters$parameter == "forward", ]
coef(hw_fit)
vcov(hw_fit)
fitted(hw_fit)
```

To express the forward path in the original phenotype units, multiply its estimate and SE by the outcome-to-predictor scale ratio:

```r
forward <- parameters[parameters$parameter == "forward", ]
scale_ratio <- hw_scaled$scaling$scale["Weight"] /
               hw_scaled$scaling$scale["Height"]

data.frame(
  path = "Height to Weight",
  estimate = forward$estimate * unname(scale_ratio),
  formal_se = forward$se * unname(scale_ratio)
)
```

```text
              path  estimate  formal_se
1 Height to Weight 0.5629883 0.03603639
```

These units are inverse-normal-transformed weight per unit of inverse-normal-transformed height. The [recorded R output](output/causal-reference.txt) includes all six fits and the software session information.

## 6. Repeat for BMI and systolic blood pressure

The [BMI/SBP script](bmi_sbp.R) contains the complete model and all three fits. Sourcing the file defines a function; calling it returns the analysis as a named list of R objects:

```r
source("tutorials/bmi_sbp.R")
bs <- bmi_sbp_example(data_dir = "tutorial-data")

bs$moments
summary(bs$forward)
```

The function builds this manifest internally:

```r
bs_files <- data.frame(
  trait = c("BMI", "SBP"),
  path = c("tutorial-data/sumstats/21001_raw.tsv.gz",
           "tutorial-data/sumstats/4080_raw.tsv.gz"),
  effect_scale = "genotype_standardized"
)
```

Its forward model uses `SBP ~ forward*BMI`, with the same common-source and trait-specific cumulant specifications as above. The actual forward summary is:

```text
GCMSEM DWLS fit
Observed traits: BMI, SBP
Moments: 12   Free parameters: 8   Jacobian rank: 8/8
Converged: yes   Uncalibrated Q: 5.349
Formal sandwich SEs: available   Inference checks: flagged
*** HEYWOOD CASE: source moment checks failed ***
Negative source variance: F
Estimates and available formal SEs retained.
Flags: source moment feasibility checks failed

Parameter estimates
                   term           label estimate formal SE           se_status      z  p
             BMI ~~ BMI      variance_x    1.204    0.1059 formal_only_flagged  11.37 NA
      k3(BMI, BMI, BMI)          skew_x    2.481     2.514 formal_only_flagged 0.9867 NA
 k4(BMI, BMI, BMI, BMI)        excess_x    440.9     270.8 formal_only_flagged  1.628 NA
             SBP ~~ SBP      variance_y    1.233    0.1431 formal_only_flagged  8.615 NA
      k3(SBP, SBP, SBP)          skew_y   -1.759     1.775 formal_only_flagged -0.991 NA
 k4(SBP, SBP, SBP, SBP)        excess_y    414.3     79.36 formal_only_flagged   5.22 NA
                 F ~~ F shared_variance  -0.2044   0.09672 formal_only_flagged -2.113 NA
              SBP ~ BMI         forward    0.319   0.08567 formal_only_flagged  3.724 NA
```

Inspect the path estimates on the original phenotype scale:

```r
bs$paths[, c("model", "from", "to", "original_estimate", "original_formal_se")]
```

```text
       model from  to original_estimate original_formal_se
1    forward  BMI SBP        0.89051544         0.23913901
2    reverse  SBP BMI       -0.02984400         0.01560565
3 reciprocal  BMI SBP        0.77326628         0.14832094
4 reciprocal  SBP BMI       -0.04287837         0.01315318
```

BMI-to-SBP paths are in mmHg per kg/m²; reverse paths are in kg/m² per mmHg. The forward common-source fit retains a negative variance and flags inference. The reverse and reciprocal fits pass the implemented necessary source-moment checks.

## 7. Fit the reverse and reciprocal height/weight models

The [height/weight script](height_weight.R) can reuse the moments already held in R:

```r
source("tutorials/height_weight.R")
hw <- height_weight_example(moments = hw_moments)

summary(hw$reverse)
summary(hw$reciprocal)
hw$paths[, c("model", "from", "to", "original_estimate", "original_formal_se")]
```

```text
       model   from     to original_estimate original_formal_se
1    forward Height Weight         0.5629883         0.03603639
2    reverse Weight Height         1.3957638         0.18554352
3 reciprocal Height Weight         0.4882977         0.01974383
4 reciprocal Weight Height         0.0500362         0.04327605
```

For the reverse model, the script changes the directed path to `Height ~ reverse*Weight` and keeps the same common-source specification. The reciprocal sensitivity model estimates both paths with independent trait-specific sources. It requests `stable = TRUE`, selecting the dynamically stable solution branch. Swapping the source assignment yields an alternative reciprocal parameterization with identical implied moments; this equivalence is examined in the [execution report](Execution-report.md).

You can inspect any fit, including a Heywood case, using the same extractors. Both scripts return `moments`, `standardized`, `forward`, `reverse`, `reciprocal` and the `paths` data frame. Passing `moments = ...` reuses an existing R object; omitting it runs moment estimation first.

## 8. Keep logs or save results when you need them

Logs remain available in the R objects:

```r
cat(tail(hw_moments$log, 8), sep = "\n")
cat(tail(hw_fit$log, 8), sep = "\n")
```

Use `verbose = TRUE` in `estimate_moments()` to print its timestamped estimation log when it completes. In `gcmsem()`, the same argument displays the optimization log as the fit runs. Save an R object with ordinary R functions:

```r
saveRDS(hw_fit, "height-weight-fit.rds")
hw_fit <- readRDS("height-weight-fit.rds")
```

For tab-separated files, export explicitly:

```r
export_gcmsem(hw_moments, "tutorial-results/height-weight-moments")
export_gcmsem(hw_fit, "tutorial-results/height-weight-fit")
```

You can also request export during estimation or fitting with `write.to.tsv = TRUE, out = "path/prefix"`. The function still returns its R object. Export paths must be new unless you explicitly request `overwrite = TRUE`.

## 9. Start from the original Neale Lab downloads

This route uses the public source files directly. You need the score file from step 1 and about 2.2 GB for the four compressed downloads. The following terminal commands save each original file under `raw-gwas/`:

```bash
mkdir -p raw-gwas prepared-gwas

wget -O raw-gwas/height.tsv.bgz \
  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget -O raw-gwas/weight.tsv.bgz \
  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21002_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget -O raw-gwas/bmi.tsv.bgz \
  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz
wget -O raw-gwas/sbp.tsv.bgz \
  https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/4080_raw.gwas.imputed_v3.both_sexes.tsv.bgz
```

The source files identify variants as `chromosome:position:REF:ALT`. They include `beta`, `se`, `n_complete_samples`, `minor_allele`, `minor_AF` and a low-confidence flag. The [Neale Lab export pipeline](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/27.export_ldsc_sumstats.py) defines `ALT` as the effect allele. The adapter below retains that orientation, calculates its GWAS allele frequency from `minor_allele` and `minor_AF`, removes low-confidence records, and matches chromosome, position and alleles to the HapMap3 score metadata:

```bash
python tutorials/prepare_neale_sumstats.py \
  --sumstats raw-gwas/height.tsv.bgz \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --out prepared-gwas/height.tsv.gz

python tutorials/prepare_neale_sumstats.py \
  --sumstats raw-gwas/weight.tsv.bgz \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --out prepared-gwas/weight.tsv.gz

python tutorials/prepare_neale_sumstats.py \
  --sumstats raw-gwas/bmi.tsv.bgz \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --out prepared-gwas/bmi.tsv.gz

python tutorials/prepare_neale_sumstats.py \
  --sumstats raw-gwas/sbp.tsv.bgz \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --out prepared-gwas/sbp.tsv.gz
```

The adapter writes `SNP`, `A1`, `A2`, `BETA`, `SE`, `N`, `EAF`, `CHR` and `BP` columns, plus adjacent `.json` provenance and `.log` files. The height and weight conversions each read 13,791,467 rows and retained 1,074,062 score-tag SNPs in about 92 seconds in the reference run.

Inspect the converted height file in R:

```r
original_head <- read.delim(gzfile("prepared-gwas/height.tsv.gz"), nrows = 3)
print(original_head, digits = 4, row.names = FALSE)
```

```text
        SNP A1 A2      BETA       SE      N    EAF CHR     BP
  rs3131969  G  A  0.002058 0.002467 360388 0.8705   1 754182
 rs12562034  A  G  0.002903 0.002679 360388 0.1060   1 768448
  rs4040617  G  A -0.002119 0.002471 360388 0.1278   1 779322
```

Here `BETA` and `SE` are per allele, so declare **`effect_scale = "allele"`** in the R manifest. `estimate_moments()` performs the genotype-scale conversion and reference-allele harmonization.

```r
hw_original_files <- data.frame(
  trait = c("Height", "Weight"),
  path = c("prepared-gwas/height.tsv.gz", "prepared-gwas/weight.tsv.gz"),
  effect_scale = "allele"
)

hw_original <- height_weight_example(sumstats = hw_original_files)
summary(hw_original$forward)

bs_original_files <- data.frame(
  trait = c("BMI", "SBP"),
  path = c("prepared-gwas/bmi.tsv.gz", "prepared-gwas/sbp.tsv.gz"),
  effect_scale = "allele"
)
bs_original <- bmi_sbp_example(sumstats = bs_original_files)
summary(bs_original$forward)
```

The executed original-source height/weight analysis took about three minutes for estimation and all three SEM fits. It retained 1,074,059 common SNPs after excluding three ambiguous palindromic SNPs, with zero allele or frequency mismatches. Every estimated cumulant, the full sampling covariance `V`, and the correlation matrix `R` were numerically identical to the prepared-file run above. All three fits reproduced the displayed path estimates and SEs. The adapter preserves original GWAS allele frequencies; the returned object's logs and QC metadata record how each SNP was handled.
