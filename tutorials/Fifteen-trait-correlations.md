# Reproduce a 15 × 15 correlation matrix with LDSC and GCMSEM

This tutorial runs all 105 distinct pairs using the same ordered common set of real GWAS SNPs. It compares the current GCMSEM package with the actual pinned LDSC regression implementation, with an independent LDSC command-line check for height/weight.

## Panel

| Order | Trait | Public source phenotype |
|---|---|---|
| 1 | Drinks per week | GSCAN 2019 public release, excluding 23andMe |
| 2 | Vitamin D | Neale `30890_raw` |
| 3 | Systolic blood pressure | Neale `4080_irnt` |
| 4 | Insomnia | Neale `1200` |
| 5 | Diagnosed depression | Neale `20544_11` |
| 6 | Direct LDL | Neale `30780_raw` |
| 7 | FEV1 | Neale `3063_irnt` |
| 8 | Heart attack | Neale `6150_1` |
| 9 | Angina | Neale `6150_2` |
| 10 | Height | Neale `50_irnt` |
| 11 | Weight | Neale `21002_irnt` |
| 12 | BMI | Neale `21001_irnt` |
| 13 | Diastolic blood pressure | Neale `4079_irnt` |
| 14 | FVC | Neale `3062_irnt` |
| 15 | Waist circumference | Neale `48_irnt` |

The panel contains **1,012,060 common SNPs**, while the LD source universe remains **1,719,440 variants**. The causal BMI/SBP tutorial uses raw phenotypes; this matrix uses their IRNT versions to match the previous panel. The two analyses must not be described as the same phenotype definition.

## Run GCMSEM

The existing GCMSEM HM3 selection is retained. GCMSEM uses the exact prealigned 15-trait common panel. Stock LDSC starts from the same released per-trait HM3-selected summary statistics and then performs its ordinary munging and pair-specific QC. It is not forced to use GCMSEM's complete-case intersection.

After installing GCMSEM, download `scores` and `panel`:

```bash
python tutorials/download_data.py --asset scores --asset panel --out tutorial-data

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tutorials/run_correlations.py \
    --panel tutorial-data/panel.npz \
    --scores tutorial-data/scores/eur_hm3.scores.npz \
    --out tutorial-results/gcmsem --stage gcmsem
```

This estimates 15 means and 120 symmetric second-order moments, refitting nuisance parameters in all 200 genomic deletions. The current package reproduces the earlier GCMSEM matrix within 1.6e-13. Its smallest correlation-matrix eigenvalue, approximately -4.73e-5, is retained rather than projected away. No higher-order tensor is requested for this comparison.

## Run absolutely standard LDSC

The primary comparator is **original, unmodified bulik/ldsc 1.0.1**, commit `2fdeeb3b44379408794154993dbd6101b8946b7e`. It runs in its own Python 2.7 environment. Both command-line scripts and all statistical code remain unchanged. The restored runtime passed 125 upstream tests. See [the environment guide](Standard-LDSC-environment.md) for installation, then verify the checkout:

```bash
python tutorials/setup_standard_ldsc.py \
  --out external/standard-ldsc --python2 /path/to/ldsc-environment/bin/python2.7

python tutorials/download_standard_ldsc.py --out external/standard-reference
```

The reference downloader fetches the **published `eur_w_ld_chr` files** without recomputing them. It uses existing DOI archives because the original Broad download links now return 404. All 48 extracted files were independently confirmed byte-identical to the Colorado GenomicSEM practical mirror. The standard files contain 1,290,028 score rows and a default M_5_50 count of 1,173,569. The same published prefix is used for `--ref-ld-chr` and `--w-ld-chr`, exactly as in the [official LDSC genetic-correlation tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

Prepare the released summary statistics and run all pairs:

```bash
python tutorials/prepare_standard_ldsc_inputs.py \
  --sumstats tutorial-data/sumstats --out tutorial-data/standard-inputs

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tutorials/run_standard_ldsc.py \
    --inputs tutorial-data/standard-inputs \
    --hm3 external/standard-reference/w_hm3.snplist \
    --ref-prefix external/standard-reference/eur_w_ld_chr \
    --ldsc-dir external/standard-ldsc \
    --python2 /path/to/ldsc-environment/bin/python2.7 \
    --out tutorial-results/standard_ldsc \
    --gcmsem-results tutorial-results/gcmsem
```

The preparation step writes `SNP A1 A2 BETA P N`. It derives two-sided normal P values from the released BETA/SE ratio; the shared genotype-scale transformation cancels in that ratio. It does not pass reference AF as study frequency or invent INFO. Extremely small P values that underflow are left for stock LDSC's normal invalid-P handling. The original GWAS phenotype definitions and upstream HM3 selection are documented in the data provenance.

The runner invokes stock `munge_sumstats.py` with `--merge-alleles` and explicit column mappings, then stock `ldsc.py --rg` with the published reference prefix for both reference LD and weighting LD. It does not override M_5_50, constrain intercepts, disable allele checking, or impose a custom two-step setting. The munging I/O chunk size is 500,000 rows to avoid a slow path in the old pandas runtime; this changes batching, not the statistical procedure. Exact subprocess commands and logs are saved.

## Inspect the results

The primary output directory contains the standard LDSC correlation and SE matrices, all 105 pair records, unmodified command-line logs, munging counts and a comparison against GCMSEM. Each method reports its own SNP counts and SEs. Input hashes, reference-file hashes, the untouched software revision and dependency versions bind a run to its provenance.

GCMSEM targets correlation of centered SNP effects pooled over its MAF-stratum distribution. Conventional LDSC uses its standard effect assumptions, reference, weighting and jackknife procedure. This is a comparison between the two procedures as used, not a claim that their estimands and intermediate fits are algebraically identical. Binary-trait results are not silently converted into liability-scale heritabilities.

The external `results` asset contains the fresh executed output. Any older matched-reference runs using recomputed custom LD weights are sensitivity analyses only; they are not the primary benchmark and are excluded from the default result archive. The scripts for such a sensitivity remain available for method development, with explicit labels.

This second-order comparison does not validate higher-order interval coverage or identify causal direction.
