# Executed tutorial checks — 5 November 2026

These are outputs from real summary statistics and actual GCMSEM and LDSC runs. Complete numerical files, input hashes, logs and figures are prepared for the external data release. Zenodo has reserved DOI `10.5281/zenodo.23977299`; the upload is incomplete and publication has not yet been verified.

## Standard LDSC comparison

All **105 of 105** genetic-correlation pairs succeeded with **untouched original bulik/ldsc 1.0.1** and its published standard `eur_w_ld_chr` files for both reference and weighting LD. Statistical defaults were retained. GCMSEM's existing HM3 SNP selection was retained upstream, followed by ordinary LDSC munging and pair-specific allele QC.

| Check | Result |
|---|---:|
| GCMSEM common panel SNPs | 1,012,060 |
| Standard LDSC SNPs after pair-specific QC | 1,012,052–1,057,882 |
| Correlation between the 105 estimates | 0.982855 |
| RMSE between methods | 0.051259 |
| Mean GCMSEM minus LDSC | −0.003636 |
| Largest absolute pair difference | 0.155582 |
| LDSC estimates outside [-1, 1] | 0 |

These summaries are descriptive: pairs share traits and both procedures use the same GWAS. The GCMSEM matrix reproduces the earlier moment-regression result within `1.6e-13`, while its small negative minimum eigenvalue (`−4.73385e-5`) is retained without projection. The standard LDSC comparison uses the ordinary published reference and defaults; the previous custom-reference sensitivity is not substituted for it.

The standard runtime passed **125 upstream tests**. Changing only stock munging's I/O chunk size from 250,000 to 500,000 produced byte-identical decompressed GSCAN outputs. The latter was used to avoid a slow string-membership path in the original pandas environment.

## Two SEM examples

Height/weight used **1,074,059 SNPs**; raw BMI/SBP used **1,074,057 SNPs**. Each run estimated 14 moments with 200 genomic deletions; the R models fitted the 12 centered second-, third- and fourth-order cumulants. Each of the six models used 32 optimizer starts.

All six fits converged with Jacobian rank **8/8** and computable formal sandwich SEs. Three Gaussian-common-source fits had a negative common-source variance. Their point estimates and formal SEs were retained, with explicit Heywood and inference flags.

| Pair and model | Path | Estimate | Formal SE | Source-moment checks |
|---|---|---:|---:|---|
| Height/weight, forward + Gaussian common source | Height → Weight | 0.562988 | 0.036036 | Heywood flagged |
| Height/weight, reverse + Gaussian common source | Weight → Height | 1.395764 | 0.185544 | Heywood flagged |
| BMI/SBP, forward + Gaussian common source | BMI → SBP | 0.890515 | 0.239139 | Heywood flagged |
| BMI/SBP, reverse + Gaussian common source | SBP → BMI | −0.029844 | 0.015606 | Passed local checks |
| Height/weight, reciprocal without common source | Height → Weight | 0.488298 | 0.019744 | Passed local checks |
| Height/weight, reciprocal without common source | Weight → Height | 0.050036 | 0.043276 | Passed local checks |
| BMI/SBP, reciprocal without common source | BMI → SBP | 0.773266 | 0.148321 | Passed local checks |
| BMI/SBP, reciprocal without common source | SBP → BMI | −0.042878 | 0.013153 | Passed local checks |

Height and weight are separately inverse-normal transformed. BMI→SBP uses mmHg per kg/m²; SBP→BMI uses kg/m² per mmHg. The Gaussian common-source equality is defined on the fixed SNP-effect-SD scale documented in the tutorial. Formal SEs condition on that scale. The reciprocal fits select the prespecified stable branch and assume no common source; the exactly source-swapped alternatives reproduce their implied cumulants within relative error `3.2e-16`.

Passing local checks is not a claim of causal identification or calibrated higher-order coverage. Q remains an uncalibrated discrepancy, and source-assignment equivalence remains visible in the outputs.

## Scores and software

All **21,481,240 signed LD-power entries** were independently rebuilt from checksum-verified reference genotypes and agreed within the documented numerical tolerance. Raw reference data are excluded from the release. All 17 exported GWAS gzip files passed complete CRC and row-count checks.

The updated Python suite passed **34 tests**, including truncated-gzip failure logging and independent score oracles. The R suite passed **73 checks**. Detailed build/check results and source hashes are included in the external provenance files.
