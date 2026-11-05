# Executed tutorial checks — 5 November 2026

These are outputs from real summary statistics and actual GCMSEM and LDSC runs. The original numerical files, input hashes and logs are public at [Zenodo, DOI 10.5281/zenodo.23977299](https://zenodo.org/records/23977299). All eight deposited files are complete, including the four checksum-verified archives. The updated R 0.3.0 tutorials recompute moment estimates and add model-fit diagnostics.

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

The refitted models have the following Q diagnostics. Q SEs use the observed residuals as a plug-in estimate of misfit; the null reference projects the full moment covariance off the fitted model Jacobian.

| Pair / model | Q | Approximate Q SE | Approximate null 95th percentile | Approximate tail probability | Source checks |
|---|---:|---:|---:|---:|---|
| Height/weight forward | 1.797 | 1.257 | 0.903 | 0.00331 | Heywood flagged |
| Height/weight reverse | 5.456 | 4.881 | 4.230 | 0.02354 | Heywood flagged |
| Height/weight reciprocal | 4.286 | 2.013 | 1.018 | 0.00000710 | Passed necessary checks |
| BMI/SBP forward | 5.349 | 4.848 | 6.199 | 0.07704 | Heywood flagged |
| BMI/SBP reverse | 16.389 | 10.801 | 10.545 | 0.00742 | Passed necessary checks |
| BMI/SBP reciprocal | 4.193 | 3.037 | 4.810 | 0.07761 | Passed necessary checks |

All three height/weight models and the reverse BMI/SBP model exceed their respective approximate reference thresholds. The tutorial therefore omits effect interpretation for those models. It illustrates conditional interpretation with the reciprocal BMI/SBP model: BMI → SBP = 0.773266 (formal SE 0.148321) mmHg per kg/m², and SBP → BMI = −0.042878 (formal SE 0.013153) kg/m² per mmHg. These estimates assume independent trait-specific sources and the selected stable branch.

The approximate reference is diagnostic; calibrated higher-order GWAS fit inference is not established. The [Q derivation](R-results.md#how-the-approximate-se-is-calculated) explains its fixed-V, local-Gaussian assumptions. The external archives retain every original parameter estimate, including Heywood and poorly fitting models, for audit and simulation aggregation. Every refitted coefficient reproduced its archived value within `1e-10`; adding Q diagnostics does not alter the optimizer or discard estimates.

The equal common-source loadings are defined on the fixed SNP-effect-SD scale. The reciprocal model omits the common source and requests a stable branch. Source-swapped reciprocal parameterizations retain identical implied cumulants within relative error `3.2e-16`. Fit and local rank checks do not resolve this causal-model equivalence.

## Scores and software

All **21,481,240 signed LD-power entries** were independently rebuilt from checksum-verified reference genotypes and agreed within the documented numerical tolerance. Raw reference data are excluded from the release. All 17 exported GWAS gzip files passed complete CRC and row-count checks.

The Python suite passed **49 tests**, including the original score and estimator checks and 15 tests of the Neale Lab tutorial adapter. The R 0.3.0 checks cover the original 73 fitting checks, 37 presentation/export checks, 25 R-to-Python integration checks, moment-transformation/export regressions, and 14 new Q-uncertainty checks. The latter include independent Gaussian simulation with 80,000 refits, fixed-model cases, singular covariance and saturated models. The bridge checks include an independent three-trait covariance calculation and a real Python virtual environment.

`R CMD check` of the built GCMSEM 0.3.0 source package completed with **0 errors, 0 warnings and 0 notes**. The Python engine was available during that check, so the R-to-Python tests executed. The [current verification record](../validation/results/tutorial_q_checks.json) records the checks and tested source hashes.

The R-native height/weight estimation was rerun from the prepared GWAS files: it retained 1,074,059 SNPs and reproduced every archived cumulant and the full sampling covariance within `1e-12`. All six fits from the new R tutorial scripts reproduced the archived coefficients and parameter covariance matrices within `1e-12`. The [causal transcript](output/causal-reference.txt) records the readable R output.

Starting from the original Neale Lab height and weight downloads also retained 1,074,059 common tags. After the tutorial adapter and R estimation, the moment table, sampling covariance and correlation matrix were exactly identical to the prepared-file route. All three height/weight refits reproduced the reference results.

The correlation tutorial now calls `estimate_moments()` on all 15 GWAS tables. This run took 858.7 seconds, retained 1,012,060 common SNPs and used 200 genomic deletions. Its correlation matrix was exactly equal to the archived matrix. The potentially affected palindromic tags were already absent from the 15-trait intersection; the earlier suggestion that the normal route necessarily uses three fewer common SNPs was incorrect. The [console transcript](output/correlation-comparison.txt) and [scatterplot](figures/correlation-scatterplot.png) now come from this executed estimation workflow.

All six model diagrams are generated directly from the corresponding R model objects by [plot_models.R](plot_models.R). The updated [causal transcript](output/causal-reference.txt) reports fit diagnostics first and the conditional reciprocal BMI/SBP paths.

The external release preserves the original numerical archives and provenance. New software, diagrams, transcripts and tests are versioned with R 0.3.0 in this repository.
