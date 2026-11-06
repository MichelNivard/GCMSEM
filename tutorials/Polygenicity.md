# Effective polygenicity in R: standalone analysis and S-LD4M validation

Start with `GCMSEM::estimate_polygenicity()` to analyse your own GWAS and inspect the results in R. Section 3 starts with a relatively precise height estimate, then uses depression to show how to diagnose a less stable result. It needs the R package, its Python engine and GCMSEM scores.

The later validation uses the same 15 GWAS studies as the [genetic-correlation tutorial](Fifteen-trait-correlations.md) to estimate effective polygenicity with each pipeline and compare the results in R. S-LD4M receives the reference scores, annotations and regression weights distributed by its author. GCMSEM keeps its HapMap3 score workflow.

**Correction to the earlier comparison:** the earlier r = 0.9651 used GCMSEM scores and weights in both programs. That is a [matched-input diagnostic](Polygenicity-matched-input.md). This page runs S-LD4M with its own reference inputs and reports its native baseline-category result.

## 1. What the papers establish

O'Connor and colleagues introduced S-LD4M in [*Extreme Polygenicity of Complex Traits Is Explained by Negative Selection* (2019)](https://doi.org/10.1016/j.ajhg.2019.07.003). Their effective number measures the concentration of genetic variance and incorporates LD. Their [supplementary mathematical note](https://github.com/lukejoconnor/SLD4M/blob/dd947104cb3cb86bb8a651bf3b5abf5e5cbf46f3/polygenicity_supp.pdf) relates it to fourth moments and to a variance-weighted LD matrix.

The recent paper is O'Connor and Sella's [*Principled measures and estimates of trait polygenicity* (2026)](https://doi.org/10.1016/j.ajhg.2026.04.009), first posted in 2025. It develops several measures, including harmonic/effective and entropy polygenicity, and estimates them with Fourier Mixture Regression (FMR). The comparison below targets the effective measure. Second and fourth moments alone cannot determine entropy polygenicity or the literal number of nonzero effects. Running FMR would be a separate distributional analysis.

An effective count discounts very small contributions. For independent variants whose expected contributions to genetic variance are $v_i$, it is

```math
M_{e,\mathrm{LE}}=\frac{(\sum_i v_i)^2}{\sum_i v_i^2}.
```

If 1,000 variants contribute equally, this count is 1,000. If a few variants dominate, it is smaller. These are expected variance contributions under the random-effects model; this interpretation matters when connecting the count to moments of signed SNP effects.

## 2. Connect the definition to GCMSEM moments

GCMSEM reports a variance $v$ and an **excess fourth cumulant** $k_4$. The ordinary fourth central moment is $m_4=k_4+3v^2$, and excess kurtosis is $\gamma_2=k_4/v^2$.

For a zero-mean Gaussian variance mixture, the conversion without LD is

```math
M_{e,\mathrm{LE}}=\frac{3Mv^2}{m_4}=\frac{3M}{\gamma_2+3}.
```

Here $M$ counts **source variants** represented by the LD scores. Our release has 1,719,440 sources; the HM3 regression rows are a subset of tags. Substituting the number of tags would change the answer incorrectly.

The factor three follows from the Gaussian-mixture assumption. For example, a point-normal model with a nonzero fraction $\pi$ has $\gamma_2=3/\pi-3$, giving $M_e=M\pi$. For arbitrary non-Gaussian source effects, the same transformation is a Gaussian-mixture-equivalent count, and need not equal the number of causal variants. Negative fourth cumulants or implausible counts remain visible as diagnostic failures of that interpretation.

LD adds another term. Write $R$ for the genotype correlation matrix, $H=\sum_i v_i$, and $\Sigma=\mathrm{diag}(v_i)$. The variance-mixture target is

```math
M_e=\frac{H^2}{\mathrm{tr}\{(\Sigma^{1/2}R\Sigma^{1/2})^2\}}
=\frac{H^2}{\sum_{i,j}R_{ij}^2v_i v_j}.
```

Perfectly correlated variants contribute less independent information. Consequently the inverse-kurtosis count alone is insufficient for a comparison with S-LD4M.

### The calculation used here

Our implementation below is a derived comparison under the zero-mean, independent-source working model. For MAF bin $c$, let $M_c$, $v_c$ and $k_c$ denote its source count, estimated variance and fourth cumulant. Let $a_c$ be its mean marginal genetic variance. We calculate

```math
H=\sum_c M_c v_c,\qquad
D=\sum_c M_c k_c+3\sum_c M_c v_c a_c,\qquad
\widehat M_e=\frac{3H^2}{D}.
```

For the count without LD, replace $a_c$ by $v_c$. This keeps the extra kurtosis caused by different variances across bins. Averaging standardized bin kurtoses would lose that information.

We estimate $a_c$ using the mean second-order score among retained **HM3 tags in bin $c$**, multiplied by the vector of estimated source variances. Applying that mean to every source in the bin assumes that the tags represent the source LD distribution. The released archive does not contain second-order scores at every source SNP, so this is an additional approximation. This approximation belongs to the GCMSEM conversion. The S-LD4M analysis below uses its own complete reference, annotation model, weights and native reporting rule. Different reference SNP universes can therefore change the target as well as the estimate.

## 3. Run a standalone GCMSEM analysis in R

We start with height to learn the workflow with a relatively precise estimate. Height had the lowest relative SE among the 15 traits in the earlier shared-SNP comparison. We rerun it here with the standalone function and its own eligible SNPs. We then apply exactly the same settings to depression to understand a harder case.

Install R package **0.4.0** and Python engine **0.2.0** using the [installation instructions](README.md#install). If updating an existing checkout, reinstall both packages. Keep the Python environment activated when starting R, or set the `python` argument below to its executable.

Download the checksum-verified scores and GWAS panel from Zenodo, from the repository root:

```bash
python tutorials/download_data.py --asset scores --asset panel --out tutorial-data
```

We use three files from those archives:

| File | Contents and purpose |
|---|---|
| `tutorial-data/scores/eur_hm3.scores.npz` | Signed LD-power scores of orders one through four, HM3 tag IDs and alleles, genomic coordinates and deletion blocks, five MAF strata and their source SNP counts. The function uses these scores to separate SNP-effect moments from LD and sampling noise. |
| `tutorial-data/sumstats/50_irnt.tsv.gz` | Prepared GWAS of inverse-normal-transformed height. Columns include `SNP`, `A1`, `A2`, `BETA`, `SE` and `N`. Effects and SEs are already in standardized-genotype units and aligned to the release's reference. |
| `tutorial-data/sumstats/20544_11.tsv.gz` | Prepared depression GWAS with the same columns and standardized-genotype effect scale. This is the second, less precise example. |

### Start with height

Open R in the repository directory:

```r
library(GCMSEM)

height <- estimate_polygenicity(
  sumstats = c(Height = "tutorial-data/sumstats/50_irnt.tsv.gz"),
  scores = "tutorial-data/scores/eur_hm3.scores.npz",
  effect_scale = "genotype_standardized",
  n_blocks = 200L
)
height
```

`Height` names the trait in the returned object. For your own GWAS reported as an effect per additional effect allele, use `effect_scale = "allele"`; GCMSEM converts both BETA and SE using the reference genotype SD. Ensure that your GWAS and scores use the same genome build. SNP IDs alone cannot verify that.

The function checks and harmonizes the GWAS, estimates all four conditional cumulants in each source stratum, recalculates the polygenicity transformation in every genomic deletion, and returns a `gcmsem_polygenicity` object. It retains every eligible tag for this trait and does not trim a tail. With several traits, each is processed independently; adding another trait does not change this trait's SNP set.

The executed height analysis prints:

```text
GCMSEM effective polygenicity
LD-adjusted Gaussian-mixture count; genomic-block jackknife SE
  trait effective_count     se delta_se  n_tags n_sources n_blocks
 Height            5756 725.08   706.58 1074059   1719440      200
Signed estimates retained. Reference-panel uncertainty is not included.
```

Height gives an effective count of **5,756 with SE 725**, or a relative SE of **12.6%**. No diagnostic flags are raised. All five estimated source variances and fourth cumulants are positive, the moment-inequality checks pass, and the estimated denominator is about **5.0 SEs above zero**. These features make height a useful first example of a reasonably precise estimate. They do not establish general confidence-interval coverage or the literal number of causal variants.

`se` is the block-jackknife SE of the complete nonlinear ratio. `delta_se` is a local linear approximation from the joint variance/fourth-cumulant covariance, holding mean LD fixed. For height they are close: 725 and 707. The returned object also contains the reciprocal count and its SE.

### Inspect the height result in R

```r
summary(height, trait = "Height")
print(height$details$Height$stages, row.names = FALSE, digits = 5)
```

`summary()` prints the estimate, all intermediate stages, the four cumulants and SEs for each MAF stratum, the five most influential deletions, and the assumptions. The second command selects just the intermediate table:

```text
            quantity   estimate         se
    genetic_variance 3.2967e-01 3.2980e-02
 fourth_cumulant_sum 5.1893e-05 1.0783e-05
             ld_term 4.7530e-06 7.1275e-07
           numerator 3.2605e-01 6.3638e-02
         denominator 5.6646e-05 1.1363e-05
     effective_count 5.7560e+03 7.2508e+02
   no_ld_denominator 5.2097e-05 1.0809e-05
         no_ld_count 6.2587e+03 8.4687e+02
     excess_kurtosis 8.2119e+02 1.1057e+02
       inverse_count 1.7373e-04 2.1748e-05
```

The numerator is $3H^2$; the denominator is the fourth-cumulant sum plus the LD term defined above. `genetic_variance` is in squared GWAS phenotype units; it is a proportion of phenotypic variance only when the phenotype has unit variance. `no_ld_count` is the conversion before accounting for LD; `inverse_count` is `1/effective_count`. Every stage SE comes from its own deletion estimates and includes covariance between the fitted components.

These are ordinary R objects:

```r
as.data.frame(height)                    # One summary row per trait
coef(height)                             # Named effective counts
height$details$Height$strata              # Source counts, cumulants and SEs
height$details$Height$conditional_covariance # Joint variance/fourth-cumulant covariance
height$details$Height$covariance          # Joint covariance of the intermediate stages
height$details$Height$ld_mean             # Tag-based mean L2, by source stratum
height$moments$Height                     # Full GCMSEM moment object
height$moments$Height$metadata$qc         # Harmonization counts and input hashes
height$moments$Height$metadata$noise_intercept # Estimated sampling-noise calibration
height$log                               # Estimation and diagnostic log
```

Inspect the genomic deletions and draw their sensitivity plot:

```r
blocks <- height$details$Height$jackknife
head(blocks[order(blocks$variance_fraction, decreasing = TRUE),
            c("chrom", "start_bp", "end_bp", "effective_count", "variance_fraction")])
plot(height, trait = "Height")
```

Each row is one leave-block-out calculation. `variance_fraction` is that deletion's contribution to the estimated jackknife variance of the count; it is not a proportion of genetic variance. In the plot, each open point represents a deletion and the red point is the full-data estimate. Height's largest single contribution is **17.1%** of the estimated count variance. Deletion estimates are concentrated near the full-data count.

![Standalone height analysis: denominator and effective count after each genomic deletion.](figures/standalone-height-deletions.png)

### Now try depression: a harder example

Use the same function, reference scores, settings and 200-block maximum:

```r
depression <- estimate_polygenicity(
  sumstats = c(Depression = "tutorial-data/sumstats/20544_11.tsv.gz"),
  scores = "tutorial-data/scores/eur_hm3.scores.npz",
  effect_scale = "genotype_standardized",
  n_blocks = 200L
)
depression
```

The executed analysis prints:

```text
GCMSEM effective polygenicity
LD-adjusted Gaussian-mixture count; genomic-block jackknife SE
      trait effective_count    se delta_se  n_tags n_sources n_blocks
 Depression           20269 35513    29999 1074026   1719440      200
*** POLYGENICITY DIAGNOSTICS: interpret flagged counts and SEs with care ***
Depression: nonpositive_bin_variance; negative_fourth_cumulant_gaussian_mixture_incompatible; source_moment_inequality_violated; denominator_not_separated_from_zero
Signed estimates retained. Reference-panel uncertainty is not included.
```

An R warning reports those same flags when the function returns. The point estimate is **20,269**, but its SE is **35,513**, or **175.2%** of the estimate. The first MAF bin has a negative estimated variance, some fourth cumulants are negative, and the estimated denominator is only **0.69 SEs above zero**. The signed estimates and all 200 deletions remain available for inspection.

Trace where the uncertainty enters:

```r
summary(depression, trait = "Depression")
print(depression$details$Depression$stages, row.names = FALSE, digits = 5)
```

```text
            quantity   estimate         se
    genetic_variance 9.0929e-03 2.1531e-03
 fourth_cumulant_sum 6.5185e-09 1.7642e-08
             ld_term 5.7193e-09 1.3275e-09
           numerator 2.4804e-04 1.1891e-04
         denominator 1.2238e-08 1.7762e-08
     effective_count 2.0269e+04 3.5513e+04
   no_ld_denominator 6.7612e-09 1.7645e-08
         no_ld_count 3.6686e+04 1.8736e+05
     excess_kurtosis 1.3761e+02 3.6752e+02
       inverse_count 4.9337e-05 7.2323e-05
```

The fourth-cumulant sum is especially uncertain: its SE is about 2.7 times its estimate. Adding the LD term gives a denominator whose SE still exceeds its estimate. Dividing by this uncertain denominator produces an unstable count. The full nonlinear jackknife SE is 35,513; the local delta approximation is 29,999. Both indicate weak precision. The function avoids attaching a symmetric confidence interval to this potentially unstable ratio.

```r
plot(depression, trait = "Depression")
```

![Standalone depression analysis: denominator and effective count after each genomic deletion.](figures/standalone-polygenicity-deletions.png)

The two plots use different axis ranges. Compare the **relative** movement around each full-data estimate: depression's deletions span a much larger fraction of its count. All deletion denominators happen to remain positive in this run; that alone does not establish precision. Their jackknife variation still leaves the full-data denominator poorly separated from zero.

### Compare precision, then decide what to report

Build a small comparison table directly from the returned R objects:

```r
comparison <- rbind(as.data.frame(height), as.data.frame(depression))
comparison$relative_se_percent <- 100 * comparison$se / abs(comparison$effective_count)
comparison[, c("trait", "effective_count", "se", "delta_se", "relative_se_percent")]
```

Rounded results:

| Trait | Effective count | Jackknife SE | Delta SE | Relative SE |
|---|---:|---:|---:|---:|
| Height | 5,756.0 | 725.1 | 706.6 | 12.6% |
| Depression | 20,268.7 | 35,513.2 | 29,999.4 | 175.2% |

Height supplies a reasonably precise point estimate for this reference and model. Depression illustrates why a returned count must be read together with its SE, denominator and source-moment diagnostics. Its uncertainty is too large to draw a firm conclusion about how its effective polygenicity compares with height. This is a deliberately chosen teaching contrast; both trait analyses retain their usual inputs and every deletion.

### What every diagnostic check means

These checks describe the `estimate_polygenicity()` workflow. They distinguish **warnings about returned estimates**, **SNP-level input filtering**, and **errors that stop computation**. SEM model-fit and identification diagnostics are described separately in [Working with R results](R-results.md).

The ten polygenicity flags below examine point estimates or their estimated uncertainty. They are not ten formal hypothesis tests. A noisy estimate can fall outside the admissible region even when the population model is valid. GCMSEM retains signed estimates and all computed deletions; it does not constrain them or select favourable deletions to make these warnings disappear.

Write $v_c$, $k_{3c}$ and $k_{4c}$ for variance, third cumulant and fourth cumulant in source bin $c$, and $M_c$ for its source SNP count. The total genetic variance is $H=\sum_c M_cv_c$, and the count is $3H^2/D$, with $D$ defined in Section 2.

- **`nonpositive_bin_variance` — a source-bin variance is zero or negative.** Trigger: any $v_c\leq0$. A negative population variance is impossible; zero is a degenerate boundary where standardized moments cannot be formed. This warning checks the estimate's sign, without a significance test. Inspect the bin's `variance` and `variance_se` in `$strata`. Depression's first bin is only 0.46 SE below zero, so its warning is compatible with sampling noise around a small positive variance. The count remains available.

- **`negative_fourth_cumulant_gaussian_mixture_incompatible` — a source-bin fourth cumulant is negative.** Trigger: any $k_{4c}<0$. For a zero-mean Gaussian variance mixture, $k_4=3\operatorname{Var}(\sigma^2)\geq0$. A negative estimate conflicts with that interpretation at the point-estimate level. Other valid distributions can have negative fourth cumulants. Inspect `fourth_cumulant` and `fourth_cumulant_se`; the first three depression bins are negative but each is less than one SE below zero. This flag does not formally reject the mixture model. An exactly zero fourth cumulant does not trigger it.

- **`source_moment_inequality_violated` — the estimated variance, skewness and kurtosis are mutually inconsistent.** For bins with positive variance, the software checks $k_{4c}\geq k_{3c}^2/v_c-2v_c^2$. Equivalently, ordinary kurtosis must be at least squared skewness plus one. Strong asymmetry requires sufficiently large fourth moments. Failure means those point estimates cannot jointly describe a probability distribution. The check does not use their joint sampling covariance to test statistical significance. Depression fails in the 5–10% and 10–20% frequency bins; its negative-variance first bin is flagged separately and skipped here. Passing is a necessary consistency check, not validation of the complete genetic model.

- **`nonpositive_total_genetic_variance` — the aggregate variance is zero or negative.** Trigger: $H\leq0$. The estimated source contributions do not sum to a positive genetic variance. Examine the bin estimates and their covariance to distinguish weak signal from possible scale or reference problems. Because the numerator squares $H$, the numerical count can still be positive when this warning is present; that positive count does not repair the invalid variance interpretation.

- **`nonpositive_denominator` — the full-data count denominator is zero or negative.** Trigger: $D\leq0$. With a positive numerator, negative $D$ produces a negative effective count, outside its population interpretation; zero $D$ makes division undefined. A finite signed result is retained. Undefined numerical results are represented as `NA` in R, with the denominator and deletion records available for diagnosis. Inspect `$stages` to see the fourth-cumulant and LD contributions.

- **`denominator_not_separated_from_zero` — uncertainty reaches the ratio's singularity.** Trigger: $D-1.96\,SE(D)\leq0$, using the block-jackknife denominator SE. This approximate normal-reference rule checks whether the denominator is sufficiently separated from zero for a stable ratio interpretation. Height has $D/SE(D)=4.99$; depression has 0.69. A failure signals potentially large, asymmetric uncertainty in the count. It is not a calibrated confidence-interval or model-rejection test. Inspect the nonlinear count SE, delta SE and deletion plot together; neither agreement between SE methods nor passing this rule proves coverage.

- **`nonpositive_denominator_in_deletions` — deleting at least one genomic block takes the denominator to or below zero.** Trigger: any leave-block-out $D_{(-b)}\leq0$. The estimate is sensitive enough to encounter the ratio's singularity under genomic deletion. Locate the affected rows in `$jackknife`, using their chromosome and coordinate columns. Negative finite deletion counts are retained; an undefined deletion can make the corresponding uncertainty unavailable. The software does not remove the offending block. Conversely, positive denominators in every deletion do not guarantee precision: depression has all-positive deletions but still fails the preceding uncertainty check.

- **`undefined_ratio_or_uncertainty` — at least one reported stage, deletion value or stage SE is nonfinite.** This checks all ten intermediate quantities, including the no-LD count, excess kurtosis and reciprocal count. A warning can therefore arise from a supplementary quantity even when the primary count is finite. Division by zero or numerical overflow are possible causes. Inspect `$stages` and `$jackknife` for `NA` to identify the affected quantity. Its uncertainty is not silently calculated after dropping problematic deletions; other finite quantities remain available.

- **`count_outside_source_range` — a finite effective count is outside its intended range.** Trigger: $\widehat M_e\leq0$ or $\widehat M_e>\sum_c M_c$. Under the intended effective-count interpretation, the count should be positive and cannot exceed the number of represented source variants. The upper bound uses source SNPs, not HM3 tags. Failure can reflect noisy moment estimates or a poor working-model approximation. The count is retained without clipping. This check does not impose a lower bound of one and does not assess the no-LD count separately.

- **`single_block_dominates_jackknife_variance` — one deletion supplies more than half the estimated count variance.** Trigger: any `variance_fraction > 0.5`. These fractions are contributions to the weighted-pseudovalue jackknife variance; they do not measure the fraction of genetic variance explained by a region. Failure means the reported uncertainty depends heavily on one genomic region. Inspect its coordinates, deletion count and denominator and investigate data quality or regional architecture before interpreting it. A dominant block is not automatically an erroneous block, and remains included. Height's largest contribution is 17.1%; depression's is 15.4% in these standalone runs.

Passing all ten checks means none of these particular problems was detected. It does not verify ancestry or genome-build compatibility, the Gaussian-mixture interpretation, the use of tag LD averages as source LD averages, or confidence-interval coverage. Reference-panel uncertainty is not propagated. The variance and fourth-cumulant sign checks do not require estimates to be significantly positive.

### Additional estimation warnings

The moment engine can also return warnings in `$moments[[trait]]$metadata$warnings` and the log:

| Check | Trigger and meaning | How to respond |
|---|---|---|
| Few genomic deletion groups | Fewer than 30 groups; the jackknife has limited information for uncertainty estimation. | Inspect the available reference blocks and group construction. More requested groups help only when enough distinct blocks exist. |
| Large deletion | A group contains more than 10% of retained tags; the unequal-group jackknife remains approximate. | Inspect block sizes and coordinates. Do not discard a block solely because it is influential. |
| Rank-deficient joint sampling covariance | The number of pooled moments is at least the number of deletion groups. The covariance cannot have full rank. | Avoid treating it as invertible. The polygenicity transformation can still have a computable scalar SE; assess its own diagnostics. |
| Nonpositive pooled SNP-effect variance | A pooled variance is zero or negative, so corresponding genetic correlations are undefined. | Inspect the pooled moment and bin estimates. Signed moment estimates remain available. |

The multivariate moment engine also checks positive semidefiniteness of the estimated genetic correlation matrix. `estimate_polygenicity()` processes each trait separately, so it does not perform a between-trait covariance or SEM identification test. These analyses are documented in [Working with R results](R-results.md).

### Input checks, SNP filtering and errors

Input QC precedes estimation. Removing a malformed or unalignable SNP at this stage is distinct from discarding an inconvenient fitted moment or jackknife result.

| Check | What happens on failure | What to inspect or fix |
|---|---|---|
| R arguments and dependencies | Invalid trait names, effect scales, thresholds or noise-intercept values, missing files, missing `jsonlite`, or an unavailable/incompatible Python engine stop the call. | Follow the error message; supply explicit effect scales and the Python executable containing the current engine. |
| GWAS table structure | Missing/ambiguous required columns, duplicate or empty SNP IDs, malformed tables, or a requested INFO/N filter without its column stop the call. | Repair the table rather than silently selecting among duplicate records. |
| SNP-level QC | Tags absent from the reference, nonfinite BETA/SE, nonpositive SE, invalid supplied EAF, invalid alleles, allele/frequency mismatches and requested INFO/N failures are excluded and counted. | Inspect `$moments[[trait]]$metadata$qc`. A/T and C/G SNPs are dropped by default; keeping them requires EAF to resolve orientation uniquely. |
| Score archive integrity | Missing fields, unsupported schema/orientation/scale, duplicate IDs, invalid coordinates/frequencies/counts, incompatible score dimensions, invalid LD-power values or inconsistent MAF bins stop the call. | Use a compatible, checksum-verified score archive. Structural validation cannot establish biological reference suitability. |
| Sufficient data and supported strata | Too few usable tags for the regressions, an empty source/tag stratum, or a deletion that leaves a stratum without representative tags stops computation. | Inspect coverage and score strata; redesign sparse strata when building scores. Do not hide a failed deletion. |
| Regression identifiability | An empty score column, deficient design rank or a condition number above `max_condition` in a full-data or deletion fit stops computation. | Inspect redundant/unsupported strata and reference coverage. Raising the condition limit alone does not resolve weak identification. |
| Genomic grouping | Blocks must occupy contiguous rows, stay within chromosomes and support at least two deletion groups, with chromosome-respecting grouping and valid fractions. Invalid grouping stops computation. | Correct the score/block metadata and requested group count. |
| Sampling-noise calibration | A supplied or estimated calibration must be finite, symmetric and positive semidefinite; a fixed univariate intercept must be nonnegative. Invalid calibration stops computation. | Inspect scales, SEs and calibration fit. A zero intercept is permitted when explicitly supplied; it asserts no sampling-noise contribution. |
| Numerical and memory limits | Nonfinite fitted cumulants, a failed moment-regression deletion, or estimated memory above `max_memory_mb` stops the analysis. | Inspect the captured error log, input scales and available memory. A failed moment fit does not produce a partial jackknife covariance. |
| Optional export protection | Invalid output arguments or existing managed files without `overwrite=TRUE` stop export. | Choose a new prefix or explicitly permit replacement. Ordinary in-memory use needs no output prefix. |

Warnings and errors can be inspected in R:

```r
depression$estimates$flags                    # Ten polygenicity checks above
depression$details$Depression$warnings        # Individual flag names
depression$details$Depression$strata          # Bin estimates and SEs
depression$details$Depression$stages          # Aggregate quantities and SEs
depression$details$Depression$jackknife       # Every genomic deletion
depression$moments$Depression$metadata$warnings # Moment-engine warnings
depression$moments$Depression$metadata$qc     # SNP filtering counts
depression$log                               # Captured estimation log
```

For a Python-engine failure, the R error also retains its log:

```r
tryCatch(estimate_polygenicity(
  c(Depression = "tutorial-data/sumstats/20544_11.tsv.gz"),
  "tutorial-data/scores/eur_hm3.scores.npz",
  effect_scale = "genotype_standardized"
), gcmsem_engine_error = function(e) {
  cat(conditionMessage(e), "\n")
  cat(e$log, sep = "\n")
})
```

### Arguments and their defaults

Run `?estimate_polygenicity` in R for the installed help page.

| Argument | Default | What it controls |
|---|---|---|
| `sumstats` | Required | A named vector of GWAS paths, a named list of GWAS data frames, or a manifest path/data frame with `trait`, `path`, `effect_scale`. Each GWAS needs `SNP`, `A1`, `A2`, `BETA`, `SE`; supported column aliases follow `estimate_moments()`. |
| `scores` | Required | GCMSEM signed LD-power NPZ, including source counts and disjoint MAF bins. Reference ancestry and genome build should match the study. |
| `effect_scale` | `NULL` | Supply `"allele"` or `"genotype_standardized"` for named files/data frames, optionally named separately for each input trait. With a manifest, omit this argument because the manifest declares the scales. |
| `traits` | `NULL` | Select and order named traits from the input. NULL analyses all independently. |
| `n_blocks` | `200L` | Maximum chromosome-respecting jackknife groups. Uses whole supplied blocks; fewer groups may be available. |
| `noise_intercept` | `NULL` | Estimate sampling-noise calibration from LD regressions. To fix it, supply one nonnegative number or a named number for each selected trait. `1` assumes the supplied GWAS SEs correctly describe sampling variance. |
| `min_info` | `NULL` | Optional minimum imputation INFO; requesting it requires an INFO column. |
| `min_n` | `NULL` | Optional minimum sample size; requesting it requires an N column. |
| `af_tolerance` | `0.1` | Maximum study/reference allele-frequency discrepancy when EAF is available. Must be below 0.5. |
| `allow_strand_flips` | `FALSE` | Allow strand-complement allele matching. |
| `keep_palindromic` | `FALSE` | Keep A/T and C/G SNPs only when EAF uniquely resolves the orientation. |
| `max_memory_mb` | `2048` | Maximum estimated engine working memory per trait, in MiB. |
| `max_condition` | `1e6` | Reject regression designs with a larger condition number. |
| `python` | `Sys.which("python")` | Executable for the environment containing the current GCMSEM Python engine. |
| `write.to.tsv` | `FALSE` | Also export tables, covariance matrices, moments, metadata and logs when TRUE. The primary return value remains the R object. |
| `out` | `NULL` | Output filename prefix, required only when requesting TSV export. |
| `overwrite` | `FALSE` | Permit replacement of the export's existing files. |
| `verbose` | `FALSE` | Display the captured QC and estimation log. The log is always returned in `height$log`. |

`print(height, digits = 5)` controls display precision. `summary(height, trait = NULL, digits = 5)` prints all traits, or a selected subset. `plot(height, trait = "Height")` selects one trait; additional graphics arguments such as `pch` and `cex` can be passed to `plot()`.

### Several traits and saving results

```r
panel <- estimate_polygenicity(
  sumstats = "tutorial-data/manifests/panel.tsv",
  scores = "tutorial-data/scores/eur_hm3.scores.npz",
  traits = c("Height", "Depression")
)
panel
summary(panel, trait = "Height")
```

The manifest supplies effect scales and paths. The function returns per-trait results, without estimating a between-trait covariance matrix for the counts. To retain the entire object:

```r
saveRDS(height, "height-polygenicity.rds")
saveRDS(depression, "depression-polygenicity.rds")
height <- readRDS("height-polygenicity.rds")

# Optional readable exports from an already computed result:
export_gcmsem(height, "results/height-polygenicity")
```

Alternatively set `write.to.tsv = TRUE, out = "results/height-polygenicity"` in the original call. The export includes a summary TSV, per-trait stages, strata, deletions, covariance and LD-mean tables, complete moment bundles, metadata and logs. `saveRDS()` preserves the complete object for a later R session.

The executable [standalone script](standalone_polygenicity.R) runs height first and depression second, saves both R objects and draws both figures. Its [console transcript](output/standalone-polygenicity.txt) contains the full executed output.

### How this connects to the comparison below

The standalone analyses use **1,074,059 height tags** and **1,074,026 depression tags**. The historical 15-trait comparison below used **1,012,000 tags shared across all 15 studies**, including its documented removal of 60 final rows. That subset yields height's earlier 6,202 (SE 849) and depression's 30,796 (SE 108,315). SNP selection matters particularly for the rarest MAF stratum. The package implementation reproduces the historical depression count and SE on those same inputs; its ordinary standalone default uses the full trait-specific eligible set. The comparison remains frozen so the S-LD4M reference correction can be evaluated against unchanged GCMSEM results.

## 4. Understand the two inputs

| Setting | GCMSEM | S-LD4M in this tutorial |
|---|---|---|
| GWAS studies | The 15-study released panel | The same studies, downloaded from their original sources |
| Regression SNPs | 1,012,000 shared HM3 tags in the existing polygenicity analysis | Each GWAS intersected with the author's reference rsIDs |
| Reference | Released GCMSEM scores; 1,719,440 sources | `baselineLD.1kg.l2l4.mat`; 6,227,113 reference rows |
| Annotation model | Five disjoint MAF bins | All 75 supplied baseline-LD annotations |
| Second-moment weights | GCMSEM's fixed inverse-L2 weights | Supplied `EW2` values |
| Fourth-moment weights | GCMSEM's fixed inverse-L2 weights | Supplied `EW4` values |
| Jackknife | Existing 200 genomic deletion groups | The example's 100 equal-row blocks |
| Reported count | GCMSEM moment-to-count conversion above | Native `ma(1)` and `ma_se(1)` |

The [upstream README](https://github.com/lukejoconnor/SLD4M) links the reference archive. Its `example.m` specifies the full `baselineLD.1kg.l2l4.mat` file and the weights above. Appendix A of the 2019 paper describes inverse-L2 weighting for the second moment and inverse-L4 weighting for the fourth moment. Here we pass the author's stored weights directly.

This executes the author's released **1000 Genomes example configuration**. The 2019 paper's main empirical analysis used UK10K. Keep that distinction when comparing numbers with the paper. The archive also includes a smaller `baselineLD.1kg.common.l2l4.mat`; this tutorial uses the full file named by `example.m`.

The earlier GCMSEM setup removed a 60-row tail from 1,012,060 shared HM3 tags. We retain that existing GCMSEM analysis so the change here is the S-LD4M pipeline. S-LD4M now applies its own native tail handling separately to each GWAS; the counts are recorded in `provenance.json`.

## 5. Install and download the data

Start in the GCMSEM checkout and activate its Python environment as described on the [installation page](README.md#install). Install GNU Octave and the HDF5 reader:

```bash
sudo apt-get install octave libopenblas0-pthread
python -m pip install h5py pandas
mkdir -p external original-data/raw tutorial-results

git clone https://github.com/lukejoconnor/SLD4M.git external/SLD4M
git -C external/SLD4M checkout dd947104cb3cb86bb8a651bf3b5abf5e5cbf46f3

wget -O original-data/sld4m-1000g.zip \
  'https://www.dropbox.com/sh/iiyftw01gdpt6un/AACU7AmWK45RxTmDJvRkdKhIa?dl=1'

python tutorials/download_data.py --asset scores --asset panel --out tutorial-data
```

`SLD4M.m` remains unchanged. The execution script checks its SHA-256 against `d8ca72212d1618abfa0a9695d27fa49f45f851cd9641674809fe6a3ed4b5917a` before starting. The reference ZIP used here has SHA-256 `285219fc29a1305732b849a53a3372dfda4188d09a8b9ee89cf8267c990c7bfc`.

The ZIP contains the author's scores, annotations, weights, example script and height example. GCMSEM's [Zenodo release](https://zenodo.org/records/23977299) contains its derived scores and prepared GWAS panel. These numerical inputs remain externally hosted. Allow approximately 40 GB of free disk space for downloads and working files.

For S-LD4M, download the original GWAS files. For example, the height download is:

```bash
wget -O original-data/raw/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz \
  'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
```

Each Neale Lab table contains a genomic variant identifier, effect estimate, SE, `tstat`, sample size and a `low_confidence_variant` flag. `variants.tsv.bgz` supplies the variant-to-rsID mapping. The GSCAN drinks-per-week table contains rsIDs and a `STAT` column that is already chi-square. Its corrected 2019 file is used here.

The [source manifest](sld4m-source-gwas.tsv) lists all 15 GWAS files plus the mapping file, their public download URLs, sizes and checksums. This command downloads the remaining files, verifies every checksum and prepares chi-square inputs:

```bash
python tutorials/prepare_sld4m_gwas.py \
  --download --raw-dir original-data/raw \
  --out tutorial-results/raw-sld4m --workers 2
```

For Neale GWAS, the adapter squares `tstat` and excludes rows marked low confidence by the publisher. It verifies the entire variant-order hash before pairing with the rsID map. For GSCAN it reads `STAT` directly. Both routes require finite nonnegative chi-square statistics and valid rsIDs. There is no significance cutoff and large associations remain included.

Next, retain the original reference order and intersect each trait's rsIDs with that reference:

```bash
python tutorials/prepare_original_sld4m.py \
  --archive original-data/sld4m-1000g.zip \
  --raw-cache tutorial-results/raw-sld4m \
  --out tutorial-results/original-sld4m
```

This creates `ell2.f64`, `ell4.f64` and `annot.f64`, containing the original double-precision matrices in disk-readable form, and `reference_metadata.mat`, containing the original weights, rsIDs and annotation totals. Each `TRAIT.input.mat` contains the matched chi-square values and their original reference row indices. The accompanying JSON records input and matched row counts.

The three full matrices would need over 11 GB in memory. [SLD4MDiskMatrix.m](SLD4MDiskMatrix.m) reads the block requested by the original regression function. It preserves the numerical values and row order. An independent small-array test compares all 14 native outputs against an ordinary in-memory call; every output agreed exactly in the executed check.

## 6. Run both analyses from R

```r
library(GCMSEM)
source("tutorials/original_polygenicity.R")

result <- run_original_polygenicity(
  data_dir = "tutorial-data",
  gcmsem_out = "tutorial-results/gcmsem-polygenicity",
  native_dir = "tutorial-results/original-sld4m",
  out = "tutorial-results/original-polygenicity",
  sld4m_dir = "external/SLD4M",
  jobs = 1L
)
result
```

GCMSEM estimates each trait's moments through order four. S-LD4M then fits each original GWAS independently with its supplied annotations and weights. The R return value contains `comparison`, `annotations`, `moments`, `agreement` and `provenance`. The small on-disk tables also allow the analysis to be reopened later.

For a machine with more memory and cores, increase `jobs`. The executed S-LD4M run used Octave 8.4.0 with up to six concurrent traits and one BLAS thread each. Each trait has a separate log, such as `tutorial-results/original-sld4m/Height.log`.

The native call inside the wrapper is:

```matlab
[ma,ma_se,h,h_se,ma_enrich,ma_enrich_se,h_enrich,h_enrich_se,ma_jk,h_jk] = ...
    SLD4M(chisq,ell2,ell4,idx_ref,annot,EW2(idx_ref),EW4(idx_ref),100);
```

This preserves the default reporting annotations, baseline category, point-estimate convention and jackknife-SE calculation. We read `ma(1)` and `ma_se(1)` directly. The full signed annotation estimates are retained in `result$annotations`. The `h` output is left in its native chi-square input scale; it is not labelled a heritability estimate in the comparison.

### Executed output

```text
GCMSEM versus S-LD4M with the author's reference inputs
15 traits; S-LD4M: 75 annotations, 100 jackknife blocks
Native baseline-category estimates; separate SNP selection and weights
         trait   gcmsem  gcmsem_se     sld4m sld4m_se
 DrinksPerWeek 15910.96   4375.655  6454.315 3917.312
      VitaminD   273.40     99.255    50.712   26.578
           SBP 11971.02   2022.819  8737.174  630.430
      Insomnia 26836.94   6143.741 14450.802 1434.939
    Depression 30796.43 108314.673  4735.252 1497.308
           LDL   121.22    489.674    86.001   43.721
          FEV1 11874.57   1958.386  7086.678 1625.868
   HeartAttack  4615.55   6364.806  2125.719  582.328
        Angina  3017.08   8022.627  1999.420  632.905
        Height  6202.17    848.620  3612.232  229.116
        Weight  9002.62   3650.955  6578.993 1191.541
           BMI  9478.37   6276.382  7030.477 2026.151
           DBP 11779.27   2057.020  9855.071  678.169
           FVC 10620.47   1590.254  6439.964 1153.949
         Waist  8857.68   6477.507  6913.366 1906.715

Pearson r: 0.66682
RMSE: 8190.36 effective variants
Traits with diagnostic flags: 12 of 15
Inspect $comparison$flags and $annotations; signed estimates are retained.
```

## 7. Inspect the R objects and draw the scatterplot

```r
result$comparison
subset(result$annotations, trait == "Height")
result$moments$Height
result$comparison[, c("trait", "gcmsem_tags", "sld4m_matched_snps", "flags")]
saveRDS(result, "original-pipeline-polygenicity.rds")

plot(result, uncertainty = FALSE)
png("original-polygenicity-scatterplot.png", width = 1600, height = 1200, res = 160)
plot(result, uncertainty = FALSE)
dev.off()
```

![All 15 traits: GCMSEM versus S-LD4M with the author's complete reference inputs.](figures/original-polygenicity-scatterplot.png)

Across all 15 traits, Pearson **r = 0.6668** and RMSE is **8,190.4 effective variants**. The earlier matched-input r = 0.9651 therefore substantially overstated the agreement obtained with independently configured reference pipelines.

For height, S-LD4M estimates 3,612 (SE 229), versus GCMSEM's 6,202 (SE 849). For drinks per week, the estimates are 6,454 (SE 3,917) and 15,911 (SE 4,376). Depression contributes a large discrepancy: 4,735 (SE 1,497) versus 30,796 (SE 108,315). Its very large GCMSEM uncertainty limits interpretation of that point. All 15 traits remain in the reported correlation.

Show the full one-SE uncertainty bars as well:

```r
plot(result, uncertainty = TRUE)
```

![The same 15 traits with full one-standard-error bars.](figures/original-polygenicity-uncertainty.png)

Orange triangles identify diagnostic flags. Negative estimates and intervals crossing zero remain visible. These ratios can be unstable when estimated fourth-moment denominators are near zero; symmetric SE bars are descriptive rather than calibrated confidence intervals. The two methods share GWAS studies, so their errors are correlated.

Their references and SNP universes differ. The scatterplot measures agreement between independently configured pipelines; its correlation alone cannot establish accuracy or equality of their population targets. GCMSEM additionally uses a Gaussian-mixture interpretation and HM3-tag LD averages in its moment-to-count transformation.

## 8. Reopen results and inspect provenance

```r
result <- read_original_polygenicity(
  "tutorial-results/original-polygenicity",
  "tutorial-results/gcmsem-polygenicity"
)
result$agreement
result$provenance$reference
result$provenance$traits
```

The [execution record](../docs/original-sld4m-validation.json) records software/reference hashes, every trait's SNP count and native output checksum. The earlier [matched-input diagnostic and simulation](Polygenicity-matched-input.md) remain available to examine the moment conversion under shared inputs and known LD. They answer a different question from the author-reference comparison above.
