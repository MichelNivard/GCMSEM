# Two causal SEM worked examples

This tutorial fits height/weight and BMI/systolic blood pressure using public Neale Lab UK Biobank round-2 summary statistics. It estimates the moments from the released SNP-level summaries and then uses the standalone GCMSEM R package. MCMSEM is not installed or called.

## Inputs

| Pair | Phenotype IDs | Phenotype scale |
|---|---|---|
| Height and weight | `50_irnt`, `21002_irnt` | Separate inverse-normal transformations; not centimetres and kilograms |
| BMI and systolic blood pressure | `21001_raw`, `4080_raw` | Original source units: kg/m² and mmHg |

The distributed BETA and SE columns are on the **genotype-standardized effect scale**, aligned to the score reference. Phenotype units are preserved. `REF_AF` is a reference frequency, not a substitute for GWAS EAF. The [rights page](DATA_RIGHTS.md) identifies the original sources.

Download the `scores` and `causal` assets after publication, or extract the corresponding supplied archives into `tutorial-data/`:

```bash
python tutorials/download_data.py --asset scores --asset causal --out tutorial-data
```

Install both GCMSEM packages as described on the [tutorial home page](README.md). Then run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tutorials/run_causal.py \
    --data-dir tutorial-data \
    --out tutorial-results/causal \
    --starts 32 --seed 20261105
```

The command runs `gcmsem estimate` separately for both pairs with orders 1–4 and up to 200 chromosome-respecting deletion groups. It writes the full joint jackknife covariance, the conditional MAF-stratum results, estimated sampling-error intercepts, QC counts and timestamped logs. It then runs `fit_causal.R` with 32 starts per model. Large outputs remain in the ignored `tutorial-results/` directory.

## Models

For each pair the script fits the same three prespecified models:

1. Forward path with an equal-direct-loading Gaussian common source.
2. Reverse path with the same common-source restriction.
3. Reciprocal paths without the common source, as a sensitivity model.

For example, the forward height/weight specification is equivalent to:

```r
model <- gcmsem_model('
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
', observed = c("Height", "Weight"))
```

Every variable is explicit inside a `k3` or `k4` tuple. Parameter labels are quoted strings in `free(...)`; a numeric right-hand side fixes a source cumulant. Fourth-order parameters are excess cumulants, so Gaussian source cumulants are zero.

The fitted common-source variance is **signed and unconstrained**. A negative value is retained with a prominent Heywood warning. Fixing both loadings at one and estimating the variance avoids implicitly constraining their shared contribution to a nonnegative squared loading.

Before fitting, the script scales each trait by its fixed estimated SNP-effect SD for numerical conditioning and applies the identical linear transformation to both axes of the moment covariance. The equal-loading restriction is imposed on this specified scale. It is a substantive restriction and is not invariant to arbitrary trait rescaling. Output contains both the effect-SD-scale paths and their transformations back to the original phenotype scale, together with the scaling constants.

Only the reciprocal sensitivity requests the dynamically stable branch. The script also evaluates the exactly equivalent source-swapped parameterization and records its implied-moment agreement. Stability selects a branch of that equivalence; it is not evidence that the selected direction is causal.

## Read the outputs

Within each pair's output directory:

| File | Read it for |
|---|---|
| `moments.moments.tsv` | Means, covariances and higher cumulants with jackknife SEs |
| `moments.sampling_cov.tsv` | Joint sampling covariance, already on the estimator scale |
| `moments.metadata.json` and `moments.log` | Source hashes, harmonization counts, overlap correction and deletion details |
| `sem/paths.tsv` | Paths and computable formal sandwich SEs on both reported scales |
| `sem/model_status.tsv` | Convergence, rank, source admissibility, formal-SE and inference status |
| `sem/*.model.txt`, `*.parameters.tsv`, `*.starts.tsv`, `*.log` | Exact model, retained parameters, every start and optimization log |
| `sem/source_swapped.equivalence.tsv` | Reciprocal-model source-assignment equivalence check |

Read `paths.tsv` together with `model_status.tsv`. A finite SE is not an assertion that all source moments are admissible or that inference is validated. No fit is removed simply because it is a higher-order Heywood case. `Q` is an uncalibrated discrepancy; comparing its values does not produce a causal-direction test.

## Executed reference results

The external `results` asset contains the complete executed outputs, not synthetic expected values. A compact execution report is linked from the tutorial home page after validation. Numerical differences can arise from optimizer and linear algebra versions; inspect the logged input hashes and fit statuses before interpreting a difference.

The reference uses rolling LD windows rather than independent block-truncated scores, so LD can cross deletion boundaries. Higher-order finite-reference corrections are unavailable: L1, L3 and L4 are conditional on the supplied reference. These limitations apply to the released examples as they do to the underlying estimator.
