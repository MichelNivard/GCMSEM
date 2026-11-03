# GWAS summary statistics to genomic cumulants

`gcmsem estimate` accepts any positive number of traits. The GWAS manifest has
three required columns: `trait`, `path`, and `effect_scale`. Paths are relative
to the manifest. `effect_scale` must explicitly be `allele` or
`genotype_standardized`; an optional `phenotype_unit` is retained in provenance.

```text
trait  path          effect_scale
BMI    bmi.tsv.gz    allele
SBP    sbp.tsv.gz    allele
LDL    ldl.tsv.gz    allele
```

Use tabs, commas, or whitespace as field delimiters; use quoted CSV/TSV if paths
contain spaces. Trait identifiers start with a letter and contain letters,
digits, periods or underscores. Input GWAS tables require `SNP A1 A2 BETA SE`.
Accepted aliases are `rsid/rs_id/variant_id`, `effect_allele/ea/alt`,
`other_allele/nea/ref`, `b/effect`, and `standard_error/stderr`, respectively.
Names are case-insensitive. Ambiguous or repeated columns cause an error.
`EAF`, `effect_allele_frequency`, `a1freq`, or `a1_freq` optionally supplies the
frequency of the reported effect allele. `INFO` and `N` are required when the
corresponding filter is requested. The software does not infer log-odds,
liability, phenotype units, ancestry or genome build.

```bash
gcmsem estimate --sumstats traits.tsv --scores european.scores.npz \
  --out results/moments --blocks 200 --max-order 4
```

The harmonizer retains a common complete-case set of tag SNPs across all
requested traits. It logs per-trait exclusion counts, sign flips, strand flips,
missingness, frequencies and retained counts. Duplicated IDs fail rather than
being arbitrarily deduplicated. Nonfinite coefficients and nonpositive SEs are
excluded. A/T and C/G variants are excluded by default. `--keep-palindromic`
requires EAF that uniquely resolves orientation; ambiguous frequencies remain
excluded. Nonpalindromic strand complements require `--allow-strand-flips`.
Per-allele effects and SEs are multiplied by the reference genotype SD
`sqrt(2*p*(1-p))`. **No significance pruning is performed.** Filtering tag
rows never changes the reference/source variant counts used for the target
distribution. All inputs must be on a compatible genome build and phenotype
scale; the software cannot establish these from SNP identifiers alone.

## Statistical operation

For reference correlation matrix R, the score of degree q for tag j and source
annotation c is `Lq[j,c] = sum(k in c) R[j,k]^q`. Odd powers retain the sign.
Different source SNPs are assumed independent conditional on disjoint
annotations; their effect cumulants are constant within each annotation. This
is a working genetic architecture, not a claim about all biological effects.

1. Regress each trait's effect estimates on L1 to estimate annotation-specific
   means. Subtract their LD-mixed prediction to obtain centered marginal
   estimates x.
2. For every trait pair a,b, regress `x[a]*x[b]` jointly on L2 and
   `SE[a]*SE[b]`. The L2 slopes estimate source-effect covariance within each
   annotation. The last coefficient estimates the sampling-error intercept
   `I[a,b]`. This is the LDSC-style component that accounts for overlap;
   neither sample overlap counts nor genetic correlation are substituted for
   this intercept. A labeled `--noise-intercept` matrix can instead supply
   externally justified calibration, held fixed in every deletion.
3. Construct third-order Gaussian noise corrections. For example,
   `H[a,a,b]=x[a]^2*x[b]-Omega[a,a]*x[b]-2*Omega[a,b]*x[a]`, where
   `Omega[a,b]=SE[a]*SE[b]*I[a,b]`. Regress each corrected response on L3.
4. For order four, subtract every single-noise-pair contribution, add the three
   noise covariance pairings, and subtract the three **fitted signal covariance**
   pairings. Regress the result on L4. For one trait the response is
   `x^4 - 6*Omega*x^2 + 3*Omega^2 - 3*C_signal^2`.
5. Pool annotations by converting conditional cumulants to raw moments,
   weighting by reference/source counts, then converting back to cumulants.
   Averaging conditional cumulants directly would be incorrect when annotation
   means or variances differ. Annotation-specific bundles are also exported.

Each regression uses the fixed redundancy weight
`1/max(sum_c L2[j,c],1)`. Predictor columns are scaled, and rank/condition checks
apply to the full data and every deletion. This is not full-covariance WLS on
the individual SNP regressions. The R SEM stage subsequently uses DWLS on the
estimated cumulant vector, with its joint sampling covariance for sandwich SEs.

The estimated sampling-intercept matrix must be positive semidefinite. An
invalid matrix stops the run, including when it occurs only in a jackknife
deletion. There is **no silent nearest-PSD projection**. This failure can reflect
weak separation of LD-dependent signal and error, an inadequate source
architecture, finite sampling variation, or incompatible inputs. Supplying a
known matrix is appropriate only when it is independently justified; it is not
a switch to manufacture a successful result.

## Uncertainty and variable LD-block lengths

The estimator aggregates consecutive whole reference LD blocks into groups
with approximately equal tag-row counts, without crossing chromosomes. It
never divides an LD block to force equal sizes. Each deletion refits means,
noise intercepts, covariances and all selected higher cumulants. The source
universe and mixture proportions remain fixed.

For equal groups the joint covariance is `(B-1)/B` times the centered deletion
cross-product matrix. For unequal groups with tag fraction f_b, the weighted
delete-group jackknife uses

```text
h_b     = 1/f_b
tau_b   = h_b*theta - (h_b-1)*theta_without_b
theta_J = sum_b f_b*tau_b
V       = (1/B) sum_b (tau_b-theta_J)(tau_b-theta_J)'/(h_b-1)
```

This is the [Busing/Patterson weighted jackknife formulation](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf).
It remains an approximation for nonlinear regressions with heterogeneous LD.
The log reports group counts, fraction range, effective group count and warnings
for dominant groups or fewer than 30 groups. Nominal higher-order coverage must
be calibrated by simulations relevant to the intended dataset. Reference-LD
sampling uncertainty is not included.

The full joint V may be rank deficient when its dimension exceeds B-1. This is
reported, not treated as a failure: robust DWLS does not invert V. Rank is
calculated from its B-by-M factor rather than a cubic decomposition of the
potentially much larger M-by-M covariance.

## Outputs and scalability

The prefix produces:

| Suffix | Content |
|---|---|
| `.moments.tsv` | Canonical tuple key, order, estimate, jackknife SE |
| `.sampling_cov.tsv` | Full labeled joint covariance of all moment estimates |
| `.jackknife.tsv` | Labeled deletion estimates |
| `.covariance.tsv` | N-by-N pooled source-effect covariance |
| `.correlation.tsv` | N-by-N pooled source-effect correlation |
| `.metadata.json` | Schema, units, traits, assumptions, QC, inputs and SHA-256 hashes |
| `.log` | Timestamped command, versions, stage progress, warnings and outputs |

When there is more than one annotation, `.stratum-1`, `.stratum-2`, etc. prefixes
contain the same files for annotation-specific moments. Those bundles can be
read directly by `GCMSEM::read_moments()`. Metadata is installed last as the
completion marker; overwriting invalidates the old marker before replacing
files. Existing results are protected unless `--overwrite` is passed.

Fourth-order values are excess cumulants, not ordinary kurtosis. Off-diagonal
third and fourth cumulants are symmetric but retain multiplicity through their
keys: `A|A|B` differs from `A|B|B`. Means are exported for diagnostics even when
the SEM fit selects only orders 2–4. Nonpositive effect variances leave their
correlations undefined; an indefinite estimated genetic covariance is reported
without silently projecting the estimates.

The number of unique moments through degree q is `choose(N+q,q)-1`. For 15
traits this is 3,875 through degree four, versus 135 through degree two.
Moment responses are evaluated one at a time; the code does not allocate an
SNP-by-fourth-order-tensor array. The joint V still costs O(M^2) memory and
all-stage jackknifing can be expensive. `--max-order 2` is the fast matrix
route. `--moments selected.txt` retains selected higher tuples plus all means,
covariances and required lower-order subtuples. The file contains one
`A|B|C|D` key per line. `--max-moments` and `--max-memory-mb` reject unexpectedly
large requests before estimation.

Python callers can use `estimate_moments(beta, se, scores, traits=..., ...)`
directly. Arrays are already harmonized, genotype-standardized and complete;
this low-level API cannot establish allele correctness from unnamed arrays.
It returns pooled estimates, V, all deletions, conditional arrays, intercepts,
canonical tuple indices and metadata. Use `write_results(result, prefix)` for
the R interchange format.

Pooled annotation mixtures can induce dependence between otherwise independent
SEM sources. Conditional-stratum fits are therefore useful sensitivity checks.
Robust SEs address correlated estimation error; they do not resolve this
architectural issue, source-label ambiguity, weak identification or causal
interpretation.

## Reproducible dimensionality smoke

`validation/scalability.py` exercises the full 15-trait, order-four route on
120 independent synthetic tag SNPs. It checks all 3,875 cumulants and the full
3,875-by-3,875 sampling covariance, including the expected rank of five from
six deletion groups. This deliberately small fixture tests execution and
dimensions; it does not establish statistical calibration or genome-wide
runtime. The approximately 230 MiB memory figure is the software's conservative
guard calculation, not measured process RSS. The recorded runtime and software
versions are in `validation/results/scalability.json`.

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=python python validation/scalability.py
```
