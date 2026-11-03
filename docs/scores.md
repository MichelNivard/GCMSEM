# Building signed LD-power scores

`gcmsem score` constructs the reference quantities needed to recover joint SNP-effect cumulants from marginal GWAS estimates. It accepts SNP-major PLINK BED/BIM/FAM references, genotype arrays, or signed LD matrices. It processes one explicitly defined LD block at a time; it does not construct a dense genome-wide LD matrix.

For tag SNP j, source stratum c, and order q = 1, 2, 3, 4:

```math
L^{(q)}_{jc}=\sum_{k\in c} R_{jk}^{q}.
```

`R` is the signed correlation of genotype-standardized, minor-allele-coded genotypes. In particular, third powers preserve negative signs. Squared or absolute-valued LD input cannot recover these scores and is not a valid substitute. A supplied matrix cannot be automatically recognized as unsigned: its provenance remains the user's responsibility.

The **tag set** contains the SNPs at which GWAS products will be regressed. The **source universe** contains all reference SNPs whose effects may contribute to those marginal associations. Restricting tags to HapMap3 does not remove other QC-passed reference SNPs from source scores or counts. Source strata are disjoint MAF bins. Keep zero-effect SNPs; neither score creation nor moment estimation performs significance filtering.

## Command-line input

Create a tab-separated manifest:

```text
block   path              format  n_reference  block_map
1       eur_chr1          bed                  chr1_blocks.tsv
2       eur_chr2          bed                  chr2_blocks.tsv
```

The displayed spaces stand for tabs. Paths are resolved relative to the manifest. For `bed` rows, `path` is a PLINK prefix, and `block_map` supplies genomic intervals. The manifest `block` field is unused for BED rows; interval IDs provide their block IDs. Each reference needs its `.bed`, `.bim`, and `.fam` files. SNP-major BED encoding is required. Genotypes are decoded as A1 dosage using the alleles in BIM. Sex chromosomes requiring non-diploid dosage handling should be preprocessed explicitly; this reader treats all dosages as diploid.

A block map is a tab-separated file with **1-based inclusive** coordinates:

```text
chrom   start       end         block
1       1           1000000     1:001
1       1000001     2500000     1:002
```

Use defensible LD blocks, with globally unique IDs. Every BIM variant must be assigned to exactly one interval. Missing or overlapping assignments produce errors. Arbitrarily chopping a chromosome into chunks does not establish that the resulting blocks are independent. Cross-block LD is assumed zero; a block boundary through a strong LD region can bias scores and compromise subsequent inference.

Run:

```bash
gcmsem score \
  --manifest reference_blocks.tsv \
  --tags hm3.snplist \
  --maf-edges 0,0.05,0.1,0.2,0.3,0.5 \
  --min-maf 0.01 \
  --finite-reference \
  --out eur_hm3.scores.npz
```

`hm3.snplist` contains one SNP identifier per line, without a header. Omit `--tags` to use every QC-passed reference SNP as a tag. An empty source bin is an error: merge bins or provide a larger reference, rather than assigning a moment to an empty population. Bin edges are lower-inclusive and upper-exclusive, except that the final bin includes MAF 0.5. Frequencies within 1e-12 of an edge are snapped to that edge for bin assignment to avoid allele-recoding roundoff differences.

`--max-distance-bp` can additionally truncate within-block LD contributions. This option changes the assumed LD model; it does not compute LD across the supplied block boundaries. `--max-block-variants` defaults to 10000 and prevents an unexpectedly large block from exhausting memory. Reference matrices remain block-dense; score intermediates are processed in chunks of at most 512 tag rows. Approximate peak memory is at least 8 m² bytes for a block of m variants, plus genotype data and work arrays. `--skip-psd-check` skips the expensive minimum-eigenvalue test for trusted matrices; the choice is recorded.

## Matrix and genotype inputs

For NPZ input, use manifest rows such as:

```text
block   path          format  n_reference
1:001   block1.npz    npz     1000
1:002   block2.npz    npz     1000
```

Each pickle-free archive contains:

| Field | Required value |
|---|---|
| `ld` **or** `genotype` | Signed m × m LD correlation matrix, or samples × m A1 dosage matrix |
| `snp` | m unique SNP IDs |
| `a1`, `a2` | m biallelic A/C/G/T allele labels |
| `af` | m A1 frequencies; required for `ld`, computed from `genotype` |
| `chrom` | m chromosome labels, all on the same chromosome |
| `bp` | m positive integer base-pair positions |
| `n_reference` | Optional scalar; required for finite-reference correction of supplied LD |
| `ld_origin` | Optional scalar description, such as `sample_reference` or `generating_truth` |

All arrays must follow the same SNP order as the LD matrix rows and columns or genotype columns. Metadata do not prove alignment: preserve ordering when exporting references. For genotype input, missing dosages are NaN. Infinities and dosages outside [0,2] are errors. Missing dosages are imputed to the observed mean before correlations are computed. All-missing and monomorphic SNPs are removed and counted. If an AF array is also supplied, a discrepancy exceeding 0.02 from genotype-derived A1 frequency produces an error.

Alleles with AF greater than 0.5 are flipped, their LD rows and columns are multiplied by -1, and their labels and AF are updated. At AF exactly 0.5, lexical allele order breaks the tie consistently. Palindromic SNPs remain in the reference scores because source variants contribute LD regardless of whether they can safely be harmonized as GWAS tags; the summary-statistic reader applies its own tag harmonization policy.

Supplied LD must be finite, symmetric, have a unit diagonal and entries in [-1,1], and be positive semidefinite unless that check is explicitly disabled. Multiple variants or blocks with the same SNP ID are rejected. An all-missing or fully filtered block produces an error requiring correction of the input.

## Finite-reference correction and its limits

With `--finite-reference`, the second-order score uses the conventional correction

```math
\widetilde{r^2}=r^2-\frac{1-r^2}{n_\mathrm{ref}-2}
```

for off-diagonal entries, retaining the self-term as exactly one. Negative corrected values are retained. This is an approximation for sample correlations, conditional on its sampling assumptions; it is not an exact correction for arbitrary reference sampling, relatedness, dosage imputation, or variant-dependent missingness. The complete reference sample count is used after mean imputation, rather than a different pairwise sample count for every LD pair.

**L1, L3, and L4 are not corrected for finite-reference bias.** Metadata state this explicitly, and downstream standard errors condition on the supplied reference. Increasing reference size or externally validating higher-order LD scores remains necessary for strong population-level claims. Genotype-derived LD is labeled `sample_genotypes`; supplied matrices are labeled `supplied_ld` unless their origin is supplied. Do not enable finite-reference correction for an LD matrix used as the exact generating covariance in a simulation. `ld_origin="generating_truth"` or `"population"` explicitly rejects that combination.

The output is a validated software product, not a prevalidated universal ancestry reference panel. Reference ancestry, genome build, allele coding, SNP coverage, block definitions, and sampling uncertainty must match the intended study.

## Python API

```python
from gcmsem.scores import build_scores, save_scores, load_scores

scores = build_scores(
    blocks,                         # iterable; one block held at a time
    maf_edges=(0, 0.05, 0.1, 0.5),
    tag_snps=hm3_ids,
    finite_reference=False,          # generating truth in this example
)
save_scores(scores, "reference.scores.npz")
restored = load_scores("reference.scores.npz")
```

Each block dictionary has the NPZ fields above and an `id` string. A scalar `block` can be used instead of `id`. The direct API accepts SciPy sparse LD matrices; these are converted to dense within a block. A `log` argument can receive a logger or a callable accepting a message.

## Output and audit information

The NPZ contains tag metadata (`snp`, `a1`, `a2`, `af`, `chrom`, `bp`, `block`), four tag × stratum arrays (`L1`–`L4`), source `counts`, and a JSON `metadata` scalar. `load_scores()` returns metadata as a dictionary. Archives use schema version `1.0`, no pickles, and atomic writes. The CLI refuses to replace an existing score archive unless `--overwrite` is supplied.

Metadata record source and tag counts, MAF definitions, input block sizes, per-block QC counts, reference sample size and LD origin, correction settings, warnings, software versions, and SHA-256 checksums of CLI input files. Requested tags absent after reference QC are counted and warned about. Blocks with no retained tags still contribute their source counts, under the shared within-stratum effect-distribution assumption. The command logs each block's progress and a final count summary. If a downstream analysis uses only a subset of tag rows, it must preserve the archive's full source counts.
