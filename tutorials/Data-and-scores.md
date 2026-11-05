# Data, score definitions and reproducibility

The release separates GCMSEM's signed LD-power scores from the **published standard LDSC reference**. GCMSEM requires signed first, third and fourth powers that conventional LDSC files do not contain. The primary LDSC comparison uses its own stock software, standard published scores/weights, and ordinary analysis procedure.

## GCMSEM reference

| Property | Released value |
|---|---|
| Population / build | 503 European individuals; GRCh37 |
| Reference source | Alan Aw's [1000 Genomes reference, Zenodo 6614170](https://zenodo.org/records/6614170) |
| Regression tags | 1,074,062 harmonized HapMap3 SNPs |
| Source universe | 1,719,440 QC-passed variants |
| MAF boundaries | 0, 0.05, 0.10, 0.20, 0.30, 0.50 |
| Source counts by bin | 464,696; 205,404; 314,059; 261,563; 473,718 |
| LD neighborhood | Same chromosome, absolute distance at most 1,000,000 bp |
| Orientation | Minor-allele-coded, genotype-standardized signed LD |
| Genomic deletions | 200 chromosome-respecting groups, with unequal-size adjustment |
| Finite-reference correction | L2 corrected once; L1, L3, L4 uncorrected |

For a tag and source MAF stratum, each Lq score is the sum of signed LD raised to power q over that stratum's source variants. Thus L3 preserves negative contributions. A conventional r²-only score cannot be used to reconstruct the odd-order scores.

The scores were originally computed with a rolling window and imported losslessly into the package's validated archive schema. All **21,481,240** L1–L4 entries were then independently recomputed from checksum-verified reference files. The rebuild passed the recorded absolute-plus-relative tolerance of `3e-5`; small float32 linear algebra differences mean this is numerical reproduction, not a claim of byte-identical recomputation. The release preserves the original numerical array bytes so previous analyses remain reproducible.

The score's `block` field records deletion groups. It does **not** claim that LD vanishes between groups. This reference differs from `gcmsem score`'s explicitly block-truncated construction. Metadata and the tutorial logs preserve that difference.

## Rebuild without redistributing reference genotypes

Ordinary tutorial users download the derived scores and do not need raw reference genotypes. A maintainer can independently rebuild them:

```bash
OPENBLAS_NUM_THREADS=1 python tutorials/rebuild_empirical_scores.py \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --frequency-exclusions tutorial-data/provenance/frequency_outliers.tsv \
  --reference-dir reference-work --download \
  --out rebuilt/eur_hm3.scores.npz
```

This command downloads the original reference from its source into `reference-work/`, checks the published file MD5 values and the reconstructed source-index SHA-256, and writes derived scores. Raw BED/BIM/FAM files stay local and are excluded from the release builder. No participant genotypes or copied HM3 reference archive are distributed by GCMSEM.

The earlier matched-reference sensitivity has its own HM3-only weighting scores, which can also be rebuilt:

```bash
OPENBLAS_NUM_THREADS=1 python tutorials/build_ldsc_weights.py \
  --scores tutorial-data/scores/eur_hm3.scores.npz \
  --reference-dir reference-work --download \
  --out rebuilt/HM3_weights.npz
```

**Those recomputed weights are only for the explicitly labeled matched-reference sensitivity. They are not used in the primary standard LDSC comparison.**

## Release contents and integrity

| Asset ID | Contents |
|---|---|
| `scores` | GCMSEM signed scores, secondary matched-reference weights, metadata and score-rebuild audit |
| `causal` | Four SNP-level summary-statistic subsets and height/weight, BMI/SBP manifests |
| `panel` | Fifteen-trait SNP-level subsets, prealigned GCMSEM common panel, trait manifest and provenance |
| `results` | Executed moment/SEM bundles, correlation matrices, stock LDSC logs, diagnostics and primary comparison outputs |

The actual archive list and exact sizes are in [data_manifest.json](data_manifest.json). Archives contain their own per-file inventories. `download_data.py` verifies archive size and SHA-256 before extraction, rejects unsafe paths, and refuses to overwrite modified local files. A DOI identifies a scholarly release; the checksum identifies its exact bytes.

The standard LDSC inputs and unchanged published reference files are documented by the standard-comparison workflow. The original Broad URLs may move; any archived mirror must retain the standard numerical files and exact checksum provenance. Our custom signed reference is never substituted when a standard LDSC download fails.

See [DATA_RIGHTS.md](DATA_RIGHTS.md) for original-source attribution, mixed source rights, and the public GSCAN release's exclusion of 23andMe. The scripts and documentation are in GitHub; all numerical archives are public on [Zenodo](https://zenodo.org/records/23977299), DOI `10.5281/zenodo.23977299`. Follow the [download instructions](README.md#data-release-status) to obtain the checksum-pinned scores and GWAS files.
