# Tutorial data: rights and provenance

**GCMSEM empirical data release 1.0.0:** [DOI 10.5281/zenodo.23977299](https://zenodo.org/records/23977299).
The numerical archives are publicly available on Zenodo, separately from GitHub.
The source DOIs below identify the upstream datasets.

The tutorial archive contains derived LD scores, harmonized public GWAS summary
statistics, analysis results and provenance. It does **not** contain the HM3/1000
Genomes reference genotypes, BED/BIM/FAM files, participant-level phenotypes or
identifiers. SNP identifiers, alleles and reference frequencies accompanying a
derived score are variant metadata. The scores describe LD moments; they are
not polygenic-score weights.

## Per-file rights

| Archive material | Rights and attribution |
|---|---|
| Derived GCMSEM score tables and newly authored result tables | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attribute Michel G Nivard/GCMSEM, the released dataset version and the upstream LD reference below. This licence covers the newly generated material, not a relicensing of third-party GWAS data. |
| `sumstats/GSCAN_DrinksPerWeek.tsv.gz` | GSCAN 2019 public release, [dataset DOI 10.13020/3b1n-ff32](https://doi.org/10.13020/3b1n-ff32), [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). Cite Liu et al. (2019), [10.1038/s41588-018-0307-5](https://doi.org/10.1038/s41588-018-0307-5). The public release excludes 23andMe. |
| Other `sumstats/*.tsv.gz` files | Neale Lab UK Biobank round 2. Its [official FAQ](https://www.nealelab.is/uk-biobank/faq) permits unrestricted use, including commercial use. Preserve that upstream policy and cite the [round-2 release](https://www.nealelab.is/uk-biobank/). These source data are not labelled CC BY by GCMSEM. |
| Repository code | The repository's [MIT licence](../LICENSE). |

The public GSCAN release used here is the corrected **23 September 2019** file.
It includes UK Biobank, so it overlaps the Neale traits. Its released effects
are reconstructed standardized-score effects, not coefficients in raw
drinks/week. Use its per-SNP sample size rather than the paper's headline
full-consortium sample size. The older Broad mirror is not interchangeable.

The LD reference is Alan Aw (2022), *1000 Genomes (Phase 3, GRCh37 — split into
five superpopulations)*, [10.5281/zenodo.6614170](https://doi.org/10.5281/zenodo.6614170),
under CC BY 4.0. The analysis uses its European panel (503 individuals), with
HM3 regression SNPs and the documented source-variant universe. Credit the
[1000 Genomes Project](https://www.internationalgenome.org/) as well as this
reference distribution. The numerical reference files are fetched upstream
only when rebuilding scores; they are excluded from our release.

## Exact phenotype inventory

The 15-trait comparison uses these codes in this order:

| Code | Trait | Source phenotype scale |
|---|---|---|
| `GSCAN_DrinksPerWeek` | Alcohol drinks/week | GSCAN transformed, reconstructed standardized score |
| `30890_raw` | Vitamin D | Raw nmol/L |
| `4080_irnt` | Systolic blood pressure | Inverse-normal transformed |
| `1200` | Insomnia | Ordinal |
| `20544_11` | Professionally diagnosed depression | Binary |
| `30780_raw` | Direct LDL cholesterol | Raw mmol/L |
| `3063_irnt` | FEV1 | Inverse-normal transformed |
| `6150_1` | Heart attack | Binary |
| `6150_2` | Angina | Binary |
| `50_irnt` | Height | Inverse-normal transformed |
| `21002_irnt` | Weight | Inverse-normal transformed |
| `21001_irnt` | BMI | Inverse-normal transformed |
| `4079_irnt` | Diastolic blood pressure | Inverse-normal transformed |
| `3062_irnt` | FVC | Inverse-normal transformed |
| `48_irnt` | Waist circumference | Inverse-normal transformed |

The causal tutorial additionally includes `21001_raw` (BMI in kg/m²) and
`4080_raw` (SBP in mmHg). Height–weight uses `50_irnt` and `21002_irnt`.
Do not interchange raw and inverse-normal transformed files. A causal model
coefficient depends on this phenotype scale.

Every exported file is a filtered, harmonized derivative, not the original
full GWAS release: effects and SEs are on a reference genotype-SD scale and
effects point to the reference minor allele. `REF_AF` is a reference-panel
frequency, not the GWAS study's allele frequency. Consult the tutorial's
manifest for columns and the genotype scale before using these files elsewhere.

## Provenance and integrity

The external archive's `provenance/traits/CODE.json` records each original
source URL, source release, phenotype transformation and original-file
checksums. `provenance/panel_manifest.json`, `height_weight_sources.json` and
`reference_manifest.json` retain the panel selection and reference provenance.
Original source-file checksums identify the upstream inputs; they do **not**
identify the filtered `.tsv.gz` outputs. The release download manifest supplies
separate SHA-256 hashes for the actual downloadable files.

Rights metadata were verified on **5 November 2026** through the following
primary APIs. These hashes identify the saved JSON evidence snapshots, not the
GWAS or genotype data:

| Evidence snapshot | Source | SHA-256 |
|---|---|---|
| `gscan-2019-datacite.json` | [DataCite registered dataset metadata](https://api.datacite.org/dois/10.13020/3b1n-ff32), `rightsList`: CC0 1.0 Universal | `cf638f0ef26e015e420fc814d5616fd3b10bb0166c7696a481289b5f3eb33bca` |
| `reference-6614170-zenodo.json` | [Zenodo record metadata](https://zenodo.org/api/records/6614170), `metadata.license.id`: `cc-by-4.0` | `f67eadef0311e837019a91fbb900d11d3f8e3b3916cbc666966e6e4776fcff55` |

Archive publication must retain these distinctions. A single default licence
selection on a hosting service must not override the per-file rights above.
