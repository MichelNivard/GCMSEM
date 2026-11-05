#!/usr/bin/env python3
"""Convert the archived real-GWAS inputs to portable tutorial inputs.

This is a provenance-preserving migration of existing rolling-window scores,
not a claim that block-truncated gcmsem score regenerated those quantities.
No participant genotypes, raw HM3 list, or credentials are copied to output.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from gcmsem.scores import save_scores, load_scores

PANEL_CODES = ['GSCAN_DrinksPerWeek', '30890_raw', '4080_irnt', '1200',
               '20544_11', '30780_raw', '3063_irnt', '6150_1', '6150_2',
               '50_irnt', '21002_irnt', '21001_irnt', '4079_irnt',
               '3062_irnt', '48_irnt']
PANEL_NAMES = ['Alcohol', 'VitaminD', 'SBP', 'Insomnia', 'Depression', 'LDL',
               'FEV1', 'HeartAttack', 'Angina', 'Height', 'Weight', 'BMI',
               'DBP', 'FVC', 'Waist']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def manifest(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['trait', 'path', 'effect_scale'])
        for name, code in rows:
            w.writerow([name, '../sumstats/' + code + '.tsv.gz', 'genotype_standardized'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, required=True,
                    help='LD_moments_15_trait_cache.npz')
    ap.add_argument('--pilot', type=Path, required=True,
                    help='Extracted height_weight_run directory')
    ap.add_argument('--walkthrough', type=Path, required=True,
                    help='Extracted walkthrough directory')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--skip-tsv', action='store_true', help='Only produce score and aligned NPZ inputs')
    ap.add_argument('--overwrite', action='store_true', help='Replace previously generated tutorial inputs')
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if not args.overwrite and any((out / p).exists() for p in ['scores/eur_hm3.scores.npz', 'panel.npz']):
        raise SystemExit('Tutorial inputs already exist; choose a fresh output directory or pass --overwrite.')
    score_path = args.pilot / 'artifacts/height_weight_LD_scores.npz'
    aligned_path = args.pilot / 'artifacts/height_weight_HM3_aligned.npz'
    original_meta = json.loads((args.pilot / 'ld_alltargets/ld_scores_metadata.json').read_text())
    deletion_meta = json.loads((args.pilot / 'results_alltargets/analysis_metadata.json').read_text())
    with np.load(args.cache, allow_pickle=False) as z:
        cache = {k: z[k] for k in z.files}
    with np.load(score_path, allow_pickle=False) as z:
        original_scores = {k: z[k] for k in z.files}
    with np.load(aligned_path, allow_pickle=False) as z:
        aligned = {k: z[k] for k in z.files}
    n = len(cache['reference_rsid'])
    for key in ['bed_index', 'chr', 'bp', 'rsid', 'a1', 'a2', 'minor_allele', 'minor_af_ref']:
        np.testing.assert_array_equal(cache['reference_' + key], aligned[key])
    np.testing.assert_array_equal(cache['reference_bed_index'], original_scores['bed_index'])
    codes = list(cache['trait_codes'].astype(str))
    for code, j in [('50_irnt', 0), ('21002_irnt', 1)]:
        np.testing.assert_array_equal(cache['beta_std'][:, codes.index(code)], aligned['beta_std'][:, j])
        np.testing.assert_array_equal(cache['se_std'][:, codes.index(code)], aligned['se_std'][:, j])
    snp, a1 = cache['reference_rsid'], cache['reference_minor_allele']
    a2 = np.where(a1 == cache['reference_a1'], cache['reference_a2'], cache['reference_a1'])
    blocks = np.full(n, '', dtype='<U20')
    deletion_groups = []
    for k, (chrom, start, end) in enumerate(deletion_meta['blocks']):
        assert start < end and np.all(cache['reference_chr'][start:end] == chrom)
        assert np.all(blocks[start:end] == '')
        label = f'chr{chrom}:delete{k + 1:03d}'
        blocks[start:end] = label
        deletion_groups.append({'id': label, 'chrom': str(chrom), 'start_row': start,
                                'end_row_exclusive': end, 'tag_variants': end - start})
    assert np.all(blocks != '') and len(deletion_groups) == 200
    source_hashes = {str(p.name): digest(p) for p in [args.cache, score_path, aligned_path]}
    score_meta = {
        'schema_version': '1.0', 'kind': 'signed_ld_power_scores',
        'effect_scale': 'genotype_standardized', 'allele_orientation': 'minor_allele',
        'source_variants': int(original_scores['target_bin_counts'].sum()),
        'tag_variants': n, 'maf_edges': [0, .05, .1, .2, .3, .5],
        'l2_finite_reference_correction': True, 'reference_sample_size': 503,
        'l1_l3_l4_finite_reference_correction': False,
        'genome_build': 'GRCh37', 'reference_ancestry': 'EUR',
        'ld_model': 'sliding_window', 'window_bp': 1000000,
        'cross_block_ld_assumed_zero': False,
        'block_semantics': 'Legacy 200 chromosome-respecting genomic deletion groups; not independent LD truncation blocks.',
        'migration': 'Exact archived L1,L2,L3,L4 arrays imported into GCMSEM schema; L2 correction is already applied once.',
        'source_hashes': source_hashes,
        'reference_source': 'https://zenodo.org/records/6614170',
        'reference_license': 'CC-BY-4.0',
        'original_score_metadata': original_meta,
        'blocks': deletion_groups,
        'warnings': [
            'Scores were computed earlier from signed 503-person EUR reference LD using a 1 Mb rolling window, not the block-truncated GCMSEM score command.',
            'Deletion boundaries can cut LD neighborhoods; jackknife uncertainty is approximate and conditions on the reference.',
            'Finite-reference correction applies to L2 only; higher-order reference uncertainty is not propagated.',
        ],
    }
    scores = dict(snp=snp, a1=a1, a2=a2, af=cache['reference_minor_af_ref'],
                  chrom=cache['reference_chr'].astype(str), bp=cache['reference_bp'], block=blocks,
                  counts=original_scores['target_bin_counts'], metadata=score_meta)
    scores.update({f'L{q}': original_scores[f'L{q}'] for q in range(1, 5)})
    save_scores(scores, out / 'scores/eur_hm3.scores.npz')
    check = load_scores(out / 'scores/eur_hm3.scores.npz')
    for q in range(1, 5):
        np.testing.assert_array_equal(check[f'L{q}'], original_scores[f'L{q}'])
    del check
    write_json(out / 'provenance/score_metadata.json', score_meta)

    panel_cols = np.array([codes.index(code) for code in PANEL_CODES])
    keep = cache['valid'][:, panel_cols].all(axis=1)
    row_index = np.flatnonzero(keep)
    assert len(row_index) == 1012060
    input_meta = {
        'kind': 'trusted_archived_harmonized_gwas', 'genome_build': 'GRCh37',
        'allele_orientation': 'reference_minor_allele', 'effect_scale': 'genotype_standardized',
        'source_hashes': source_hashes,
        'alignment': 'Exact original coordinate/allele-pair matches; no strand-complement guesses; source-reference AF discrepancy <=0.1 checked at acquisition.',
        'EAF': 'Original study effect-allele frequencies are not retained. Reference AF is not relabeled as study EAF.',
        'palindromic_policy': 'Direct-array reproduction preserves trusted archived strand orientation. Generic TSV CLI harmonization drops palindromic SNPs because study EAF is unavailable.',
        'panel_rows': int(keep.sum()), 'original_rows': n,
    }
    np.savez_compressed(out / 'panel.npz',
                        beta_std=cache['beta_std'][np.ix_(row_index, panel_cols)],
                        se_std=cache['se_std'][np.ix_(row_index, panel_cols)],
                        n=cache['n'][np.ix_(row_index, panel_cols)],
                        row_index=row_index, trait_codes=np.array(PANEL_CODES),
                        trait_names=np.array(PANEL_NAMES), chrom=cache['reference_chr'][keep],
                        bp=cache['reference_bp'][keep], snp=snp[keep], a1=a1[keep], a2=a2[keep],
                        af=cache['reference_minor_af_ref'][keep], block=blocks[keep],
                        metadata=np.array(json.dumps(input_meta, sort_keys=True)))
    write_json(out / 'provenance/harmonization.json', input_meta)
    for code in codes:
        p = args.walkthrough / 'expanded_analysis/source_research' / (code + '.json')
        if p.exists():
            write_json(out / 'provenance/traits' / (code + '.json'), json.loads(p.read_text()))
    for name, path in [
        ('panel_manifest.json', args.walkthrough / 'expanded_analysis/source_research/panel_manifest.json'),
        ('reference_manifest.json', args.pilot / 'reference/reference_manifest.json'),
        ('height_weight_sources.json', args.pilot / 'source_research/source_metadata.json'),
    ]:
        write_json(out / 'provenance' / name, json.loads(path.read_text()))

    manifests = out / 'manifests'
    manifest(manifests / 'height_weight.tsv', [('Height', '50_irnt'), ('Weight', '21002_irnt')])
    manifest(manifests / 'bmi_sbp.tsv', [('BMI', '21001_raw'), ('SBP', '4080_raw')])
    manifest(manifests / 'bmi_sbp_irnt.tsv', [('BMI', '21001_irnt'), ('SBP', '4080_irnt')])
    manifest(manifests / 'panel.tsv', list(zip(PANEL_NAMES, PANEL_CODES)))
    print(f'Prepared scores and panel ({len(row_index):,} rows); exporting summary tables.', flush=True)
    if not args.skip_tsv:
        (out / 'sumstats').mkdir(exist_ok=True)
        priority = ['50_irnt', '21002_irnt', '21001_raw', '4080_raw']
        export_codes = priority + [code for code in codes if code not in priority]
        for code in export_codes:
            j = codes.index(code)
            path = out / 'sumstats' / (code + '.tsv.gz')
            temporary = path.with_suffix(path.suffix + '.partial')
            with temporary.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding='utf-8', newline='') as text:
                    w = csv.writer(text, delimiter='\t', lineterminator='\n')
                    w.writerow(['SNP', 'A1', 'A2', 'BETA', 'SE', 'N', 'CHR', 'BP', 'REF_AF'])
                    for i in np.flatnonzero(cache['valid'][:, j]):
                        w.writerow([snp[i], a1[i], a2[i], format(cache['beta_std'][i, j], '.17g'),
                                    format(cache['se_std'][i, j], '.17g'), int(cache['n'][i, j]),
                                    int(cache['reference_chr'][i]), int(cache['reference_bp'][i]),
                                    format(cache['reference_minor_af_ref'][i], '.17g')])
            temporary.replace(path)
            print(f'{code}: {cache["valid"][:, j].sum():,} rows; {path.stat().st_size:,} bytes', flush=True)
    files = []
    for path in sorted(out.rglob('*')):
        if path.is_file() and path.name != 'checksums.json':
            files.append({'path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'sha256': digest(path)})
    write_json(out / 'checksums.json', {'schema_version': 1, 'files': files})
    print(f'Completed {len(files)} checksummed files at {out}', flush=True)


if __name__ == '__main__':
    main()
