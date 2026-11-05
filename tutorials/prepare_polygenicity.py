#!/usr/bin/env python3
"""Estimate univariate GCMSEM moments and prepare matched original S-LD4M inputs.

Large intermediate files stay in --out. No upstream implementation is copied
or modified. Run from a checkout with the GCMSEM Python package installed.
"""
import argparse
import json
import logging
from pathlib import Path
import platform
import time

import numpy as np
from scipy.io import savemat
from scipy.sparse import csc_matrix
from gcmsem.scores import load_scores
from gcmsem.io import harmonize_sumstats, read_manifest, sha256_file, write_results
from gcmsem.moments import estimate_moments


def prepare(data_dir, out, n_blocks=200):
    data_dir, out = Path(data_dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'analysis.json').exists():
        raise FileExistsError('Choose a new output directory; analysis.json already exists')
    log = logging.getLogger('gcmsem')
    log.setLevel(logging.INFO)
    handler = logging.FileHandler(out / 'preparation.log')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    log.addHandler(handler)
    start = time.monotonic()
    score_path = data_dir / 'scores/eur_hm3.scores.npz'
    scores = load_scores(score_path)
    entries = read_manifest(data_dir / 'manifests/panel.tsv')
    for entry in entries:
        if entry['trait'] == 'Alcohol':
            entry['trait'] = 'DrinksPerWeek'
    beta, se, sc, traits, qc = harmonize_sumstats(entries, scores, logger=log)
    # Upstream uses floor(n/B) and silently ignores a final partial block.
    # Remove that tail explicitly from BOTH methods and record the count.
    original_n = len(beta)
    n = original_n - original_n % n_blocks
    if n < n_blocks * 2:
        raise ValueError('Too few tags for the requested jackknife blocks')
    beta, se = beta[:n], se[:n]
    for key in ('snp', 'a1', 'a2', 'af', 'chrom', 'bp', 'block', 'L1', 'L2', 'L3', 'L4'):
        sc[key] = np.asarray(sc[key])[:n]
    counts = np.asarray(sc['counts'], dtype=int)
    # The released score archive has these five disjoint MAF strata.
    edges = np.array([0, .05, .1, .2, .3, .5])
    if len(counts) != 5:
        raise ValueError('This worked example requires the released five-stratum scores')
    bins = np.minimum(np.searchsorted(edges, sc['af'], side='right') - 1, 4)
    nt = np.bincount(bins, minlength=5)
    remainder = counts - nt
    if np.any(remainder < 0) or np.any(nt == 0):
        raise ValueError('Invalid tag/source stratum counts')
    # SLD4M needs annotations for tags and the SOURCE counts. Put the aligned
    # tag annotation rows first and append count-only rows for unobserved sources.
    # idxldsc selects only the real tag rows; no dummy GWAS observations are added.
    source_bins = np.r_[bins, np.repeat(np.arange(5), remainder)]
    annot = csc_matrix((np.ones(len(source_bins)), (np.arange(len(source_bins)), source_bins)),
                       shape=(int(counts.sum()), 5))
    weights = 1 / np.maximum(sc['L2'].sum(axis=1), 1)
    savemat(out / 'sld4m_input.mat', dict(beta=beta, se=se, ell2=sc['L2'], ell4=sc['L4'],
        annot=annot, weights=weights[:, None], idxldsc=np.arange(1, n + 1, dtype=float)[:, None],
        no_blocks=n_blocks, counts=counts[:, None], tagbins=(bins + 1)[:, None]))
    np.savez(out / 'design.npz', L2=sc['L2'], bins=bins, counts=counts,
             chrom=sc['chrom'], block=sc['block'])
    log.info('MATCHED_INPUT tags=%d removed_tail=%d source_variants=%d traits=%d',
             n, original_n - n, counts.sum(), len(traits))
    results = []
    for i, trait in enumerate(traits):
        log.info('TRAIT %s (%d/%d)', trait, i + 1, len(traits))
        result = estimate_moments(beta[:, i:i+1], se[:, i:i+1], sc, traits=[trait],
                                  max_order=4, n_blocks=n_blocks, logger=log)
        prefix = out / trait
        write_results(result, prefix, overwrite=False)
        np.savez(out / (trait + '.estimates.npz'), conditional=result['conditional'],
                 conditional_jackknife=result['conditional_jackknife'])
        results.append(dict(trait=trait, noise_intercept=float(result['noise_intercept'][0, 0])))
    record = dict(traits=traits, n_tags_before_tail_removal=original_n, n_tags=n,
        removed_tail=original_n - n, source_counts=counts.tolist(), tag_counts=nt.tolist(),
        maf_edges=edges.tolist(), n_blocks=n_blocks, qc=qc, score_sha256=sha256_file(score_path),
        score_metadata=sc['metadata'], trait_results=results, python=platform.python_version(),
        elapsed_seconds=time.monotonic() - start,
        comparison='same tags, source counts, L2/L4 and fixed inverse-L2 weights',
        source_ld_average='HM3-tag mean L2 within each MAF bin; representative-tag approximation')
    (out / 'analysis.json').write_text(json.dumps(record, indent=2) + '\n')
    log.info('PREPARATION_DONE elapsed_seconds=%.1f', record['elapsed_seconds'])
    handler.close(); log.removeHandler(handler)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', default='tutorial-data')
    parser.add_argument('--out', default='tutorial-results/polygenicity')
    parser.add_argument('--n-blocks', type=int, default=200)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    prepare(args.data_dir, args.out, args.n_blocks)
