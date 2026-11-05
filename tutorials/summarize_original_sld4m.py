#!/usr/bin/env python3
"""Compare native S-LD4M base output with the unchanged GCMSEM analysis."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from gcmsem.moments import jackknife_cov
from polygenicity_helpers import polygenicity_metrics, TagLDAverages
from prepare_original_sld4m import sha256


def table(path, rows):
    with path.open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def summarize(gcmsem, native, upstream, out):
    out.mkdir(parents=True,exist_ok=True)
    a = json.loads((gcmsem/'analysis.json').read_text())
    r = json.loads((native/'reference-provenance.json').read_text())
    d = np.load(gcmsem/'design.npz'); counts = d['counts']
    ld = TagLDAverages(d['L2'],d['bins'],len(counts))
    rows = []; annotations = []; inputs = []; hashes = {}
    for trait in a['traits']:
        g = np.load(gcmsem/(trait+'.estimates.npz'))
        meta = json.loads((gcmsem/(trait+'.metadata.json')).read_text())
        theta = g['conditional']; point = polygenicity_metrics(theta,counts,ld.mean)
        deletions = np.array([polygenicity_metrics(t,counts,ld.delete(b['start_row'],b['stop_row']))
            for t,b in zip(g['conditional_jackknife'],meta['jackknife_blocks'])])
        se = np.sqrt(np.diag(jackknife_cov(deletions,point,meta['jackknife_fractions'])))
        s = loadmat(native/(trait+'.result.mat'),squeeze_me=True)
        if s['ma_jk'].dtype.kind != 'f' or s['ma_jk'].shape != (75,100):
            raise ValueError(f'{trait}: expected double-precision native deletion estimates')
        inputs.append(json.loads((native/(trait+'.input.json')).read_text()))
        flags = []
        if np.any(theta[:,1] <= 0): flags.append('nonpositive_GCMSEM_bin_variance')
        if np.any(theta[:,3] < 0): flags.append('negative_GCMSEM_bin_fourth_cumulant')
        if point[4]-1.96*se[4] <= 0: flags.append('GCMSEM_denominator_not_separated_from_zero')
        if np.any(deletions[:,4] <= 0): flags.append('GCMSEM_denominator_nonpositive_in_deletions')
        native_ma = float(s['ma'][0]); native_se = float(s['ma_se'][0])
        if native_ma <= 0: flags.append('nonpositive_SLD4M_base_count')
        if np.any(s['ma_jk'][0] <= 0): flags.append('SLD4M_base_count_nonpositive_in_deletions')
        if native_ma-1.96*native_se <= 0: flags.append('SLD4M_base_count_wald_interval_includes_zero')
        rows.append(dict(trait=trait,gcmsem=point[3],gcmsem_se=se[3],
            sld4m=native_ma,sld4m_se=native_se,
            gcmsem_no_ld=point[2],gcmsem_no_ld_se=se[2],
            excess_kurtosis=point[1],excess_kurtosis_se=se[1],
            gcmsem_tags=a['n_tags'],sld4m_matched_snps=inputs[-1]['matched_snps'],
            sld4m_native_rows_used=inputs[-1]['native_rows_used'],
            flags=';'.join(flags) or 'none'))
        for c,name in enumerate(r['annotations']):
            annotations.append(dict(trait=trait,annotation=name,index=c+1,
                ma=float(s['ma'][c]),ma_se=float(s['ma_se'][c])))
        hashes[trait] = dict(input_sha256=sha256(native/(trait+'.input.mat')),
                            output_sha256=sha256(native/(trait+'.result.mat')))
    x = np.array([row['sld4m'] for row in rows]); y = np.array([row['gcmsem'] for row in rows])
    stats = dict(n_traits=len(rows),pearson=float(np.corrcoef(x,y)[0,1]),
        rmse=float(np.sqrt(np.mean((x-y)**2))),mean_difference=float(np.mean(y-x)),
        sld4m_input='Raw GWAS chi-square',sld4m_reporting='Native ma(1), ma_se(1); no custom aggregation',
        reference='Author-released 1000 Genomes baselineLD.1kg.l2l4.mat',
        reference_rows=r['reference_rows'],n_annotations=len(r['annotations']),n_blocks=100,
        source_modified=False,sld4m_sha256=sha256(upstream/'SLD4M.m'),
        prior_matched_input_pearson=0.9651203420917207,
        interpretation='Independent pipeline comparison with different reference SNP universes; descriptive agreement, not a truth validation')
    table(out/'comparison.tsv',rows); table(out/'annotations.tsv',annotations)
    (out/'comparison.json').write_text(json.dumps(stats,indent=2)+'\n')
    provenance = dict(reference=r,traits=inputs,hashes=hashes,gcmsem_analysis=a,
        sld4m_commit='dd947104cb3cb86bb8a651bf3b5abf5e5cbf46f3',
        sld4m_sha256=stats['sld4m_sha256'])
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(stats,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['gcmsem','native','upstream','out']: p.add_argument('--'+name,type=Path,required=True)
    a = p.parse_args(); summarize(a.gcmsem,a.native,a.upstream,a.out)
