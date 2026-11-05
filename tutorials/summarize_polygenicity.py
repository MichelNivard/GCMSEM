#!/usr/bin/env python3
"""Read freshly executed GCMSEM/S-LD4M results and write small comparison tables."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from gcmsem.moments import jackknife_cov
from polygenicity_helpers import polygenicity_metrics, combine_sld4m_bins, TagLDAverages


def write_table(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def summarize(directory, upstream):
    out, upstream = Path(directory), Path(upstream)
    a = json.loads((out / 'analysis.json').read_text())
    d = np.load(out / 'design.npz'); counts = d['counts']; M = counts.sum()
    ld = TagLDAverages(d['L2'], d['bins'], len(counts))
    beta_fit = loadmat(out / 'sld4m.mat')
    z_path = out / 'sld4m_zscore.mat'
    if not z_path.exists() and 'ld_reference' not in a:
        raise FileNotFoundError('Run the primary Z-squared analysis: sld4m_zscore.mat is missing')
    s = loadmat(z_path) if z_path.exists() else beta_fit
    B = a['n_blocks']; n = len(d['bins'])
    # Upstream averages each equal-row block's within-bin means equally;
    # that differs from pooling rows when bin frequencies vary by block.
    group_C = np.array([TagLDAverages(d['L2'][b*(n//B):(b+1)*(n//B)],
                        d['bins'][b*(n//B):(b+1)*(n//B)], len(counts)).mean for b in range(B)])
    native_C = (group_C.sum(axis=0) - group_C) / (B - 1)
    summary, strata = [], []
    for i, trait in enumerate(a['traits']):
        g = np.load(out / (trait + '.estimates.npz'))
        meta = json.loads((out / (trait + '.metadata.json')).read_text())
        theta = g['conditional']; point = polygenicity_metrics(theta, counts, ld.mean)
        j = np.array([polygenicity_metrics(t, counts, ld.delete(b['start_row'], b['stop_row']))
                      for t, b in zip(g['conditional_jackknife'], meta['jackknife_blocks'])])
        se = np.sqrt(np.diag(jackknife_cov(j, point, meta['jackknife_fractions'])))
        hj, mj = s['variance_jk'][:, :, i], s['ma_jk'][:, :, i]
        sj = combine_sld4m_bins(hj, mj)
        # Preserve the original code's convention: std(..., ddof=1)*sqrt(B+1).
        native_se = sj.std(ddof=1) * np.sqrt(B + 1)
        beta_j = combine_sld4m_bins(beta_fit['variance_jk'][:, :, i], beta_fit['ma_jk'][:, :, i])
        sd = 3 * np.sum(hj**2 / mj, axis=0)
        sd_se = sd.std(ddof=1) * np.sqrt(B + 1)
        def diagnostics(native_counts, denominator):
            flags = []
            if np.any(theta[:, 1] <= 0): flags.append('nonpositive_GCMSEM_bin_variance')
            if np.any(theta[:, 3] < 0): flags.append('negative_GCMSEM_bin_fourth_cumulant')
            if np.any(native_counts <= 0): flags.append('nonpositive_SLD4M_bin_count')
            if point[4] - 1.96 * se[4] <= 0: flags.append('GCMSEM_denominator_not_separated_from_zero')
            denominator_se = denominator.std(ddof=1) * np.sqrt(B + 1)
            if denominator.mean() - 1.96 * denominator_se <= 0: flags.append('SLD4M_denominator_not_separated_from_zero')
            if np.any(j[:, 4] <= 0) or np.any(denominator <= 0): flags.append('denominator_sign_change_in_deletions')
            return ';'.join(flags) or 'none'
        flags = diagnostics(s['ma'][:, i], sd)
        beta_sd = 3 * np.sum(beta_fit['variance_jk'][:, :, i]**2 / beta_fit['ma_jk'][:, :, i], axis=0)
        summary.append(dict(trait=trait, gcmsem=point[3], gcmsem_se=se[3],
            sld4m=sj.mean(), sld4m_se=native_se,
            sld4m_beta=beta_j.mean(), sld4m_beta_se=beta_j.std(ddof=1)*np.sqrt(B+1),
            gcmsem_no_ld=point[2], no_ld_se=se[2],
            excess_kurtosis=point[1], excess_kurtosis_se=se[1],
            genetic_variance=point[0], genetic_variance_se=se[0],
            gcmsem_ld_denominator=point[4], gcmsem_ld_denominator_se=se[4],
            sld4m_ld_denominator=sd.mean(), sld4m_ld_denominator_se=sd_se,
            n_tags=n, n_sources=int(M), flags=flags,
            flags_beta=diagnostics(beta_fit['ma'][:, i], beta_sd)))
        # Reconstruct upstream K_c from native per-bin Ma/H outputs. This is
        # its fourth-regression coefficient, which need not equal a causal
        # fourth cumulant when the two methods' assumptions fail.
        beta_hj = beta_fit['variance_jk'][:, :, i]; beta_mj = beta_fit['ma_jk'][:, :, i]
        vj = beta_hj / counts[:, None]
        native_k = np.column_stack([3 * counts * vj[:, b]**2 / beta_mj[:, b]
                                   - 3 * vj[:, b] * (native_C[b] @ vj[:, b]) for b in range(B)])
        for c in range(len(counts)):
            v, k = theta[c, 1], theta[c, 3]
            gp = np.array([v, k])
            gj = g['conditional_jackknife'][:, c, :][:, [1, 3]]
            gse = np.sqrt(np.diag(jackknife_cov(gj, gp, meta['jackknife_fractions'])))
            strata.append(dict(trait=trait, stratum=c + 1, source_count=int(counts[c]),
                gcmsem_variance=v, gcmsem_variance_se=gse[0], gcmsem_fourth_cumulant=k,
                gcmsem_fourth_cumulant_se=gse[1], sld4m_beta_fourth_coefficient=native_k[c].mean(),
                sld4m_beta_fourth_coefficient_se=native_k[c].std(ddof=1) * np.sqrt(B + 1),
                sld4m_effective_count=s['ma'][c, i], sld4m_effective_count_se=s['ma_se'][c, i],
                sld4m_beta_effective_count=beta_fit['ma'][c, i],
                sld4m_beta_effective_count_se=beta_fit['ma_se'][c, i]))
    write_table(out / 'comparison.tsv', summary); write_table(out / 'strata.tsv', strata)
    x = np.array([r['sld4m'] for r in summary]); y = np.array([r['gcmsem'] for r in summary])
    stats = dict(n_traits=len(x), pearson=float(np.corrcoef(x, y)[0, 1]),
                 rmse=float(np.sqrt(np.mean((x-y)**2))),
                 mean_difference=float(np.mean(y-x)),
                 sld4m_input='Z squared' if z_path.exists() else 'BETA squared (constant simulation SE)',
                 beta_input_pearson=float(np.corrcoef(y, [r['sld4m_beta'] for r in summary])[0,1]),
                 beta_input_rmse=float(np.sqrt(np.mean((y-np.array([r['sld4m_beta'] for r in summary]))**2))),
                 sld4m_sha256=hashlib.sha256((upstream / 'SLD4M.m').read_bytes()).hexdigest(),
                 upstream='https://github.com/lukejoconnor/SLD4M',
                 estimand='LD-adjusted Gaussian-mixture effective count; source-weighted bin aggregation; ' +
                          ('exact source LD averages' if 'ld_reference' in a else 'tag LD proxy'),
                 uncertainty='method-specific block jackknife; conditional on scores and source counts')
    (out / 'comparison.json').write_text(json.dumps(stats, indent=2) + '\n')
    print(json.dumps(stats, indent=2))
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', default='tutorial-results/polygenicity')
    p.add_argument('--sld4m', default='external/SLD4M')
    a = p.parse_args(); summarize(a.out, a.sld4m)
