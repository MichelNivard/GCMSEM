"""Effective polygenicity from univariate, annotation-specific cumulants.

No comparator software or tutorial module is required. The LD average is a
tag-based proxy for source LD; uncertainty is conditional on the reference.
"""
import logging
import numpy as np
from .moments import jackknife_cov

METRICS = ('genetic_variance', 'fourth_cumulant_sum', 'ld_term', 'numerator',
           'denominator', 'effective_count', 'no_ld_denominator',
           'no_ld_count', 'excess_kurtosis', 'inverse_count')


def _metrics(theta, counts, mean_ld):
    v, k = theta[:, 1], theta[:, 3]
    h = counts @ v
    fourth = counts @ k
    ld = 3 * (counts * v) @ (mean_ld @ v)
    numerator = 3 * h**2
    denominator = fourth + ld
    no_ld = counts @ (k + 3 * v**2)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.array([h, fourth, ld, numerator, denominator,
            numerator / denominator, no_ld, numerator / no_ld,
            counts.sum() * no_ld / h**2 - 3, denominator / numerator])


def _covariance(values, point, fractions):
    # Singular ratios remain visible, with unavailable uncertainty, rather
    # than silently deleting or clipping individual jackknife estimates.
    valid = np.isfinite(point) & np.isfinite(values).all(axis=0)
    out = np.full((len(point), len(point)), np.nan)
    if valid.any():
        out[np.ix_(valid, valid)] = jackknife_cov(values[:, valid], point[valid], fractions)
    return out


def summarize_polygenicity(result, scores):
    """Return audit-ready summaries using every original genomic deletion."""
    theta = np.asarray(result['conditional'], float)
    reps = np.asarray(result['conditional_jackknife'], float)
    counts = np.asarray(scores['counts'], float)
    p = len(counts)
    if len(result['metadata']['traits']) != 1 or theta.shape != (p, 4):
        raise ValueError('Polygenicity requires all four univariate conditional cumulants')
    edges = np.asarray(scores['metadata']['maf_edges'], float)
    bins = np.minimum(np.searchsorted(edges, scores['af'], side='right') - 1, p - 1)
    L2 = np.asarray(scores['L2'], float)
    # Cache totals once; each deletion only reads its own rows.
    nb = np.bincount(bins, minlength=p)
    sums = np.array([np.bincount(bins, weights=L2[:, i], minlength=p) for i in range(p)]).T
    if np.any(nb == 0):
        raise ValueError('Every source MAF stratum needs representative tags')
    C = sums / nb[:, None]
    blocks = result['metadata']['jackknife_blocks']
    f = np.asarray(result['metadata']['jackknife_fractions'])
    Cj = []
    for block in blocks:
        start, stop = block['start_row'], block['stop_row']
        b = bins[start:stop]
        n = nb - np.bincount(b, minlength=p)
        if np.any(n == 0):
            raise ValueError('A deletion empties a MAF stratum; merge sparse score strata')
        removed = np.array([np.bincount(b, weights=L2[start:stop, i], minlength=p) for i in range(p)]).T
        Cj.append((sums - removed) / n[:, None])
    point = _metrics(theta, counts, C)
    deletions = np.array([_metrics(t, counts, c) for t, c in zip(reps, Cj)])
    V = _covariance(deletions, point, f)
    se = np.sqrt(np.maximum(np.diag(V), 0))
    # Joint covariance of all variance/fourth-cumulant estimates, including
    # between-stratum covariance. R also retains all four source cumulants.
    keys = [f'stratum-{c+1}:{q}' for c in range(p) for q in ('variance', 'fourth_cumulant')]
    Q = jackknife_cov(reps[:, :, [1, 3]].reshape(len(f), -1), theta[:, [1, 3]].ravel(), f)
    h, _, _, a, d = point[:5]
    W = np.diag(counts) @ C
    ga = np.zeros(2*p); ga[::2] = 6*h*counts
    gd = np.empty(2*p); gd[::2] = 3*(W+W.T) @ theta[:, 1]; gd[1::2] = counts
    delta_se = np.nan
    if np.isfinite(d) and d != 0:
        gradient = (ga*d - a*gd) / d**2
        delta_se = float(np.sqrt(max(gradient @ Q @ gradient, 0)))
    influence = np.full(len(f), np.nan)
    if np.isfinite(deletions[:, 5]).all() and np.isfinite(point[5]):
        pseudo = point[5] + (1/f-1)*(point[5]-deletions[:, 5])
        contribution = (pseudo-f@pseudo)**2*f/(1-f)/len(f)
        if contribution.sum() > 0:
            influence = contribution/contribution.sum()
    flags = []
    if np.any(theta[:, 1] <= 0): flags.append('nonpositive_bin_variance')
    if np.any(theta[:, 3] < 0): flags.append('negative_fourth_cumulant_gaussian_mixture_incompatible')
    positive = theta[:, 1] > 0
    if np.any(theta[positive, 3] < theta[positive, 2]**2/theta[positive, 1] - 2*theta[positive, 1]**2):
        flags.append('source_moment_inequality_violated')
    if h <= 0: flags.append('nonpositive_total_genetic_variance')
    if d <= 0: flags.append('nonpositive_denominator')
    if d - 1.96*se[4] <= 0: flags.append('denominator_not_separated_from_zero')
    if np.any(deletions[:, 4] <= 0): flags.append('nonpositive_denominator_in_deletions')
    if not np.isfinite(point).all() or not np.isfinite(deletions).all() or not np.isfinite(se).all(): flags.append('undefined_ratio_or_uncertainty')
    if np.isfinite(point[5]) and (point[5] <= 0 or point[5] > counts.sum()): flags.append('count_outside_source_range')
    if np.isfinite(influence).any() and np.nanmax(influence) > .5: flags.append('single_block_dominates_jackknife_variance')
    for flag in flags: logging.getLogger('gcmsem').warning('POLYGENICITY %s: %s', result['metadata']['traits'][0], flag)
    strata = []
    for c in range(p):
        cov = jackknife_cov(reps[:, c, :], theta[c], f)
        row = dict(stratum=c+1, maf_min=float(edges[c]), maf_max=float(edges[c+1]), source_snps=int(counts[c]), tag_snps=int(nb[c]))
        for j, name in enumerate(('mean', 'variance', 'third_cumulant', 'fourth_cumulant')):
            row[name] = float(theta[c, j]); row[name+'_se'] = float(np.sqrt(max(cov[j, j], 0)))
        strata.append(row)
    rows = []
    for b, block in enumerate(blocks):
        row = dict(block=b+1, chrom=str(block['chrom']),
            start_bp=int(scores['bp'][block['start_row']]), end_bp=int(scores['bp'][block['stop_row']-1]),
            n_tags=int(block['stop_row']-block['start_row']), fraction=float(f[b]), variance_fraction=float(influence[b]))
        row.update(zip(METRICS, map(float, deletions[b])))
        rows.append(row)
    summary = dict(trait=result['metadata']['traits'][0], effective_count=float(point[5]), se=float(se[5]),
        delta_se=delta_se, inverse_count=float(point[9]), inverse_count_se=float(se[9]),
        no_ld_count=float(point[7]), no_ld_se=float(se[7]), n_tags=len(bins), n_sources=int(counts.sum()),
        n_blocks=len(f), flags='; '.join(flags) if flags else 'none')
    logging.getLogger('gcmsem').info('POLYGENICITY count=%.8g SE=%.8g denominator=%.8g denominator_SE=%.8g', point[5], se[5], d, se[4])
    return dict(summary=summary, stages=[dict(quantity=k, estimate=float(x), se=float(s)) for k,x,s in zip(METRICS,point,se)],
        strata=strata, covariance=V.tolist(), covariance_keys=list(METRICS),
        conditional_covariance=Q.tolist(), conditional_covariance_keys=keys,
        ld_mean=C.tolist(), jackknife=rows, warnings=flags,
        assumptions=['Zero-mean independent-source Gaussian variance mixture interpretation',
            'Mean L2 among tags in each MAF bin approximates the mean among all source variants',
            'Reference scores and source counts are treated as fixed; LD-reference uncertainty is not propagated',
            'Unequal-group block jackknife; nonlinear-ratio coverage remains uncalibrated',
            'Signed estimates retained without clipping or dropping deletions',
            'delta_se uses the joint conditional cumulant covariance with mean LD held fixed'])


def json_safe(value):
    if isinstance(value, dict): return {k: json_safe(v) for k,v in value.items()}
    if isinstance(value, list): return [json_safe(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value): return None
    return value
