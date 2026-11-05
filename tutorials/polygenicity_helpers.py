"""Small, independently testable transformations for the polygenicity tutorial."""
import numpy as np


def polygenicity_metrics(theta, counts, mean_ld):
    """Zero-mean Gaussian-variance-mixture interpretation of conditional cumulants.

    theta[:, 1] is variance and theta[:, 3] excess fourth cumulant. mean_ld[c,d]
    is the mean L2 to source bin d for variants in bin c. In empirical runs the
    latter is approximated by HM3-tag averages, not known for every source SNP.
    Signed algebraic estimates are retained. No clipping/feasibility projection.
    """
    theta, counts, mean_ld = map(np.asarray, (theta, counts, mean_ld))
    p = len(counts)
    if theta.shape != (p, 4) or mean_ld.shape != (p, p):
        raise ValueError('Expected four univariate cumulants per source stratum')
    if np.any(counts <= 0) or not all(np.isfinite(x).all() for x in (theta, counts, mean_ld)):
        raise ValueError('Counts must be positive and inputs finite')
    v, k = theta[:, 1], theta[:, 3]
    M = counts.sum(); H = counts @ v
    no_ld_denominator = counts @ (k + 3 * v**2)
    ld_denominator = counts @ k + 3 * (counts * v) @ (mean_ld @ v)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.array([H, M * no_ld_denominator / H**2 - 3,
                         3 * H**2 / no_ld_denominator, 3 * H**2 / ld_denominator,
                         ld_denominator])


def combine_sld4m_bins(genetic_variance, effective_count):
    """Source-weighted total, computed for EACH upstream deletion, then averaged.

    H_c = M_c v_c and D_c = 3 H_c^2 / M_e,c. Thus M_e,total =
    (sum_c H_c)^2 / sum_c(H_c^2/M_e,c). This is an external transformation
    of native per-bin outputs, not the upstream 'all SNPs' reporting shortcut.
    """
    h, m = np.asarray(genetic_variance), np.asarray(effective_count)
    if h.shape != m.shape:
        raise ValueError('Variance and effective-count arrays must match')
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.sum(h, axis=0)**2 / np.sum(h**2 / m, axis=0)


class TagLDAverages:
    """Within-bin tag means with exact leave-group-out recomputation."""
    def __init__(self, L2, bins, p):
        self.L2, self.bins, self.p = np.asarray(L2), np.asarray(bins), p
        self.n = np.bincount(bins, minlength=p)
        self.sums = np.array([np.bincount(bins, weights=L2[:, j], minlength=p)
                              for j in range(p)]).T
        if np.any(self.n == 0):
            raise ValueError('Each source bin needs representative tags')
        self.mean = self.sums / self.n[:, None]

    def delete(self, start, stop):
        b = self.bins[start:stop]; L = self.L2[start:stop]
        n = self.n - np.bincount(b, minlength=self.p)
        if np.any(n == 0):
            raise ValueError('Deletion empties a tag bin')
        s = np.array([np.bincount(b, weights=L[:, j], minlength=self.p)
                       for j in range(self.p)]).T
        return (self.sums - s) / n[:, None]
