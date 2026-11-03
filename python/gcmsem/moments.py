"""LD-deconvolved cumulants for an arbitrary number of traits.

The working model is b_hat[j] = sum_k R[j,k] beta[k] + error[j].
Conditional source-variant cumulants are constant within disjoint annotations;
different source variants are independent. Sampling errors are jointly Gaussian,
independent of effects, with covariance se[j,a]*se[j,b]*intercept[a,b].

Means and sampling intercepts are fitted, not assumed zero. Higher-order
responses use centered multivariate Hermite polynomials. Fourth-order responses
also subtract the three marginal signal-covariance pairings. All stages are
refitted in each genomic deletion. This is a consistent plug-in procedure under
its working model, not a finite-sample unbiased cumulant estimator. Reference-LD
uncertainty and misspecified architectures are not covered by the jackknife.
"""
from __future__ import annotations

import itertools
import logging
import math
import re
import time
from functools import lru_cache

import numpy as np
from scipy.linalg import cho_factor, cho_solve


def moment_indices(n_traits, min_order=1, max_order=4):
    """Canonical, zero-based trait tuples, ordered by degree then lexically."""
    if not isinstance(n_traits, (int, np.integer)) or n_traits < 1:
        raise ValueError("n_traits must be a positive integer")
    if not 1 <= min_order <= max_order <= 4:
        raise ValueError("orders must satisfy 1 <= min_order <= max_order <= 4")
    return [t for q in range(min_order, max_order + 1)
            for t in itertools.combinations_with_replacement(range(n_traits), q)]


@lru_cache(None)
def partitions(n):
    """Set partitions of occurrence positions; duplicate indices retain multiplicity."""
    if n == 0:
        return ((),)
    answer = []
    for previous in partitions(n - 1):
        answer.append(previous + ((n - 1,),))
        for b in range(len(previous)):
            answer.append(previous[:b] + (previous[b] + (n - 1,),) + previous[b + 1:])
    return tuple(answer)


def raw_from_cumulants(t, cumulants, omit_full=False):
    """Multivariate cumulant-to-raw-moment identity for one canonical tuple."""
    total = 0.0
    for partition in partitions(len(t)):
        if omit_full and len(partition) == 1:
            continue
        term = 1.0
        for positions in partition:
            term = term * cumulants[tuple(t[k] for k in positions)]
        total = total + term
    return total


def _keys(n_traits, max_order, selected, max_moments):
    if not isinstance(max_order, (int, np.integer)) or not 2 <= max_order <= 4:
        raise ValueError("max_order must be 2, 3, or 4")
    if selected is None:
        count = math.comb(n_traits + max_order, max_order) - 1
        if count > max_moments:
            raise ValueError(f"{count:,} moments exceed max_moments={max_moments:,}; "
                             "select higher-order tuples or reduce max_order")
        return moment_indices(n_traits, max_order=max_order)
    keep = set(moment_indices(n_traits, max_order=2))
    for raw in selected:
        t = tuple(raw)
        if not 1 <= len(t) <= max_order or any(not isinstance(i, (int, np.integer))
                or i < 0 or i >= n_traits for i in t):
            raise ValueError("moment_tuples must contain valid zero-based trait tuples")
        t = tuple(sorted(t))
        # Pooling requires all lower-order subtuples, including triples of quartets.
        for q in range(1, len(t) + 1):
            keep.update(itertools.combinations(t, q))
    if len(keep) > max_moments:
        raise ValueError("selected tuples and required lower orders exceed max_moments")
    return sorted(keep, key=lambda t: (len(t), t))


def validate_traits(traits, p):
    if traits is None:
        traits = [f"T{i + 1}" for i in range(p)]
    traits = list(traits)
    if len(traits) != p or len(set(traits)) != p:
        raise ValueError("traits must be unique names matching beta columns")
    if any(not isinstance(t, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", t)
           for t in traits):
        raise ValueError("trait names must start with a letter and use letters, digits, _ or .")
    return traits


def _noise_matrix(value, p, name):
    value = np.asarray(value, dtype=float)
    if value.shape != (p, p) or not np.isfinite(value).all():
        raise ValueError(f"{name} must be a finite {p} by {p} matrix")
    if not np.allclose(value, value.T, rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} must be symmetric")
    eig = np.linalg.eigvalsh(value)
    if eig[0] < -1e-8 * max(float(np.max(np.abs(eig))), 1e-12):
        raise ValueError(f"{name} is not positive semidefinite (min eigenvalue {eig[0]:.6g}); "
                         "no automatic projection is performed. Check scale, LD, overlap, "
                         "or supply an independently justified noise_intercept.")
    return value


class _Stage:
    """Scaled weighted normal equations with cached genomic-block contributions."""
    def __init__(self, design, weights, groups, name, max_condition=1e6):
        self.name = name
        design = np.asarray(design, dtype=float)
        self.scale = np.sqrt(np.sum(weights[:, None] * design ** 2, axis=0))
        if np.any(self.scale == 0) or not np.isfinite(self.scale).all():
            raise ValueError(f"{name}: empty/invalid score columns; merge unsupported annotations")
        z = design / self.scale
        self.zw = z * weights[:, None]
        self.gram = z.T @ self.zw
        self.block_gram = [z[g].T @ self.zw[g] for g in groups]
        self.groups = groups
        self.factors = []
        self.conditions = []
        for i in range(-1, len(groups)):
            gram = self.gram if i < 0 else self.gram - self.block_gram[i]
            eig = np.linalg.eigvalsh(gram)
            cond = np.inf if eig[0] <= 0 else np.sqrt(eig[-1] / eig[0])
            if not np.isfinite(cond) or cond > max_condition:
                which = "full data" if i < 0 else f"deletion {i}"
                raise ValueError(f"{name}, {which}: rank deficient or ill-conditioned design "
                                 f"(condition={cond:.4g}); merge annotations or supply noise intercepts")
            self.conditions.append(float(cond))
            self.factors.append(cho_factor(gram, lower=True, check_finite=False))

    def solve(self, response, deletion=None):
        rhs = self.zw.T @ response
        if deletion is not None:
            g = self.groups[deletion]
            rhs -= self.zw[g].T @ response[g]
        ans = cho_solve(self.factors[0 if deletion is None else deletion + 1], rhs,
                        check_finite=False)
        return ans / self.scale if ans.ndim == 1 else ans / self.scale[:, None]


class _OverlapStage:
    """Append an SE-product predictor without duplicating n-by-C LD arrays."""
    def __init__(self, base, se, pair, weights, groups, name, max_condition):
        self.base, self.se, self.pair = base, se, pair
        self.weights, self.groups = weights, groups
        a, b = pair
        product = se[:, a] * se[:, b]
        self.noise_scale = np.sqrt(np.sum(weights * product ** 2))
        znoise = product / self.noise_scale
        cross = base.zw.T @ znoise
        gram = np.block([[base.gram, cross[:, None]], [cross[None, :], np.ones((1, 1))]])
        self.factors, self.conditions = [], []
        for i in range(-1, len(groups)):
            g = gram.copy()
            if i >= 0:
                s = groups[i]
                cross_block = base.zw[s].T @ znoise[s]
                g[:-1, :-1] -= base.block_gram[i]
                g[:-1, -1] -= cross_block
                g[-1, :-1] -= cross_block
                g[-1, -1] -= np.sum(weights[s] * znoise[s] ** 2)
            eig = np.linalg.eigvalsh(g)
            cond = np.inf if eig[0] <= 0 else np.sqrt(eig[-1] / eig[0])
            if not np.isfinite(cond) or cond > max_condition:
                which = "full data" if i < 0 else f"deletion {i}"
                raise ValueError(f"{name}, {which}: rank deficient or ill-conditioned design "
                                 f"(condition={cond:.4g}); overlap and genetic covariance cannot be separated")
            self.conditions.append(float(cond))
            self.factors.append(cho_factor(g, lower=True, check_finite=False))

    def solve(self, response, deletion=None):
        a, b = self.pair
        wz = self.weights * self.se[:, a] * self.se[:, b] / self.noise_scale
        rhs = np.r_[self.base.zw.T @ response, wz @ response]
        if deletion is not None:
            s = self.groups[deletion]
            rhs[:-1] -= self.base.zw[s].T @ response[s]
            rhs[-1] -= wz[s] @ response[s]
        ans = cho_solve(self.factors[0 if deletion is None else deletion + 1], rhs,
                        check_finite=False)
        return ans / np.r_[self.base.scale, self.noise_scale]


def make_jackknife_groups(block, chrom, n_groups):
    """Aggregate consecutive whole LD blocks without crossing chromosomes.

    Groups aim for equal tag-row counts. Exact equality is generally impossible
    because an LD block must never be split. Returned fractions are used in the
    unequal-delete-group jackknife rather than silently assuming equal sizes.
    """
    block, chrom = np.asarray(block).astype(str), np.asarray(chrom).astype(str)
    n = len(block)
    if chrom.shape != (n,) or n < 2:
        raise ValueError("block and chrom must have one entry per tag SNP")
    boundaries = np.r_[0, 1 + np.flatnonzero(block[1:] != block[:-1]), n]
    ids = [block[a] for a in boundaries[:-1]]
    if len(ids) != len(set(ids)):
        raise ValueError("each LD block must form one contiguous group of rows")
    whole = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        if len(set(chrom[a:b])) != 1:
            raise ValueError("an LD block cannot span chromosomes")
        whole.append((str(chrom[a]), int(a), int(b)))
    chroms = list(dict.fromkeys(chrom.tolist()))
    if len(chroms) != len(list(dict.fromkeys(c for c, _, _ in whole))):
        raise ValueError("chromosome rows must be ordered and contiguous")
    for c in chroms:
        ix = np.flatnonzero(chrom == c)
        if ix[-1] - ix[0] + 1 != len(ix):
            raise ValueError("chromosome rows must be contiguous")
    if not isinstance(n_groups, (int, np.integer)) or n_groups < 2:
        raise ValueError("n_blocks must be at least 2")
    target = min(n_groups, len(whole))
    if target < max(2, len(chroms)):
        raise ValueError("need at least two independent LD blocks and at least one group per chromosome")
    available = np.array([sum(bc == c for bc, _, _ in whole) for c in chroms])
    size = np.array([sum(b - a for bc, a, b in whole if bc == c) for c in chroms])
    allocation = np.ones(len(chroms), dtype=int)
    while allocation.sum() < target:
        deficit = target * size / n - allocation
        deficit[allocation >= available] = -np.inf
        allocation[int(np.argmax(deficit))] += 1
    groups, description = [], []
    for c, count in zip(chroms, allocation):
        entries = [(a, b) for bc, a, b in whole if bc == c]
        # Greedy closest boundary to each target, retaining at least one block
        # for every remaining group. This cannot split an LD-dependence block.
        start = 0
        cum = np.cumsum([b - a for a, b in entries])
        for k in range(int(count)):
            if k == count - 1:
                stop = len(entries)
            else:
                candidates = np.arange(start + 1, len(entries) - (count - k - 1) + 1)
                goal = (k + 1) * cum[-1] / count
                stop = int(candidates[np.argmin(abs(cum[candidates - 1] - goal))])
            lo, hi = entries[start][0], entries[stop - 1][1]
            groups.append(slice(lo, hi))
            description.append({"chrom": c, "start_row": lo, "stop_row": hi,
                                "n_tags": hi - lo, "n_ld_blocks": stop - start})
            start = stop
    fractions = np.array([g.stop - g.start for g in groups], dtype=float) / n
    return groups, fractions, description


def jackknife_cov(replicates, point=None, fractions=None):
    """Delete-group covariance, supporting unequal whole-LD-block groups.

    For equal sizes this is (B-1)/B times centered deletion cross-products.
    Unequal sizes use weighted pseudovalues (Busing/Patterson formulation):
    h_b=1/f_b, tau_b=h_b*theta-(h_b-1)*theta_-b, theta_J=sum f_b*tau_b,
    V=(1/B)*sum (tau_b-theta_J)(tau_b-theta_J)'/(h_b-1).
    Unequal-block and nonlinear-estimator calibration remains approximate.
    """
    values = np.asarray(replicates, dtype=float)
    if values.ndim != 2 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("jackknife replicates must be a finite B by M matrix, B>=2")
    b = len(values)
    if fractions is None:
        delta = values - values.mean(axis=0)
        return (b - 1) / b * delta.T @ delta
    fractions = np.asarray(fractions, dtype=float)
    if fractions.shape != (b,) or np.any(fractions <= 0) or np.any(fractions >= 1) \
            or not np.isclose(fractions.sum(), 1):
        raise ValueError("jackknife fractions must be positive, sum to one, and be below one")
    if np.allclose(fractions, 1 / b, atol=1e-12, rtol=1e-10):
        return jackknife_cov(values)
    point = np.asarray(point, dtype=float)
    if point.shape != (values.shape[1],) or not np.isfinite(point).all():
        raise ValueError("unequal jackknife requires the finite full-sample point estimate")
    pseudo = point + (1 / fractions[:, None] - 1) * (point - values)
    mean = fractions @ pseudo
    delta = pseudo - mean
    return delta.T @ (delta * (fractions / (1 - fractions))[:, None]) / b


def pooled_cumulants(conditional, weights, indices):
    """Mix annotation distributions through raw moments, then recenter."""
    output = {}
    for t in indices:
        output[t] = float(weights @ raw_from_cumulants(t, conditional)) \
                    - raw_from_cumulants(t, output, omit_full=True)
    return np.array([output[t] for t in indices])


def _hermite_response(t, x, se, intercept):
    q = len(t)
    response = np.prod(x[:, t], axis=1)
    if q >= 3:
        for a, b in itertools.combinations(range(q), 2):
            rest = [t[i] for i in range(q) if i != a and i != b]
            response -= se[:, t[a]] * se[:, t[b]] * intercept[t[a], t[b]] \
                        * np.prod(x[:, rest], axis=1)
    if q == 4:
        for pairs in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
            term = np.ones(len(x))
            for a, b in pairs:
                term *= se[:, t[a]] * se[:, t[b]] * intercept[t[a], t[b]]
            response += term
    return response


def estimate_moments(beta, se, scores, traits=None, max_order=4, moment_tuples=None,
                     noise_intercept=None, n_blocks=200, max_moments=5000,
                     max_memory_mb=2048, max_condition=1e6, logger=None):
    """Estimate pooled and annotation-specific cumulants and joint uncertainty.

    Parameters
    ----------
    beta, se : ndarray, shape (n_tags, n_traits)
        Aligned genotype-standardized effects and standard errors. Allele
        orientation must match the score reference. No missing values allowed.
    scores : mapping
        L1..L4 (n_tags by n_annotations), counts (reference/source counts),
        block and chrom tag arrays. API assumes these already passed score QC.
    moment_tuples : iterable of integer tuples, optional
        Higher-order tuples to estimate. Means, all covariances and necessary
        lower-order subtuples are retained. None means every tuple through
        max_order; memory guards reject infeasible requests explicitly.
    noise_intercept : ndarray, optional
        Known sampling-error calibration matrix; defaults to estimation jointly
        with L2 slopes. A supplied matrix is held fixed in every deletion.

    Returns
    -------
    dict
        estimates, sampling_cov, jackknife, conditional,
        conditional_jackknife, noise_intercept, intercept_jackknife,
        indices, traits, metadata. Fourth-order entries are excess cumulants.
    """
    started = time.monotonic()
    log = logger or logging.getLogger("gcmsem")
    if not isinstance(max_moments, (int, np.integer)) or max_moments < 1:
        raise ValueError("max_moments must be a positive integer")
    if not np.isfinite(max_memory_mb) or max_memory_mb <= 0:
        raise ValueError("working memory limit max_memory_mb must be positive and finite")
    if not np.isfinite(max_condition) or max_condition < 1:
        raise ValueError("max_condition must be finite and at least one")
    y, se = np.asarray(beta, dtype=float), np.asarray(se, dtype=float)
    if y.ndim != 2 or y.shape != se.shape or min(y.shape) < 1:
        raise ValueError("beta and se must have matching nonempty n_tags by n_traits shapes")
    if not np.isfinite(y).all() or not np.isfinite(se).all() or np.any(se <= 0):
        raise ValueError("beta and se must be finite; standard errors must be strictly positive")
    n, p = y.shape
    traits = validate_traits(traits, p)
    indices = _keys(p, max_order, moment_tuples, max_moments)
    m = len(indices)
    counts = np.asarray(scores["counts"], dtype=float)
    if counts.ndim != 1 or not len(counts) or not np.isfinite(counts).all() or np.any(counts <= 0):
        raise ValueError("scores counts must be positive finite reference/source annotation counts")
    c = len(counts)
    L = {}
    for q in range(1, max_order + 1):
        L[q] = np.asarray(scores[f"L{q}"], dtype=float)
        if L[q].shape != (n, c) or not np.isfinite(L[q]).all():
            raise ValueError(f"L{q} must be a finite n_tags by n_annotations matrix")
    groups, fractions, descriptions = make_jackknife_groups(scores["block"], scores["chrom"], n_blocks)
    bcount = len(groups)
    # Estimate peak working memory, including full joint covariance. Explicit
    # selected moments avoid fourth-order combinatorial allocations.
    memory = 8 * (6 * n * p + 8 * n * c + 2 * m * m
                  + bcount * m * (c + 3) + 2 * bcount * p * p)
    if memory > max_memory_mb * 1024 ** 2:
        raise ValueError(f"estimated working memory {memory / 1024 ** 2:.0f} MiB exceeds "
                         f"max_memory_mb={max_memory_mb}; select tuples or raise the explicit limit")
    if n <= c + 2:
        raise ValueError("too few tag SNPs for annotation and overlap regressions")
    log.info("ESTIMATION tags=%d traits=%d annotations=%d moments=%d max_order=%d", n, p, c, m, max_order)
    log.info("JACKKNIFE groups=%d fraction_min=%.5g fraction_max=%.5g effective_groups=%.2f",
             bcount, fractions.min(), fractions.max(), 1 / np.sum(fractions ** 2))
    warnings = []
    if bcount < 30:
        warnings.append("Fewer than 30 genomic jackknife groups; uncertainty is weakly calibrated")
    if fractions.max() > 0.1:
        warnings.append("A deletion contains more than 10% of tag SNPs; unequal-block jackknife is approximate")
    if bcount <= m:
        warnings.append("Joint sampling covariance is necessarily rank deficient (moments >= jackknife groups); "
                        "DWLS sandwich inference does not invert this matrix")
    for warning in warnings:
        log.warning(warning)
    weights = 1 / np.maximum(L[2].sum(axis=1), 1.0)
    stages = {q: _Stage(L[q], weights, groups, f"order {q}", max_condition)
              for q in range(1, max_order + 1)}
    pairs = moment_indices(p, 2, 2)
    fixed_noise = None if noise_intercept is None else _noise_matrix(noise_intercept, p, "noise_intercept")
    second = {} if fixed_noise is not None else {
        t: _OverlapStage(stages[2], se, t, weights, groups,
                         f"overlap {traits[t[0]]}/{traits[t[1]]}", max_condition) for t in pairs}

    def fit(deletion=None):
        theta = {}
        means = stages[1].solve(y, deletion)
        x = y - L[1] @ means
        for a in range(p):
            theta[(a,)] = means[:, a]
        intercept = np.zeros((p, p)) if fixed_noise is None else fixed_noise.copy()
        for a, b in pairs:
            response = x[:, a] * x[:, b]
            if fixed_noise is None:
                solution = second[(a, b)].solve(response, deletion)
                theta[(a, b)] = solution[:-1]
                intercept[a, b] = intercept[b, a] = solution[-1]
            else:
                response -= se[:, a] * se[:, b] * intercept[a, b]
                theta[(a, b)] = stages[2].solve(response, deletion)
        _noise_matrix(intercept, p, f"estimated noise intercept (deletion={deletion})")
        # Only coefficient vectors are retained across moments, avoiding an
        # n_tags by number-of-fourth-order-moments array.
        for t in indices:
            if len(t) < 3:
                continue
            response = _hermite_response(t, x, se, intercept)
            if len(t) == 4:
                for pairings in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
                    term = np.ones(n)
                    for a, b in pairings:
                        term *= L[2] @ theta[tuple(sorted((t[a], t[b])))]
                    response -= term
            theta[t] = stages[len(t)].solve(response, deletion)
        pooled = pooled_cumulants(theta, counts / counts.sum(), indices)
        conditional = np.column_stack([theta[t] for t in indices])
        if not np.isfinite(pooled).all() or not np.isfinite(conditional).all():
            raise FloatingPointError("non-finite cumulant estimate; inspect input scales")
        return pooled, conditional, intercept

    point, conditional, intercept = fit()
    log.info("NOISE_INTERCEPT %s", np.array2string(intercept, precision=6, separator=","))
    replicates = np.empty((bcount, m))
    conditional_jk = np.empty((bcount, c, m))
    intercept_jk = np.empty((bcount, p, p))
    for b in range(bcount):
        try:
            replicates[b], conditional_jk[b], intercept_jk[b] = fit(b)
        except Exception as exc:
            raise ValueError(f"jackknife deletion {b} failed; no partial covariance will be reported: {exc}") from exc
        if (b + 1) % max(1, bcount // 10) == 0 or b + 1 == bcount:
            log.info("JACKKNIFE completed=%d/%d elapsed_seconds=%.2f", b + 1, bcount, time.monotonic() - started)
    covariance = jackknife_cov(replicates, point, fractions)
    variance_keys = {t: i for i, t in enumerate(indices)}
    variances = np.array([point[variance_keys[(a, a)]] for a in range(p)])
    if np.any(variances <= 0):
        warnings.append("One or more estimated SNP-effect variances are nonpositive; corresponding correlations are undefined")
        log.warning(warnings[-1])
    genetic_cov = np.array([[point[variance_keys[tuple(sorted((a, b)))]] for b in range(p)] for a in range(p)])
    covariance_min_eigenvalue = float(np.linalg.eigvalsh(genetic_cov)[0])
    if np.all(variances > 0):
        genetic_corr = genetic_cov / np.sqrt(np.outer(variances, variances))
        if np.linalg.eigvalsh(genetic_corr)[0] < -1e-8:
            warnings.append("Estimated genetic covariance is not positive semidefinite; estimates are retained without projection")
            log.warning(warnings[-1])
    diagonal = np.diag(covariance)
    active = diagonal > 0
    rank = 0
    if active.any():
        scale = np.sqrt(diagonal[active])
        # The B-by-M factor is much smaller than a full M-by-M correlation
        # eigendecomposition when thousands of moments share ~200 deletions.
        pseudo = point + (1 / fractions[:, None] - 1) * (point - replicates)
        centered = pseudo - fractions @ pseudo
        factor = centered[:, active] * np.sqrt(fractions / (1 - fractions) / bcount)[:, None] / scale
        rank = int(np.linalg.matrix_rank(factor, tol=1e-4))
    metadata = {
        "schema_version": "1.0", "traits": traits,
        "effect_scale": "genotype_standardized", "moment_type": "cumulant",
        "n_tags": n, "n_reference_variants": float(counts.sum()),
        "n_annotations": c, "source_counts": counts.tolist(),
        "mixture_weights": (counts / counts.sum()).tolist(),
        "moment_indices": [list(t) for t in indices], "n_moments": m,
        "max_order": max_order, "noise_mode": "estimated" if fixed_noise is None else "supplied_fixed",
        "noise_intercept": intercept.tolist(),
        "noise_intercept_se": np.sqrt(np.maximum(np.diag(jackknife_cov(
            intercept_jk.reshape(bcount, -1), intercept.reshape(-1), fractions)), 0)).reshape(p, p).tolist(),
        "noise_min_eigenvalue_all_fits": min(float(np.linalg.eigvalsh(a)[0]) for a in [intercept, *intercept_jk]),
        "jackknife_blocks": descriptions, "jackknife_fractions": fractions.tolist(),
        "jackknife_method": "equal-delete-group" if np.allclose(fractions, 1 / bcount)
                            else "unequal-delete-group weighted pseudovalues (Busing/Patterson)",
        "sampling_covariance_rank": rank,
        "genetic_covariance_min_eigenvalue": covariance_min_eigenvalue,
        "regression_weight": "1 / max(sum_c L2_jc, 1), fixed by reference",
        "design_condition_max": {str(q): max(stages[q].conditions) for q in stages},
        "intercept_design_condition_max": None if fixed_noise is not None
                    else max(max(stage.conditions) for stage in second.values()),
        "estimated_peak_memory_mib": memory / 1024 ** 2,
        "elapsed_seconds": time.monotonic() - started,
        "warnings": warnings,
        "assumptions": ["Independent source-variant effects conditional on disjoint annotations",
            "Gaussian sampling errors independent of effects; supplied or LD-regression-estimated covariance",
            "Reference LD, allele orientation, genotype scaling and source universe are appropriate",
            "Whole LD blocks are approximately independent across jackknife groups"],
        "limitations": ["Finite-reference LD uncertainty is not propagated; inspect score metadata for power debiasing",
            "Plug-in higher moments are not exactly unbiased in finite samples",
            "Unequal-block nonlinear jackknife and higher-order interval coverage require calibration",
            "Pooled annotation mixtures may induce dependence between otherwise independent SEM sources",
            "No liability transformation, binary-GWAS ascertainment correction, or automatic causal identification"],
    }
    log.info("DONE moments=%d sampling_cov_rank=%d/%d elapsed_seconds=%.2f", m, rank, m, metadata["elapsed_seconds"])
    return {"indices": indices, "traits": traits, "estimates": point,
            "sampling_cov": covariance, "jackknife": replicates,
            "conditional": conditional, "conditional_jackknife": conditional_jk,
            "noise_intercept": intercept, "intercept_jackknife": intercept_jk,
            "metadata": metadata}
