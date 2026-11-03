import itertools
import logging

import numpy as np
import pytest

from gcmsem.moments import (
    _hermite_response, estimate_moments, jackknife_cov, make_jackknife_groups,
    moment_indices, pooled_cumulants, raw_from_cumulants,
)


def independence_scores(n, blocks=20):
    return {**{f"L{q}": np.ones((n, 1)) for q in range(1, 5)},
            "counts": np.array([n]), "block": np.repeat(np.arange(blocks), n // blocks),
            "chrom": np.ones(n, dtype=int)}


def empirical_cumulants(y):
    mean = y.mean(axis=0)
    centered = y - mean
    covariance = centered.T @ centered / len(y)
    out = []
    for t in moment_indices(y.shape[1]):
        if len(t) == 1:
            answer = mean[t[0]]
        else:
            answer = np.prod(centered[:, t], axis=1).mean()
            if len(t) == 4:
                for (a, b), (c, d) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
                    answer -= covariance[t[a], t[b]] * covariance[t[c], t[d]]
        out.append(answer)
    return np.array(out)


def test_exact_three_trait_empirical_cumulants_and_full_refits():
    rng = np.random.default_rng(7701)
    source = np.column_stack((rng.gamma(2, 0.2, 800), rng.standard_t(7, 800), rng.normal(size=800)))
    y = source @ np.array([[1, .4, -.2], [0, 1, .3], [0, 0, 1]])
    scores = independence_scores(800)
    result = estimate_moments(y, np.ones_like(y), scores, noise_intercept=np.zeros((3, 3)), n_blocks=20)
    assert len(result["indices"]) == 34
    np.testing.assert_allclose(result["estimates"], empirical_cumulants(y), rtol=1e-11, atol=1e-11)
    for i in (0, 8, 19):
        np.testing.assert_allclose(result["jackknife"][i], empirical_cumulants(y[scores["block"] != i]), rtol=1e-10, atol=1e-10)
    assert result["sampling_cov"].shape == (34, 34)
    assert np.linalg.eigvalsh(result["sampling_cov"]).min() > -1e-10
    assert result["metadata"]["sampling_covariance_rank"] <= 19


def test_hermite_multiplicities_and_mixed_fourth():
    x = np.array([[.3, -.7], [1.5, .9]])
    se = np.array([[.2, .4], [.3, .1]])
    I = np.array([[1.2, -.4], [-.4, .9]])
    v, w, c = se[:, 0] ** 2 * I[0, 0], se[:, 1] ** 2 * I[1, 1], se[:, 0] * se[:, 1] * I[0, 1]
    a, b = x.T
    np.testing.assert_allclose(_hermite_response((0, 0, 0), x, se, I), a**3 - 3*v*a)
    np.testing.assert_allclose(_hermite_response((0, 0, 1), x, se, I), a*a*b-v*b-2*c*a)
    np.testing.assert_allclose(_hermite_response((0, 0, 0, 0), x, se, I), a**4-6*v*a*a+3*v*v)
    np.testing.assert_allclose(_hermite_response((0, 0, 1, 1), x, se, I), a*a*b*b-v*b*b-w*a*a-4*c*a*b+v*w+2*c*c)


def test_mixture_pooling_uses_raw_moments():
    # Two point masses: their within-bin variance/kurtosis vanish, but mixture
    # has positive variance and negative fourth cumulant.
    indices = moment_indices(1)
    conditional = {(0,): np.array([-1., 2.]), (0, 0): np.zeros(2),
                   (0, 0, 0): np.zeros(2), (0, 0, 0, 0): np.zeros(2)}
    actual = pooled_cumulants(conditional, np.array([.25, .75]), indices)
    expected = empirical_cumulants(np.array([[-1.], [2.], [2.], [2.]]))
    np.testing.assert_allclose(actual, expected)


def test_unequal_jackknife_formula_and_equal_reduction():
    rng = np.random.default_rng(11)
    reps = rng.normal(size=(5, 3)); point = rng.normal(size=3)
    equal = np.full(5, .2)
    np.testing.assert_allclose(jackknife_cov(reps, point, equal), .8*(reps-reps.mean(0)).T@(reps-reps.mean(0)))
    f = np.array([.1, .2, .15, .25, .3]); h = 1/f
    pseudo = h[:, None]*point-(h-1)[:, None]*reps
    center = f@pseudo
    expected = sum(np.outer(t-center, t-center)/(hh-1) for t,hh in zip(pseudo,h))/5
    np.testing.assert_allclose(jackknife_cov(reps, point, f), expected)


def test_whole_ld_blocks_never_split_or_cross_chromosome():
    lengths = [2, 7, 3, 10, 5, 9]
    block = np.repeat(np.arange(6), lengths)
    chrom = np.repeat([1, 1, 1, 2, 2, 2], lengths)
    groups, fractions, info = make_jackknife_groups(block, chrom, 4)
    assert len(groups) == 4
    assert sum(g.stop-g.start for g in groups) == len(block)
    for g in groups:
        assert len(set(chrom[g])) == 1
        assert g.start == 0 or block[g.start] != block[g.start-1]
        assert g.stop == len(block) or block[g.stop] != block[g.stop-1]
    assert np.isclose(fractions.sum(),1)
    assert len(info) == 4


def test_fifteen_traits_second_order_and_selected_quartet():
    rng = np.random.default_rng(88)
    y = rng.normal(size=(400, 15))
    score = independence_scores(400)
    result = estimate_moments(y, np.ones_like(y), score, max_order=2,
                             noise_intercept=np.zeros((15,15)), n_blocks=20)
    assert len(result["indices"]) == 135
    assert result["sampling_cov"].shape == (135,135)
    selected = estimate_moments(y, np.ones_like(y), score, moment_tuples=[(0,1,2,3)],
                               noise_intercept=np.zeros((15,15)), n_blocks=20)
    assert (0,1,2,3) in selected["indices"]
    assert (1,2,3) in selected["indices"]
    assert (0,0,0,0) not in selected["indices"]


def test_rank_and_noise_failures_are_explicit():
    y = np.random.default_rng(1).normal(size=(100,2))
    scores = independence_scores(100)
    with pytest.raises(ValueError, match="rank deficient"):
        estimate_moments(y,np.ones_like(y),scores,n_blocks=20)
    with pytest.raises(ValueError, match="positive semidefinite"):
        estimate_moments(y,np.ones_like(y),scores,noise_intercept=[[1,2],[2,1]],n_blocks=20)
    with pytest.raises(ValueError, match="strictly positive"):
        estimate_moments(y,np.zeros_like(y),scores,n_blocks=20)
    with pytest.raises(ValueError, match="max_moments"):
        estimate_moments(y,np.ones_like(y),scores,max_moments=2)
    with pytest.raises(ValueError, match="working memory"):
        estimate_moments(y,np.ones_like(y),scores,max_memory_mb=0)


def test_noise_overlap_is_estimated_in_joint_regression():
    # Exact moment-design fixture: paired plus/minus rows force known mean and
    # covariance at every heterogeneous L2/SE design point. This validates
    # coefficient separation without relying on a random convergence tolerance.
    cov = np.array([[.4,.12],[.12,.3]])
    I = np.array([[1.1,.35],[.35,.8]])
    records=[]; errors=[]; l2=[]; blocks=[]
    for block in range(20):
        for j in range(6):
            ld=1+.2*j; s=np.array([.1+.03*((j+block)%5), .12+.02*((j+2*block)%6)])
            total=ld*cov + s[:,None]*s[None,:]*I
            root=np.linalg.cholesky(total)
            # The four signed sqrt(2)*Cholesky columns have covariance total.
            for k in range(2):
                for sign in (-1,1):
                    records.append(sign*np.sqrt(2)*root[:,k]);errors.append(s);l2.append(ld);blocks.append(block)
    y=np.array(records); se=np.array(errors); n=len(y)
    scores=independence_scores(n,20);scores['L2']=np.array(l2)[:,None];scores['block']=np.array(blocks)
    result=estimate_moments(y,se,scores,max_order=2,n_blocks=20)
    np.testing.assert_allclose(result['noise_intercept'],I,atol=1e-10)
    lookup=dict(zip(result['indices'],result['estimates']))
    np.testing.assert_allclose([lookup[(0,0)],lookup[(0,1)],lookup[(1,1)]],[.4,.12,.3],atol=1e-11)
    np.testing.assert_allclose(result['intercept_jackknife'],np.broadcast_to(I,(20,2,2)),atol=1e-9)
