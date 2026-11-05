"""Independent identities for the tutorial's effective-polygenicity conversion."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location('poly_helpers', Path(__file__).parents[2] / 'tutorials/polygenicity_helpers.py')
p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)


def theta(v, k):
    v, k = np.atleast_1d(v), np.atleast_1d(k)
    return np.column_stack([np.zeros(len(v)), v, np.zeros(len(v)), k])


def test_gaussian_and_point_normal_normalizations():
    assert p.polygenicity_metrics(theta(.01, 0), [1000], [[1]])[2] == pytest.approx(1000)
    # Gaussian slab, 10% causal: v=pi*sigma2 and k4=3*pi*(1-pi)*sigma2^2.
    v = .1*.03; k = 3*.1*.9*.03**2
    assert p.polygenicity_metrics(theta(v, k), [1000], [[1]])[2] == pytest.approx(100)


def test_perfect_ld_reduces_uniform_gaussian_count_to_one():
    assert p.polygenicity_metrics(theta(.01, 0), [1000], [[1000]])[3] == pytest.approx(1)


def test_heterogeneous_strata_match_direct_variance_participation_ratio():
    counts = np.array([100, 200]); v = np.array([.01, .03])
    expected = (counts@v)**2 / (counts@(v*v))
    result = p.polygenicity_metrics(theta(v, [0, 0]), counts, np.eye(2))
    assert result[2] == pytest.approx(expected)
    assert result[3] == pytest.approx(expected)
    assert p.combine_sld4m_bins(counts*v, counts) == pytest.approx(expected)
    assert result[1] > 0  # Pooling Gaussian bins creates excess kurtosis.


def test_effect_scale_invariance():
    t = theta([.01, .03], [.002, .001]); n = [100, 200]; C = [[2, 3], [1.5, 4]]
    scaled = t*np.array([7, 7**2, 7**3, 7**4])
    a, b = p.polygenicity_metrics(t,n,C), p.polygenicity_metrics(scaled,n,C)
    np.testing.assert_allclose(a[1:4], b[1:4])


def test_signed_algebraic_estimates_are_retained():
    r = p.polygenicity_metrics(theta(.01, -.0004), [100], [[1]])
    assert r[2] < 0 and r[3] < 0


def test_tag_averages_recompute_after_deletion():
    L = np.array([[1.,2.], [3.,4.], [5.,6.], [7.,8.]])
    x = p.TagLDAverages(L, np.array([0,1,0,1]), 2)
    np.testing.assert_allclose(x.mean, [[3,4],[5,6]])
    np.testing.assert_allclose(x.delete(0,2), [[5,6],[7,8]])
    with pytest.raises(ValueError): x.delete(0,3)
