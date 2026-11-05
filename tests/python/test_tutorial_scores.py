"""Independent tiny dense references for the empirical score rebuild scripts."""
import importlib.util
from pathlib import Path
import sys

import numpy as np


TUTORIALS = Path(__file__).resolve().parents[2] / "tutorials"
spec = importlib.util.spec_from_file_location("build_ldsc_weights", TUTORIALS / "build_ldsc_weights.py")
weights = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = weights
spec.loader.exec_module(weights)
spec = importlib.util.spec_from_file_location("rebuild_empirical_scores", TUTORIALS / "rebuild_empirical_scores.py")
moments = importlib.util.module_from_spec(spec)
spec.loader.exec_module(moments)


def reference():
    rng = np.random.default_rng(5115)
    x = rng.normal(size=(31, 503))
    x -= x.mean(axis=1, keepdims=True)
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    return x.astype(np.float32), np.arange(31) * 100000, np.arange(31) % 5


def dense_oracle(x, bp, bins, rows, order):
    result = np.zeros((len(rows), 5))
    for i, tag in enumerate(rows):
        for source in range(len(x)):
            if abs(int(bp[tag]) - int(bp[source])) > 400000:
                continue
            r = 1. if tag == source else float(x[tag].astype(float) @ x[source].astype(float))
            value = r ** order
            if order == 2:
                value -= (1 - r ** 2) / 501
            result[i, bins[source]] += value
    return result


def test_hm3_weights_match_pairwise_finite_reference_oracle():
    x, bp, bins = reference()
    actual, _, _ = weights.score_chromosome(x, bp, bins, 503, 400000, 7)
    expected = dense_oracle(x, bp, bins, np.arange(len(x)), 2)
    np.testing.assert_allclose(actual, expected, atol=1e-7, rtol=1e-7)
    # Unbiased r-squared terms may be negative and must not be truncated.
    assert np.any(actual < 0)


def test_rectangular_scores_match_signed_pairwise_oracle():
    x, bp, bins = reference()
    rows = np.array([0, 2, 4, 7, 9, 13, 25, 30])
    actual, _ = moments.score_rect(x, bp, bins, rows, 503, 400000, 7)
    for order in range(1, 5):
        np.testing.assert_allclose(actual[order - 1], dense_oracle(x, bp, bins, rows, order),
                                   atol=1e-7, rtol=1e-7)
    assert np.any(actual[0] < 0)
    assert np.any(actual[2] < 0)
