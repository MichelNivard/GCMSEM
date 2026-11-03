#!/usr/bin/env python3
"""Reproduce the full 15-trait/order-four computational smoke test.

This is a dimension and execution check, not interval-coverage evidence. The
tiny independent-SNP fixture deliberately has fewer deletion groups than
moments; its sampling covariance must be singular. The separately supplied
GWASBrewer validation evaluates realistic LD and statistical behavior.

Run from the repository root:
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=python python validation/scalability.py
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from gcmsem import __version__
from gcmsem.moments import estimate_moments


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "validation" / "results" / "scalability.json")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    seed, n, p, blocks = 902, 120, 15, 6
    rng = np.random.default_rng(seed)
    beta = rng.normal(size=(n, p))
    scores = {
        **{f"L{q}": np.ones((n, 1)) for q in range(1, 5)},
        "counts": np.array([n]),
        "block": np.repeat(np.arange(blocks), n // blocks),
        "chrom": np.ones(n, dtype=int),
    }
    started = time.monotonic()
    result = estimate_moments(beta, np.ones_like(beta), scores,
                              noise_intercept=np.zeros((p, p)), n_blocks=blocks)
    elapsed = time.monotonic() - started
    expected_moments = 3875
    assert len(result["indices"]) == expected_moments
    assert result["sampling_cov"].shape == (expected_moments, expected_moments)
    assert result["metadata"]["sampling_covariance_rank"] == blocks - 1
    assert np.isfinite(result["estimates"]).all()
    assert np.isfinite(result["sampling_cov"]).all()
    np.testing.assert_allclose(result["sampling_cov"], result["sampling_cov"].T,
                               rtol=1e-12, atol=1e-12)
    summary = {
        "purpose": "Computational dimensionality smoke; not statistical or interval calibration",
        "status": "passed",
        "seed": seed,
        "n_traits": p,
        "n_tags": n,
        "jackknife_groups": blocks,
        "max_order": 4,
        "n_moments": expected_moments,
        "sampling_covariance_shape": list(result["sampling_cov"].shape),
        "sampling_covariance_rank": result["metadata"]["sampling_covariance_rank"],
        "elapsed_seconds": elapsed,
        "estimated_peak_memory_mib": result["metadata"]["estimated_peak_memory_mib"],
        "memory_measurement": "Conservative estimator guard calculation, not measured process RSS",
        "effect_fixture": "Independent Gaussian draws, identity LD, supplied zero sampling-noise matrix",
        "warnings": result["metadata"]["warnings"],
        "versions": {"gcmsem": __version__, "python": platform.python_version(),
                     "numpy": np.__version__, "scipy": scipy.__version__},
        "platform": platform.platform(),
        "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        "command": "OPENBLAS_NUM_THREADS=1 PYTHONPATH=python python validation/scalability.py",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    run()
