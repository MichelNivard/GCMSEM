#!/usr/bin/env python3
"""Independent fourth-cumulant oracle check using the generated real-LD blocks.

This isolates the LD-power identity from GWAS sampling error and fitted nuisance
quantities. It does not validate the complete estimator or its confidence
intervals. The source distribution and B = 0.4*A + e_B match the chain scenario
in generate_gwasbrewer.R. Run after that generator and run_moments.py.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import mmread


def run(workdir, replicates, seed):
    if replicates < 2:
        raise ValueError("At least two independent replications are required")
    with np.load(workdir / "scores.npz", allow_pickle=False) as scores:
        if scores["L2"].shape[1] != 1:
            raise ValueError("This scenario check requires one source annotation")
        l2 = scores["L2"][:, 0]
        l4 = scores["L4"][:, 0]
    with (workdir / "ld/blocks.tsv").open() as handle:
        blocks = list(csv.DictReader(handle, delimiter="\t"))
    pattern_size = sum(int(row["size"]) for row in blocks)
    J = len(l2)
    if J != 2 * pattern_size:
        raise ValueError("Expected two independent copies of the generating LD pattern")
    matrices = [mmread(workdir / f"ld/block_{int(row['block']):02d}.mtx").tocsr()
                for row in blocks]
    variance_a, variance_residual = 0.2 / J, 0.168 / J
    path, pi = 0.4, 0.15
    variance_b = path**2 * variance_a + variance_residual
    truth = (path**4 * variance_a**2 + variance_residual**2) * (6 / pi - 3)
    weights = 1 / np.maximum(l2, 1)
    u = weights * l4
    denominator = u @ l4
    rng = np.random.default_rng(seed)
    numerator = np.zeros(replicates)
    source_fourth = np.zeros(replicates)
    start = time.monotonic()
    offset = 0
    for matrix in matrices * 2:
        n = matrix.shape[0]
        # Each source is a centered, unit-variance Gamma(2) conditional on being
        # nonzero. Bernoulli selection retains the zero-effect source variants.
        a = ((rng.gamma(2, size=(n, replicates)) - 2) / np.sqrt(2)
             * (rng.random((n, replicates)) < pi) * np.sqrt(variance_a / pi))
        b = ((rng.gamma(2, size=(n, replicates)) - 2) / np.sqrt(2)
             * (rng.random((n, replicates)) < pi) * np.sqrt(variance_residual / pi))
        beta = path * a + b
        marginal = matrix @ beta
        numerator += u[offset:offset + n] @ marginal**4
        source_fourth += np.sum(beta**4, axis=0)
        offset += n
    estimates = (numerator - 3 * (u @ (l2 * variance_b)**2)) / denominator
    # The zero mean and variance are population quantities in this oracle.
    # This source check is likewise E(beta^4)-3*Var(beta)^2, without fitted
    # centering or a finite-sample empirical-cumulant correction.
    source_cumulants = source_fourth / J - 3 * variance_b**2
    sd = float(estimates.std(ddof=1))
    probabilities = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    summary = {
        "experiment": "noiseless scalar B fourth-cumulant LD-regression oracle",
        "replicates": replicates,
        "seed": seed,
        "source_variants": J,
        "independent_ld_blocks": len(matrices) * 2,
        "source_distribution": "15% nonzero centered Gamma(shape=2), 85% zero",
        "path_A_to_B": path,
        "known_nuisance": "population source mean zero and variance 0.2/J; no GWAS noise",
        "target": truth,
        "mean_estimate": float(estimates.mean()),
        "empirical_sd": sd,
        "monte_carlo_se_mean": sd / np.sqrt(replicates),
        "mean_to_target": float(estimates.mean() / truth),
        "sd_to_target": sd / truth,
        "monte_carlo_se_to_target": sd / (np.sqrt(replicates) * truth),
        "median_to_target": float(np.median(estimates) / truth),
        "quantile_probabilities": probabilities,
        "quantile_to_target": (np.quantile(estimates, probabilities) / truth).tolist(),
        "mean_realized_source_oracle_to_target": float(source_cumulants.mean() / truth),
        "elapsed_seconds": time.monotonic() - start,
        "scope": "Independent fresh source draws and known generating LD; checks a population moment identity, not the complete fitted-nuisance pipeline or interval coverage",
        "numpy_version": np.__version__,
    }
    return estimates, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=Path("validation/work"))
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=276351)
    parser.add_argument("--out", type=Path, default=Path("validation/results/fourth_order_oracle"))
    args = parser.parse_args()
    estimates, summary = run(args.workdir, args.replicates, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with Path(str(args.out) + ".tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["replicate", "estimate"])
        writer.writerows(enumerate(estimates, 1))
    Path(str(args.out) + ".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
