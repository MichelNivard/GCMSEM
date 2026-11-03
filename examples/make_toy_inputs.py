#!/usr/bin/env python3
"""Small deterministic CLI example; real-LD GWASBrewer tests are in validation/.

Writes three per-allele GWAS TSVs and signed generating-LD block NPZs.
No external data, downloads, or participant information are used.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="toy")
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "ld").mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)
    sizes = rng.integers(25, 71, size=80)
    n = int(sizes.sum())
    beta_parts, se_parts, af_parts, snp_parts = [], [], [], []
    ld_rows = []
    B = np.array([[0., 0., 0.], [.4, 0., 0.], [0., .25, 0.]])
    transfer = np.linalg.inv(np.eye(3) - B)
    overlap = np.full((3, 3), .1)
    np.fill_diagonal(overlap, 1.)
    noise_factor = np.linalg.cholesky(overlap)
    offset = 0
    for b, size in enumerate(sizes):
        distance = np.abs(np.arange(size)[:, None] - np.arange(size))
        signs = rng.choice([-1., 1.], size=size)
        ld = rng.uniform(.25, .75) ** distance * np.outer(signs, signs)
        af = rng.uniform(.08, .49, size=size)
        snp = np.array([f"toy{j+1}" for j in range(offset, offset + size)])
        source = ((rng.gamma(2., 1., (size, 3)) - 2.) / np.sqrt(2.)
                  * (rng.random((size, 3)) < .3) * np.sqrt(.15 / (.3 * n)))
        joint = source @ transfer.T
        noise = np.linalg.cholesky(ld) @ rng.normal(size=(size, 3)) @ noise_factor.T / np.sqrt(500000.)
        beta = ld @ joint + noise
        sx = np.sqrt(2 * af * (1 - af))
        beta_parts.append(beta / sx[:, None])
        se_parts.append(np.full((size, 3), 1 / np.sqrt(500000.)) / sx[:, None])
        snp_parts.extend(snp)
        af_parts.extend(af)
        path = out / "ld" / f"block_{b:03d}.npz"
        np.savez_compressed(path, ld=ld, snp=snp, a1=np.repeat("A", size),
                            a2=np.repeat("C", size), af=af, chrom=np.repeat("1", size),
                            bp=np.arange(offset + 1, offset + size + 1) * 1000,
                            ld_origin=np.array("generating_truth"))
        ld_rows.append([f"b{b:03d}", str(path), "npz"])
        offset += size
    beta, se = np.vstack(beta_parts), np.vstack(se_parts)
    with open(out / "ld_manifest.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["block", "path", "format"])
        writer.writerows(ld_rows)
    with open(out / "sumstats.tsv", "w", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t")
        writer.writerow(["trait", "path", "effect_scale"])
        for a, trait in enumerate(["A", "B", "C"]):
            path = out / f"{trait}.tsv"
            writer.writerow([trait, str(path), "allele"])
            with open(path, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["SNP", "A1", "A2", "BETA", "SE", "EAF"])
                for j in range(n):
                    w.writerow([snp_parts[j], "A", "C", beta[j, a], se[j, a], af_parts[j]])
    with open(out / "noise.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", "A", "B", "C"])
        for i, trait in enumerate(["A", "B", "C"]):
            w.writerow([trait, *overlap[i]])
    (out / "truth.json").write_text(json.dumps({
        "seed": args.seed, "n_snps": n, "blocks": len(sizes),
        "paths": {"A_to_B": .4, "B_to_C": .25},
        "block_size_range": [int(min(sizes)), int(max(sizes))],
        "noise_intercept": overlap.tolist(),
        "note": "Toy AR LD; the separate validation uses real HapMap3 LD with GWASBrewer."
    }, indent=2) + "\n")
    print(out)


if __name__ == "__main__":
    main()
