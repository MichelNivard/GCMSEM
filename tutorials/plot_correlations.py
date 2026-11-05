#!/usr/bin/env python3
"""Plot the primary comparison with stock published-reference LDSC."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True, help="run_standard_ldsc.py output directory")
    p.add_argument("--gcmsem-results", type=Path, required=True)
    args = p.parse_args()
    g = pd.read_csv(args.gcmsem_results / "gcmsem_rg.tsv", sep="\t", index_col=0)
    l = pd.read_csv(args.results / "ldsc_standard_rg.tsv", sep="\t", index_col=0)
    if not g.index.equals(l.index) or not g.columns.equals(l.columns):
        raise ValueError("Correlation matrix labels differ")
    if not np.isfinite(g.values).all() or not np.isfinite(l.values).all():
        raise ValueError("Incomplete matrix; inspect failed fits before plotting")
    labels = list(g.index)
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
    bound = max(1., np.max(abs(g.values)), np.max(abs(l.values)))
    diff = g.values - l.values
    dmax = max(.01, np.max(abs(diff)))
    for ax, a, title, v in zip(axes, [g.values, l.values, diff],
                             ["GCMSEM", "Standard LDSC", "GCMSEM − LDSC"], [bound, bound, dmax]):
        im = ax.imshow(a, cmap="RdBu_r", vmin=-v, vmax=v)
        ax.set_xticks(range(15), labels, rotation=65, ha="right", fontsize=8)
        ax.set_yticks(range(15), labels, fontsize=8)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=.6)
    fig.suptitle("Real HM3-selected GWAS: 15 traits, 105 pairwise correlations", fontsize=14)
    fig.savefig(args.results / "correlation_matrices.png", dpi=180)
    fig.savefig(args.results / "correlation_matrices.pdf")
    plt.close(fig)
    d = pd.read_csv(args.results / "comparison.tsv", sep="\t")
    s = json.loads((args.results / "comparison_summary.json").read_text())
    fig, ax = plt.subplots(figsize=(5.8, 5.5), constrained_layout=True)
    ax.scatter(d.ldsc_standard_rg, d.gcmsem_rg, s=25, alpha=.75, color="#28688f", edgecolors="none")
    lo = min(-.1, d.ldsc_standard_rg.min(), d.gcmsem_rg.min()) - .05
    hi = max(.1, d.ldsc_standard_rg.max(), d.gcmsem_rg.max()) + .05
    ax.plot([lo, hi], [lo, hi], color="#777", linestyle="--", linewidth=1)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel="Standard LDSC genetic correlation", ylabel="GCMSEM SNP-effect correlation")
    ax.set_title(f"Stock LDSC software, published scores and defaults\nPearson r = {s['pearson']:.4f}; RMSE = {s['rmse']:.4f}")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=.15)
    fig.supxlabel("Descriptive concordance; shared GWAS and dependent trait pairs", fontsize=9)
    fig.savefig(args.results / "correlation_concordance.png", dpi=180)
    fig.savefig(args.results / "correlation_concordance.pdf")


if __name__ == "__main__":
    main()
