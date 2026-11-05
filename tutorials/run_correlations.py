#!/usr/bin/env python3
"""Recompute GCMSEM correlations; optionally run a matched-reference sensitivity.

The primary standard-LDSC comparison is run_standard_ldsc.py. LDSC stages in
THIS script intentionally use custom matched scores and are sensitivity analyses
only; they must never be described as standard LDSC.

The released panel is already harmonized to the signed-score allele orientation.
It is NOT a shortcut for harmonizing new GWAS, and includes trusted palindromes.
No reference genotypes are needed by this script. Large outputs belong outside
this source repository. See the accompanying wiki for model/scale assumptions.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import logging
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

PIN = "1f09cf0c5e61a38a281016fefcc1e0583cc78402"
NAMES = ["DrinksPerWeek", "VitaminD", "SBP", "Insomnia", "Depression", "LDL",
         "FEV1", "HeartAttack", "Angina", "Height", "Weight", "BMI", "DBP", "FVC", "Waist"]
CODES = ["GSCAN_DrinksPerWeek", "30890_raw", "4080_irnt", "1200", "20544_11",
         "30780_raw", "3063_irnt", "6150_1", "6150_2", "50_irnt", "21002_irnt",
         "21001_irnt", "4079_irnt", "3062_irnt", "48_irnt"]
MODES = ("unstratified", "maf5", "unstratified_twostep30")
SHARED = {}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def json_write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def finite(x):
    y = float(x)
    return y if np.isfinite(y) else None


def matrix_write(path, a):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", *NAMES])
        for name, row in zip(NAMES, a):
            w.writerow([name, *row])


def load_data(args):
    from gcmsem.scores import load_scores
    z = np.load(args.panel, allow_pickle=False)
    scores = load_scores(args.scores)
    if list(map(str, z["trait_codes"])) != CODES:
        raise ValueError("Panel trait order does not match the documented 15-trait release")
    rows = np.asarray(z["row_index"], np.int64)
    if rows.ndim != 1 or len(rows) == 0 or np.any(rows < 0) or np.any(np.diff(rows) <= 0) or rows[-1] >= len(scores["snp"]):
        raise ValueError("Invalid score-row mapping")
    for key in ("snp", "a1", "a2", "chrom", "bp"):
        if not np.array_equal(z[key].astype(str), np.asarray(scores[key])[rows].astype(str)):
            raise ValueError(f"Preharmonized panel/score mismatch in {key}")
    y, se, n = (np.asarray(z[k], float) for k in ("beta_std", "se_std", "n"))
    if y.shape != (len(rows), 15) or se.shape != y.shape or n.shape != y.shape:
        raise ValueError("Expected SNP by 15 arrays for beta, SE and N")
    if not all(np.isfinite(x).all() for x in (y, se, n)) or np.any(se <= 0) or np.any(n <= 0):
        raise ValueError("Invalid common panel values")
    aligned = {k: (v[rows] if isinstance(v, np.ndarray) and v.ndim and
                  len(v) == len(scores["snp"]) and k != "counts" else v)
               for k, v in scores.items()}
    # Source counts are the full source universe, never the filtered tag count.
    aligned["counts"] = scores["counts"]
    return y, se, n, aligned, rows


def run_gcmsem(args, y, se, scores, log):
    from gcmsem.moments import estimate_moments, jackknife_cov
    from gcmsem.io import write_results
    path = args.out / "gcmsem_arrays.npz"
    if path.exists():
        log.info("RESUME existing GCMSEM output (input/configuration hashes checked)")
        return
    fit = estimate_moments(y, se, scores, traits=NAMES, max_order=2,
                           n_blocks=args.blocks, max_memory_mb=args.max_memory_mb, logger=log)
    fit["metadata"].update({"tutorial_input": "pre-harmonized real GWAS panel",
        "score_metadata": scores["metadata"],
        "palindrome_policy": "retained from archived source-allele harmonization; no invented GWAS EAF",
        "deletion_warning": "Score reference uses rolling 1 Mb LD, so deletion group boundaries can have crossing LD",
        "panel_sha256": sha256(args.panel), "scores_sha256": sha256(args.scores)})
    write_results(fit, str(args.out / "gcmsem"), conditional=True)
    ids = {t: k for k, t in enumerate(fit["indices"])}
    cov = np.array([[fit["estimates"][ids[tuple(sorted((i, j)))]] for j in range(15)] for i in range(15)])
    jk_cov = np.array([[[v[ids[tuple(sorted((i, j)))]] for j in range(15)] for i in range(15)] for v in fit["jackknife"]])
    diag = np.diag(cov)
    if np.any(diag <= 0) or np.any(np.diagonal(jk_cov, axis1=1, axis2=2) <= 0):
        raise ValueError("Nonpositive SNP-effect variance: retain fit outputs, but cannot form all correlations")
    rg = cov / np.sqrt(np.outer(diag, diag))
    jdiag = np.diagonal(jk_cov, axis1=1, axis2=2)
    jk_rg = jk_cov / np.sqrt(jdiag[:, :, None] * jdiag[:, None, :])
    pairs = list(itertools.combinations(range(15), 2))
    a, b = np.asarray(pairs).T
    v = jackknife_cov(jk_rg[:, a, b], rg[a, b], np.asarray(fit["metadata"]["jackknife_fractions"]))
    rse = np.zeros((15, 15))
    rse[a, b] = rse[b, a] = np.sqrt(np.diag(v))
    matrix_write(args.out / "gcmsem_rg.tsv", rg)
    matrix_write(args.out / "gcmsem_rg_se.tsv", rse)
    np.savez_compressed(path, covariance=cov, rg=rg, rg_se=rse, rg_sampling_cov=v,
                        rg_jackknife=jk_rg, traits=np.array(NAMES))
    log.info("GCMSEM_FINISHED pairs=%d min_rg_eigenvalue=%g", len(pairs), np.linalg.eigvalsh(rg)[0])


def ldsc_fit(task):
    from ldscore.regressions import RG
    mode, a, b = task
    d, start = SHARED, time.monotonic()
    row = {"specification": mode, "i": a, "j": b, "trait1": NAMES[a], "trait2": NAMES[b],
           "status": "succeeded", "n_snps": len(d["z"])}
    try:
        q = d["q"] if mode == "maf5" else d["qsum"]
        m = d["counts"] if mode == "maf5" else d["countsum"]
        f = RG(d["z"][:, a, None], d["z"][:, b, None], q, d["weights"],
               d["n"][:, a, None], d["n"][:, b, None], m, n_blocks=d["blocks"],
               twostep=30 if mode == "unstratified_twostep30" else None)
        row.update(rg=finite(f.rg_ratio), se=finite(f.rg_se),
                   intercept=finite(f.gencov.intercept), intercept_se=finite(f.gencov.intercept_se),
                   h2_1=finite(f.hsq1.tot), h2_2=finite(f.hsq2.tot),
                   negative_h2=bool(f._negative_hsq),
                   outside_unit_interval=bool(np.isfinite(f.rg_ratio) and abs(f.rg_ratio) > 1))
        deletes = np.column_stack([f.hsq1.tot_delete_values.ravel(), f.hsq2.tot_delete_values.ravel(),
                                   f.gencov.tot_delete_values.ravel()])
        with np.errstate(invalid="ignore", divide="ignore"):
            rg = deletes[:, 2] / np.sqrt(deletes[:, 0] * deletes[:, 1])
        row["nonpositive_variance_deletions"] = int(np.sum(np.any(deletes[:, :2] <= 0, axis=1)))
        np.savez_compressed(d["out"] / "ldsc_deletions" / f"{mode}_{a:02d}_{b:02d}.npz",
                            estimates=deletes, rg=rg)
        if row["rg"] is None or row["se"] is None:
            row["status"] = "invalid_estimate"
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    row["seconds"] = time.monotonic() - start
    return row


def run_ldsc(args, y, se, n, scores, rows, log):
    if args.ldsc_dir is None or args.weights is None:
        raise ValueError("LDSC requires --ldsc-dir and --weights (HM3-only regression weighting LD)")
    commit = subprocess.check_output(["git", "-C", str(args.ldsc_dir), "rev-parse", "HEAD"], text=True).strip()
    if commit != PIN:
        raise ValueError(f"LDSC must be checked out at {PIN}; got {commit}")
    # Verify the statistical implementation is unmodified; a CLI pandas syntax patch is allowed.
    subprocess.run(["git", "-C", str(args.ldsc_dir), "diff", "--exit-code", "HEAD", "--",
                    "ldscore/regressions.py", "ldscore/irwls.py", "ldscore/jackknife.py"], check=True)
    ldsc_config = {"weights_sha256": sha256(args.weights), "commit": commit,
                   "modes": list(MODES), "blocks": args.blocks,
                   "statistical_source_sha256": {f: sha256(args.ldsc_dir / "ldscore" / f)
                       for f in ("regressions.py", "irwls.py", "jackknife.py")}}
    ldcp = args.out / "ldsc_configuration.json"
    if ldcp.exists() and json.loads(ldcp.read_text()) != ldsc_config:
        raise ValueError("Existing LDSC output uses different weights or implementation; use a fresh --out")
    if (args.out / "ldsc_fits.jsonl").exists() and not ldcp.exists():
        raise ValueError("Cannot resume LDSC journal without its configuration; use a fresh --out")
    json_write(ldcp, ldsc_config)
    sys.path.insert(0, str(args.ldsc_dir))
    from ldscore.regressions import RG  # fail before submitting workers if dependencies missing
    wz = np.load(args.weights, allow_pickle=False)
    if "snp" not in wz.files or "L2" not in wz.files:
        raise ValueError("Weight archive requires snp and L2 arrays")
    if np.array_equal(wz["snp"].astype(str), scores["snp"].astype(str)):
        weights = np.asarray(wz["L2"], float)
    elif np.array_equal(wz["snp"][rows].astype(str), scores["snp"].astype(str)):
        weights = np.asarray(wz["L2"][rows], float)
    else:
        raise ValueError("LDSC weighting score SNP order mismatch")
    if weights.ndim == 2:
        weights = weights.sum(axis=1)
    if weights.shape != (len(y),) or not np.isfinite(weights).all():
        raise ValueError("Invalid weighting LD values")
    weights = weights[:, None]
    q = np.asarray(scores["L2"], float)
    (args.out / "ldsc_deletions").mkdir(exist_ok=True)
    SHARED.update(z=y / se, n=np.asarray(n, np.float64), q=q, qsum=q.sum(axis=1)[:, None],
                  counts=np.asarray(scores["counts"], float)[None],
                  countsum=np.array([[np.sum(scores["counts"])]], float), weights=weights,
                  blocks=args.blocks, out=args.out)
    journal = args.out / "ldsc_fits.jsonl"
    records = [json.loads(x) for x in journal.read_text().splitlines()] if journal.exists() else []
    done = {(r["specification"], r["i"], r["j"]) for r in records if r["status"] == "succeeded"}
    tasks = [(mode, a, b) for mode in MODES for a, b in itertools.combinations(range(15), 2)
             if (mode, a, b) not in done]
    with journal.open("a") as f, ProcessPoolExecutor(max_workers=args.workers,
                mp_context=multiprocessing.get_context("fork")) as pool:
        for future in as_completed([pool.submit(ldsc_fit, t) for t in tasks]):
            r = future.result()
            f.write(json.dumps(r, allow_nan=False) + "\n")
            f.flush()
            log.info("LDSC %s %s/%s %s rg=%s se=%s %.2fs", r["specification"],
                     r["trait1"], r["trait2"], r["status"], r.get("rg"), r.get("se"), r["seconds"])
    json_write(args.out / "ldsc_provenance.json", {
        "repository": "https://github.com/CBIIT/ldsc", "commit": commit,
        "package_version": "3.0.2", "cli_banner_version": "3.0.1",
        "statistics": "original RG class; regression/IRWLS/jackknife source unmodified",
        "weights_sha256": sha256(args.weights), "target_counts": scores["counts"].tolist(),
        "N_dtype": "float64 (prevents upstream N1*N2 int32 overflow)",
        "modes": list(MODES), "blocks": args.blocks,
        "chi_square_filter": "none; all common complete-case rows retained",
        "liability_conversion": False,
        "one_stage_note": "twostep=None matches the executed pinned CLI default; explicit twostep=30 sensitivity is separate"})


def summarize(args):
    fit = np.load(args.out / "gcmsem_arrays.npz", allow_pickle=False)
    records = {}
    for line in (args.out / "ldsc_fits.jsonl").read_text().splitlines():
        r = json.loads(line)
        records[(r["specification"], r["i"], r["j"])] = r
    pairs, summary = [], {}
    for mode in MODES:
        mat, semat = np.eye(15), np.zeros((15, 15))
        x, y, failed = [], [], []
        for a, b in itertools.combinations(range(15), 2):
            r = records.get((mode, a, b))
            if r is None or r["status"] != "succeeded":
                mat[a, b] = mat[b, a] = np.nan
                semat[a, b] = semat[b, a] = np.nan
                failed.append([NAMES[a], NAMES[b]])
                continue
            mat[a, b] = mat[b, a] = r["rg"]
            semat[a, b] = semat[b, a] = r["se"]
            gx = float(fit["rg"][a, b]); x.append(gx); y.append(r["rg"])
            pairs.append({"specification": mode, "trait1": NAMES[a], "trait2": NAMES[b],
                          "gcmsem_rg": gx, "gcmsem_se": float(fit["rg_se"][a, b]),
                          "ldsc_rg": r["rg"], "ldsc_se": r["se"], "difference": gx-r["rg"]})
        x, y = np.array(x), np.array(y)
        summary[mode] = {"pairs_expected": 105, "pairs_succeeded": len(x), "failed_pairs": failed,
            "pearson": finite(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else None,
            "rmse": finite(np.sqrt(np.mean((x-y)**2))) if len(x) else None,
            "mean_difference": finite(np.mean(x-y)) if len(x) else None,
            "max_absolute_difference": finite(np.max(abs(x-y))) if len(x) else None}
        matrix_write(args.out / f"ldsc_{mode}_rg.tsv", mat)
        matrix_write(args.out / f"ldsc_{mode}_rg_se.tsv", semat)
    with (args.out / "comparison.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["specification", "trait1", "trait2", "gcmsem_rg", "gcmsem_se", "ldsc_rg", "ldsc_se", "difference"], delimiter="\t")
        w.writeheader(); w.writerows(pairs)
    summary["interpretation"] = [
        "Matched-reference sensitivity only; this is not the primary standard-LDSC comparison.",
        "Descriptive agreement, not an independent validation of accuracy; pairs share traits.",
        "GCMSEM estimates centered SNP-effect correlation; LDSC estimates its genetic correlation model with iterative weights.",
        "MAF5 specifies the same five disjoint source bins, but does not make regression weights or mean treatment identical.",
        "GCMSEM uses chromosome-respecting weighted unequal-deletion SEs; LDSC uses native equal-count blocks that can cross chromosomes.",
        "No independent-method SE for the difference is reported: methods share data and deletion boundaries differ.",
        "Reference and score uncertainty, higher-order interval coverage and causal identification are not tested here."]
    json_write(args.out / "comparison_summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--weights", type=Path)
    p.add_argument("--ldsc-dir", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--stage", choices=["all", "gcmsem", "ldsc", "summarize"], default="gcmsem")
    p.add_argument("--blocks", type=int, default=200)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--max-memory-mb", type=float, default=4096)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    config = {"panel_sha256": sha256(args.panel), "scores_sha256": sha256(args.scores),
              "blocks": args.blocks, "traits": NAMES, "trait_codes": CODES}
    cp = args.out / "configuration.json"
    if cp.exists() and json.loads(cp.read_text()) != config:
        raise ValueError("Output folder belongs to different inputs/settings; use a fresh --out")
    json_write(cp, config)
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(level=logging.INFO, format="%(asctime)sZ %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.out / "run.log")])
    log = logging.getLogger("gcmsem")
    log.info("COMMAND %s", " ".join(sys.argv))
    if args.stage != "summarize":
        y, se, n, scores, rows = load_data(args)
        log.info("INPUT common_SNPs=%d traits=%d source_variants=%d", len(y), y.shape[1], np.sum(scores["counts"]))
        if args.stage in ("all", "gcmsem"):
            run_gcmsem(args, y, se, scores, log)
        if args.stage in ("all", "ldsc"):
            run_ldsc(args, y, se, n, scores, rows, log)
    if args.stage in ("all", "summarize"):
        summarize(args)


if __name__ == "__main__":
    main()
