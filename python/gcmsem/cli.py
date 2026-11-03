"""Command line front end with auditable LDSC-style stage logging."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import platform
import shlex
import sys
import time

import numpy as np
import scipy

from . import __version__
from .io import harmonize_sumstats, sha256_file, write_results, _rows
from .moments import estimate_moments
from .scores import add_score_parser, load_scores


def _read_noise(path, traits):
    rows = list(_rows(path))
    expected = [t.lower() for t in traits]
    if len(set(expected)) != len(expected):
        raise ValueError("noise-intercept table requires trait names distinct ignoring case")
    found = {}
    for number, row in rows:
        if "trait" not in row or set(row) != {"trait", *expected}:
            raise ValueError("noise-intercept file must be a labeled square TSV: trait followed by all trait names")
        name = row["trait"]
        if name in found:
            raise ValueError(f"duplicate noise-intercept row {name}")
        found[name] = [float(row[t]) for t in expected]
    if set(found) != set(traits):
        raise ValueError("noise-intercept row names do not match trait manifest")
    return np.array([found[t] for t in traits])


def _read_tuples(path, traits):
    indices = {t: i for i, t in enumerate(traits)}
    answer = []
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            line = line.split("#", 1)[0].strip()
            if not line or line == "key":
                continue
            names = line.split("|")
            if any(t not in indices for t in names):
                raise ValueError(f"moment selection line {number} contains an unknown trait")
            answer.append(tuple(indices[t] for t in names))
    if not answer:
        raise ValueError("moment selection file contains no tuples")
    return answer


def estimate_cli(args):
    log = logging.getLogger("gcmsem")
    if not args.overwrite and Path(args.out + ".metadata.json").exists():
        raise FileExistsError("completed output already exists; choose a new --out or use --overwrite")
    scores = load_scores(args.scores)
    log.info("SCORES tags=%d reference_variants=%g annotations=%d", len(scores["snp"]),
             np.sum(scores["counts"]), len(scores["counts"]))
    beta, se, aligned, traits, qc = harmonize_sumstats(
        args.sumstats, scores, min_info=args.min_info, min_n=args.min_n,
        af_tolerance=args.af_tolerance, allow_strand_flips=args.allow_strand_flips,
        keep_palindromic=args.keep_palindromic, logger=log)
    result = estimate_moments(beta, se, aligned, traits=traits, max_order=args.max_order,
        moment_tuples=None if args.moments is None else _read_tuples(args.moments, traits),
        noise_intercept=None if args.noise_intercept is None else _read_noise(args.noise_intercept, traits),
        n_blocks=args.blocks, max_moments=args.max_moments, max_memory_mb=args.max_memory_mb,
        max_condition=args.max_condition, logger=log)
    provenance = {
        "gcmsem_version": __version__, "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "command": getattr(args, "command", None),
        "sumstats_manifest": str(Path(args.sumstats).resolve()),
        "sumstats_manifest_sha256": sha256_file(args.sumstats),
        "score_file": str(Path(args.scores).resolve()), "score_sha256": sha256_file(args.scores),
        "score_metadata": scores.get("metadata", {}),
        "noise_intercept_file": None if args.noise_intercept is None else {
            "path": str(Path(args.noise_intercept).resolve()), "sha256": sha256_file(args.noise_intercept)},
        "moment_selection_file": None if args.moments is None else {
            "path": str(Path(args.moments).resolve()), "sha256": sha256_file(args.moments)},
        "qc": qc,
    }
    result["metadata"].update(provenance)
    paths = write_results(result, args.out, overwrite=args.overwrite,
                          conditional=not args.no_conditional)
    for path in paths:
        log.info("OUTPUT %s", path)
    return 0


def parser():
    ap = argparse.ArgumentParser(prog="gcmsem", description="Genomic comoment estimation from GWAS summary statistics")
    ap.add_argument("--version", action="version", version=f"gcmsem {__version__}")
    sub = ap.add_subparsers(dest="subcommand", required=True)
    add_score_parser(sub)
    e = sub.add_parser("estimate", help="harmonize N GWAS and estimate moments through order four")
    e.add_argument("--sumstats", required=True, help="TSV trait/path/effect_scale manifest")
    e.add_argument("--scores", required=True, help="validated signed LD-power score NPZ")
    e.add_argument("--out", required=True, help="output filename prefix")
    e.add_argument("--max-order", type=int, choices=(2, 3, 4), default=4)
    e.add_argument("--moments", help="optional file of trait|trait|... tuples; lower orders retained")
    e.add_argument("--blocks", type=int, default=200, help="maximum genomic jackknife groups (whole LD blocks)")
    e.add_argument("--noise-intercept", help="labeled known noise calibration matrix; default estimated from LD regressions")
    e.add_argument("--min-info", type=float)
    e.add_argument("--min-n", type=float)
    e.add_argument("--af-tolerance", type=float, default=0.1)
    e.add_argument("--allow-strand-flips", action="store_true")
    e.add_argument("--keep-palindromic", action="store_true", help="require uniquely resolving EAF; ambiguous palindromes still dropped")
    e.add_argument("--max-moments", type=int, default=5000)
    e.add_argument("--max-memory-mb", type=float, default=2048)
    e.add_argument("--max-condition", type=float, default=1e6)
    e.add_argument("--no-conditional", action="store_true", help="omit separate annotation-specific output bundles")
    e.add_argument("--overwrite", action="store_true")
    e.set_defaults(func=estimate_cli)
    return ap


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    args = parser().parse_args(argv)
    args.command = shlex.join(["gcmsem", *argv])
    log = logging.getLogger("gcmsem")
    log.setLevel(logging.INFO)
    log.propagate = False
    for handler in list(log.handlers):
        handler.close()
        log.removeHandler(handler)
    formatter = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", datefmt="%Y-%m-%dT%H:%M:%SZ")
    formatter.converter = time.gmtime
    stream = logging.StreamHandler(sys.stderr)
    stream.setFormatter(formatter)
    log.addHandler(stream)
    started = time.monotonic()
    try:
        if getattr(args, "out", None):
            logfile = Path(str(args.out) + ".log")
            logfile.parent.mkdir(parents=True, exist_ok=True)
            handle = logging.FileHandler(logfile, mode="a")
            handle.setFormatter(formatter)
            log.addHandler(handle)
        log.info("GCMSEM %s | Python %s | NumPy %s | SciPy %s", __version__, platform.python_version(), np.__version__, scipy.__version__)
        log.info("COMMAND %s", args.command)
        code = args.func(args)
        log.info("SUCCESS elapsed_seconds=%.3f", time.monotonic() - started)
        return 0 if code is None else code
    except (ValueError, OSError, KeyError, FloatingPointError, np.linalg.LinAlgError) as exc:
        log.error("FAILED %s: %s", type(exc).__name__, exc)
        return 2
    finally:
        for handler in list(log.handlers):
            handler.flush()


if __name__ == "__main__":
    raise SystemExit(main())
