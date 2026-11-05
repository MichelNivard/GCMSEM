#!/usr/bin/env python3
"""Run standard upstream LDSC munging and genetic-correlation commands.

Primary comparator: original bulik/ldsc, unmodified Python 2 source, canonical
published eur_w_ld_chr reference AND weighting scores, official w_hm3 list,
and stock statistical defaults. No custom LD scores, source counts, intercept
settings, chi-square thresholds, allele exemptions, or common-SNP mask.
The documented --chunksize 500000 controls I/O only and avoids the old pandas
quadratic string-membership path above one million rows per chunk.

Inputs retain GCMSEM's existing HM3 selection, as requested, separately for each
trait. prepare_standard_ldsc_inputs.py formats those public summary statistics for
stock munging and derives two-sided normal p-values from the archived beta/SE.
No reference-panel allele frequency is represented as a GWAS study frequency.
LDSC then applies its ordinary per-trait and pairwise QC; there is no forced
15-trait complete-case mask in the LDSC run.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import gzip
import tempfile
import shutil
import itertools
import json
import logging
import re
from pathlib import Path
import shlex
import subprocess
import time
import numpy as np
from run_correlations import NAMES, CODES, sha256, json_write, matrix_write

PIN = "2fdeeb3b44379408794154993dbd6101b8946b7e"


def execute(command, stdout_path, log):
    log.info("COMMAND %s", shlex.join(map(str, command)))
    with stdout_path.open("w") as f:
        subprocess.run(list(map(str, command)), stdout=f, stderr=subprocess.STDOUT, check=True)


def run_munge(args, log):
    def one(code):
        source = args.inputs / (code + ".tsv.gz")
        target = args.out / "munged" / code
        work = Path(tempfile.mkdtemp(prefix=code + "-", dir=str(args.out / "work")))
        produced = work / code
        command = [args.python2, args.ldsc_dir / "munge_sumstats.py", "--sumstats", source,
                   "--merge-alleles", args.hm3, "--snp", "SNP", "--a1", "A1", "--a2", "A2",
                   "--p", "P", "--N-col", "N", "--signed-sumstats", "BETA,0", "--chunksize", "500000", "--out", produced]
        if target.with_suffix(".sumstats.gz").exists() and target.with_suffix(".success.json").exists():
            old = json.loads(target.with_suffix(".success.json").read_text())
            if old["input_sha256"] != sha256(source) or old["output_sha256"] != sha256(target.with_suffix(".sumstats.gz")):
                raise ValueError("Munged input changed since its recorded successful run")
            return code, "resumed"
        execute(command, produced.with_suffix(".stdout.log"), log)
        product = produced.with_suffix(".sumstats.gz")
        if not product.is_file():
            raise RuntimeError(f"Munging did not produce {product}")
        with gzip.open(product, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        output_sha = sha256(product)
        for suffix in (".sumstats.gz", ".log", ".stdout.log"):
            produced.with_suffix(suffix).replace(target.with_suffix(suffix))
        if sha256(target.with_suffix(".sumstats.gz")) != output_sha:
            raise RuntimeError("Published munged file checksum differs")
        json_write(target.with_suffix(".success.json"), {"command": list(map(str, command)),
                   "input_sha256": sha256(source), "output_sha256": output_sha})
        return code, "completed"
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(one, c) for c in CODES]):
            log.info("MUNGE %s %s", *future.result())


def parse_summary(path, first):
    text = path.read_text()
    marker = "Summary of Genetic Correlation Results"
    if marker not in text:
        raise RuntimeError(f"No genetic correlation summary in {path}")
    lines = text.rsplit(marker, 1)[1].strip().splitlines()
    header = lines[0].split()
    if header[:4] != ["p1", "p2", "rg", "se"]:
        raise RuntimeError(f"Unexpected LDSC summary header in {path}: {header}")
    names = {code + ".sumstats.gz": i for i, code in enumerate(CODES)}
    snp_counts = {}
    for match in re.finditer(r"Computing rg for phenotype (\d+)/(\d+)\n([\s\S]*?)(?=Computing rg for phenotype|Summary of Genetic Correlation Results|\Z)", text):
        j = first + int(match.group(1)) - 1
        values = re.findall(r"(\d+) SNPs with valid alleles", match.group(3))
        if values:
            snp_counts[j] = int(values[-1])
    records = []
    for line in lines[1:]:
        cols = line.split()
        if len(cols) != len(header) or Path(cols[0]).name not in names:
            continue
        record = dict(zip(header, cols))
        i, j = names[Path(record["p1"]).name], names[Path(record["p2"]).name]
        if i != first or j <= i:
            raise RuntimeError("Unexpected phenotype ordering in LDSC log")
        out = {"i": i, "j": j, "trait1": NAMES[i], "trait2": NAMES[j], "code1": CODES[i], "code2": CODES[j]}
        for key in header[2:]:
            try:
                value = float(record[key])
            except ValueError:
                value = np.nan
            out[key] = value if np.isfinite(value) else None
        out["status"] = "succeeded" if out["rg"] is not None and out["se"] is not None else "invalid_estimate"
        out["outside_unit_interval"] = bool(out["rg"] is not None and abs(out["rg"]) > 1)
        out["source_log"] = str(path.name)
        out["n_snps"] = snp_counts.get(j)
        if out["status"] == "succeeded" and out["n_snps"] is None:
            raise RuntimeError("Successful standard LDSC fit lacks its allele-QC SNP count")
        records.append(out)
    expected = {(first, j) for j in range(first + 1, 15)}
    if {(r["i"], r["j"]) for r in records} != expected:
        raise RuntimeError(f"Incomplete LDSC summary in {path}")
    return records


def run_fit(args, log):
    def one(i):
        target = args.out / "fits" / CODES[i]
        command = [args.python2, args.ldsc_dir / "ldsc.py", "--rg",
                   ",".join(str(args.out / "munged" / (c + ".sumstats.gz")) for c in CODES[i:]),
                   "--ref-ld-chr", args.ref_prefix, "--w-ld-chr", args.ref_prefix, "--out", target]
        done = target.with_suffix(".success.json")
        if not done.exists():
            execute(command, target.with_suffix(".stdout.log"), log)
            records = parse_summary(target.with_suffix(".log"), i)
            json_write(done, {"command": list(map(str, command)), "pairs": records,
                              "log_sha256": sha256(target.with_suffix(".log"))})
        else:
            old = json.loads(done.read_text())
            if old["log_sha256"] != sha256(target.with_suffix(".log")):
                raise ValueError("LDSC log changed since its recorded successful run")
            records = parse_summary(target.with_suffix(".log"), i)
        log.info("LDSC %s pairs=%d succeeded=%d", CODES[i], len(records), sum(r["status"] == "succeeded" for r in records))
        return records
    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(one, i) for i in range(14)]):
            records.extend(future.result())
    records.sort(key=lambda r: (r["i"], r["j"]))
    if len(records) != 105:
        raise RuntimeError("Expected exactly 105 standard LDSC pairs")
    json_write(args.out / "all_pairs.json", records)
    with (args.out / "all_pairs.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t")
        w.writeheader(); w.writerows(records)
    rg, se = np.eye(15), np.zeros((15, 15))
    for r in records:
        i, j = r["i"], r["j"]
        rg[i, j] = rg[j, i] = np.nan if r["rg"] is None else r["rg"]
        se[i, j] = se[j, i] = np.nan if r["se"] is None else r["se"]
    matrix_write(args.out / "ldsc_standard_rg.tsv", rg)
    matrix_write(args.out / "ldsc_standard_rg_se.tsv", se)
    json_write(args.out / "status.json", {"pairs_expected": 105, "pairs_run": len(records),
               "pairs_succeeded": sum(r["status"] == "succeeded" for r in records),
               "n_snps_min": min(r["n_snps"] for r in records if r["n_snps"] is not None),
               "n_snps_max": max(r["n_snps"] for r in records if r["n_snps"] is not None),
               "outside_unit_interval": sum(r["outside_unit_interval"] for r in records),
               "values": "Original CLI summary precision (four decimal places); unrounded underlying outputs not substituted",
               "policy": "Invalid or out-of-range fits retained and flagged; no clipping/projection"})
    if any(r["status"] != "succeeded" for r in records):
        raise RuntimeError("Some standard LDSC pairs failed; outputs and complete logs retained")


def compare(args):
    if args.gcmsem_results is None:
        raise ValueError("Comparison requires --gcmsem-results")
    g = np.load(args.gcmsem_results / "gcmsem_arrays.npz", allow_pickle=False)
    if list(map(str, g["traits"])) != NAMES:
        raise ValueError("GCMSEM correlation trait order differs from documented panel")
    rows = json.loads((args.out / "all_pairs.json").read_text())
    pairs = []
    for r in rows:
        i, j = r["i"], r["j"]
        if r["status"] == "succeeded":
            pairs.append({"trait1": NAMES[i], "trait2": NAMES[j], "gcmsem_rg": float(g["rg"][i,j]),
                          "gcmsem_se": float(g["rg_se"][i,j]), "ldsc_standard_rg": r["rg"],
                          "ldsc_standard_se": r["se"], "difference": float(g["rg"][i,j]-r["rg"])})
    with (args.out / "comparison.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["trait1", "trait2", "gcmsem_rg", "gcmsem_se", "ldsc_standard_rg", "ldsc_standard_se", "difference"], delimiter="\t")
        w.writeheader(); w.writerows(pairs)
    x = np.array([r["gcmsem_rg"] for r in pairs]); y = np.array([r["ldsc_standard_rg"] for r in pairs])
    summary = {"comparison": "GCMSEM vs standard published-reference upstream LDSC CLI defaults",
        "gcmsem_arrays_sha256": sha256(args.gcmsem_results / "gcmsem_arrays.npz"),
        "gcmsem_metadata_sha256": sha256(args.gcmsem_results / "gcmsem.metadata.json"),
        "standard_configuration_sha256": sha256(args.out / "configuration.json"),
        "standard_all_pairs_sha256": sha256(args.out / "all_pairs.json"),
        "n_pairs": len(pairs), "pearson": float(np.corrcoef(x,y)[0,1]) if len(x)>1 else None,
        "rmse": float(np.sqrt(np.mean((x-y)**2))) if len(x) else None,
        "mean_difference": float(np.mean(x-y)) if len(x) else None,
        "max_absolute_difference": float(np.max(abs(x-y))) if len(x) else None,
        "interpretation": [
            "Primary LDSC comparator uses original unmodified bulik/ldsc with its statistical defaults and official published reference/weight scores.",
            "Each trait retains the requested upstream GCMSEM HM3 selection, then is independently munged with stock LDSC; LDSC selects its ordinary pairwise SNP intersections.",
            "GCMSEM uses its fixed common HM3 panel and five-bin signed-score working model. LDSC keeps the requested upstream HM3 selection but uses its own ordinary QC and published score model.",
            "The earlier common-reference/common-row comparison is a separately labelled sensitivity, not standard LDSC.",
            "Pairs share traits and the methods share GWAS: concordance is descriptive, not independent accuracy or causal validation.",
            "No difference SE is reported because outputs are statistically dependent and use different SNP sets/deletion groups.",
            "LDSC results are quoted at stock CLI reporting precision; no clipping or PSD projection."]}
    json_write(args.out / "comparison_summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--hm3", type=Path, required=True)
    ap.add_argument("--ref-prefix", type=Path, required=True, help="Official eur_w_ld_chr/ chromosome-prefix directory")
    ap.add_argument("--ldsc-dir", type=Path, required=True)
    ap.add_argument("--python2", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--gcmsem-results", type=Path)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--stage", choices=["all", "munge", "fit", "compare"], default="all")
    args = ap.parse_args()
    for key in ("inputs", "hm3", "ref_prefix", "ldsc_dir", "python2", "out"):
        setattr(args, key, getattr(args, key).resolve())
    # A Path removes the trailing separator, while LDSC appends chromosome numbers.
    args.ref_prefix = str(args.ref_prefix) + "/"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "munged").mkdir(exist_ok=True); (args.out / "fits").mkdir(exist_ok=True)
    (args.out / "work").mkdir(exist_ok=True)
    commit = subprocess.check_output(["git", "-C", str(args.ldsc_dir), "rev-parse", "HEAD"], text=True).strip()
    if commit != PIN:
        raise ValueError(f"Expected official upstream commit {PIN}, got {commit}")
    subprocess.run(["git", "-C", str(args.ldsc_dir), "diff", "--exit-code", "HEAD", "--", "*.py", "ldscore"], check=True)
    refs = sorted(Path(args.ref_prefix).glob("*"))
    refs = [p for p in refs if p.is_file()]
    config = {"upstream": "https://github.com/bulik/ldsc", "commit": commit,
        "statistical_settings": "stock defaults; identical official eur_w_ld_chr prefix for --ref-ld-chr and --w-ld-chr",
        "python_version": subprocess.check_output([str(args.python2), "--version"], stderr=subprocess.STDOUT, text=True).strip(),
        "numeric_environment": json.loads(subprocess.check_output([str(args.python2), "-c",
            "import json, pkg_resources; print(json.dumps(dict((p,pkg_resources.get_distribution(p).version) for p in ['numpy','scipy','pandas','bitarray'])))"], text=True)),
        "inputs_sha256": {c: sha256(args.inputs / (c + ".tsv.gz")) for c in CODES},
        "reference_files_sha256": {p.name: sha256(p) for p in refs}, "hm3_sha256": sha256(args.hm3)}
    cp = args.out / "configuration.json"
    if cp.exists() and json.loads(cp.read_text()) != config:
        raise ValueError("Output belongs to different sources/reference/software; use a new --out")
    json_write(cp, config)
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(level=logging.INFO, format="%(asctime)sZ %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.out / "run.log")])
    log = logging.getLogger("standard_ldsc")
    if args.stage in ("all", "munge"): run_munge(args, log)
    if args.stage in ("all", "fit"): run_fit(args, log)
    if args.stage in ("all", "compare"): compare(args)


if __name__ == "__main__":
    main()
