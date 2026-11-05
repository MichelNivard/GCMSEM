#!/usr/bin/env python3
"""Format the requested GCMSEM-selected HM3 inputs for stock LDSC munging.

No cross-trait complete-case mask, significance filtering, allele recoding or
new reference selection is applied here. P is derived from the supplied Z=BETA/SE
using the two-sided normal test. These are not claimed to be original GWAS P
strings. Underflow P=0 is retained for ordinary LDSC P-value QC; no clipping.
REF_AF is a reference frequency, so it is deliberately NOT sent as study FRQ.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import io
import json
import os
import numpy as np
import pandas as pd
from scipy.special import ndtr
from run_correlations import CODES, NAMES, sha256, json_write


def prepare(task):
    code, name, source_dir, out = task
    source = source_dir / (code + ".tsv.gz")
    target = out / (code + ".tsv.gz")
    source_sha = sha256(source)
    success = out / (code + ".json")
    if success.exists() and target.exists():
        old = json.loads(success.read_text())
        if old["input_sha256"] != source_sha or old["output_sha256"] != sha256(target):
            raise ValueError("Existing converted input differs; choose a fresh output directory")
        return old
    d = pd.read_csv(source, sep="\t", usecols=["SNP", "A1", "A2", "BETA", "SE", "N"])
    z = d["BETA"].to_numpy(float) / d["SE"].to_numpy(float)
    if not np.isfinite(z).all() or np.any(d["SE"] <= 0) or np.any(d["N"] <= 0):
        raise ValueError("Invalid source BETA/SE/N")
    d["P"] = 2 * ndtr(-abs(z))
    d = d[["SNP", "A1", "A2", "BETA", "P", "N"]]
    temp = target.with_suffix(target.suffix + ".part")
    with temp.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=3, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                d.to_csv(text, sep="\t", index=False)
        raw.flush(); os.fsync(raw.fileno())
    with gzip.open(temp, "rb") as f:
        for _ in iter(lambda: f.read(1024 * 1024), b""):
            pass
    temp.replace(target)
    result = {"code": code, "name": name, "path": target.name,
              "input": source.name, "input_sha256": source_sha, "output_sha256": sha256(target),
              "rows": len(d), "zero_p_underflow": int(np.sum(d["P"] == 0)),
              "p_value": "2*Phi(-abs(BETA/SE)); normal approximation derived from retained beta and SE",
              "selection": "existing GCMSEM HM3 selection per trait; no extra filter or cross-trait intersection",
              "study_frequency_supplied": False, "imputation_INFO_supplied": False}
    json_write(success, result)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sumstats", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(prepare, (c,n,args.sumstats,args.out)) for c,n in zip(CODES,NAMES)]):
            r = future.result(); rows.append(r)
            print(r["code"], r["rows"], "zero_P", r["zero_p_underflow"], flush=True)
    rows.sort(key=lambda r: CODES.index(r["code"]))
    json_write(args.out / "manifest.json", {"traits": rows,
        "instructions": "Run stock munge_sumstats.py with official w_hm3.snplist and standard QC defaults",
        "source": "Public GWAS summary statistics; existing GCMSEM HM3-selected rows retained at user request",
        "caution": "The upstream HM3 selection is shared by request; study FRQ/INFO were not archived, so they are not invented or replaced by reference AF."})


if __name__ == "__main__":
    main()
