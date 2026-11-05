#!/usr/bin/env python3
"""Check the pinned official LDSC command line against the panel RG API, height/weight.

Uses the same common SNP rows, reference LD, HM3 weighting LD and source count.
The CLI rounds rg and SE to four decimal places, so agreement is checked at that
reporting precision. Actual commands, inputs and stdout are retained externally.
"""
from pathlib import Path
import argparse
import json
import gzip
import io
import os
import re
import subprocess
import sys
import numpy as np
import pandas as pd
from run_correlations import load_data, sha256, PIN, json_write


def table_write(table, path):
    temporary = path.with_suffix(path.suffix + ".part")
    with open(temporary, "wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=3, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                table.to_csv(text, sep="\t", index=False)
        raw.flush()
        os.fsync(raw.fileno())
    with gzip.open(temporary, "rb") as handle:
        for _ in iter(lambda: handle.read(1024 * 1024), b""):
            pass
    temporary.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--weights", type=Path, required=True)
    ap.add_argument("--ldsc-dir", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    args = ap.parse_args()
    config = json.loads((args.results / "configuration.json").read_text())
    ldconfig = json.loads((args.results / "ldsc_configuration.json").read_text())
    if config["panel_sha256"] != sha256(args.panel) or config["scores_sha256"] != sha256(args.scores):
        raise ValueError("CLI inputs do not match the panel run being checked")
    if ldconfig["weights_sha256"] != sha256(args.weights):
        raise ValueError("CLI weights do not match the panel run being checked")
    y, se, n, scores, rows = load_data(args)
    commit = subprocess.check_output(["git", "-C", str(args.ldsc_dir), "rev-parse", "HEAD"], text=True).strip()
    if commit != PIN:
        raise ValueError("Unexpected LDSC source commit")
    w = np.load(args.weights, allow_pickle=False)
    if not np.array_equal(w["snp"][rows], scores["snp"]):
        raise ValueError("Weight SNP mismatch")
    out = (args.results / "cli_check").resolve()
    out.mkdir(exist_ok=True)
    base = pd.DataFrame({"CHR": scores["chrom"], "SNP": scores["snp"], "BP": scores["bp"],
                         "L2": scores["L2"].sum(axis=1)})
    table_write(base, out / "reference.l2.ldscore.gz")
    base["L2"] = w["L2"][rows]
    table_write(base, out / "weights.l2.ldscore.gz")
    (out / "reference.l2.M").write_text(str(int(np.sum(scores["counts"]))) + "\n")
    for index, name in [(9, "height"), (10, "weight")]:
        table_write(pd.DataFrame({"SNP": scores["snp"], "A1": scores["a1"], "A2": scores["a2"],
                      "Z": y[:, index] / se[:, index], "N": n[:, index]}),
                    out / (name + ".sumstats.gz"))
    cmd = [sys.executable, str((args.ldsc_dir / "ldsc.py").resolve()), "--rg",
           str(out / "height.sumstats.gz") + "," + str(out / "weight.sumstats.gz"),
           "--ref-ld", str(out / "reference"), "--w-ld", str(out / "weights"),
           "--not-M-5-50", "--n-blocks", str(config["blocks"]), "--no-check-alleles", "--out", str(out / "height_weight")]
    (out / "command.json").write_text(json.dumps(cmd, indent=2) + "\n")
    with (out / "stdout.log").open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    text = (out / "height_weight.log").read_text()
    matches = re.findall(r"^Genetic Correlation: ([0-9.eE+-]+) \(([0-9.eE+-]+)\)$", text, re.M)
    if len(matches) != 1:
        raise ValueError("Expected exactly one genetic-correlation result in CLI log")
    cli_rg, cli_se = map(float, matches[0])
    target = None
    for line in (args.results / "ldsc_fits.jsonl").read_text().splitlines():
        row = json.loads(line)
        if (row["specification"], row["i"], row["j"]) == ("unstratified", 9, 10):
            target = row
    if target is None or target["status"] != "succeeded":
        raise ValueError("Run the successful unstratified panel first")
    passed = abs(cli_rg-target["rg"]) <= 0.00005001 and abs(cli_se-target["se"]) <= 0.00005001
    result = {"passed": passed, "n_snps": len(y), "LDSC_commit": commit,
              "cli_rg": cli_rg, "api_rg": target["rg"], "cli_se": cli_se, "api_se": target["se"],
              "tolerance": "0.00005001, reflecting four-decimal CLI report rounding",
              "panel_sha256": sha256(args.panel), "scores_sha256": sha256(args.scores),
              "weights_sha256": sha256(args.weights),
              "alleles": "trusted preharmonized tutorial panel; --no-check-alleles retains same rows"}
    json_write(out / "agreement.json", result)
    print(json.dumps(result, indent=2))
    if not passed:
        raise AssertionError("LDSC command line and API disagree")


if __name__ == "__main__":
    main()
