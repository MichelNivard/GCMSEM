#!/usr/bin/env python3
"""Re-estimate real-data cumulants and fit the two causal tutorial examples.

Run after downloading the separately hosted tutorial data and installing both
GCMSEM packages. No inputs, large outputs or reference genotypes enter git.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--starts", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument("--blocks", type=int, default=200)
    parser.add_argument("--skip-moments", action="store_true",
                        help="reuse completed moment outputs at --out; their checksums are recorded")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.starts < 1 or args.blocks < 2:
        parser.error("starts must be positive and blocks must be at least two")
    data = args.data_dir.resolve()
    out = (args.out or data / "results" / "causal").resolve()
    score_path = data / "scores" / "eur_hm3.scores.npz"
    r_script = Path(__file__).resolve().with_name("fit_causal.R")
    if not score_path.is_file():
        parser.error(f"missing downloaded score archive: {score_path}")
    pairs = ("height_weight", "bmi_sbp")
    # Reject accidental reruns before touching an existing completion record.
    for pair in pairs:
        manifest = data / "manifests" / f"{pair}.tsv"
        if not manifest.is_file():
            parser.error(f"missing downloaded manifest: {manifest}")
        if (out / pair).exists() and not (args.overwrite or args.skip_moments):
            parser.error(f"output exists: {out / pair}; select a new --out or --overwrite")
        with manifest.open() as handle:
            for entry in csv.DictReader(handle, delimiter="\t"):
                path = (manifest.parent / entry["path"]).resolve()
                if not path.is_file():
                    parser.error(f"missing downloaded GWAS: {path}")
    provenance = {"started_utc": datetime.now(timezone.utc).isoformat(),
                  "status": "running", "score_sha256": sha256(score_path),
                  "commands": [], "pairs": {}}

    def validate_bundle(pair: str) -> list[Path]:
        manifest = data / "manifests" / f"{pair}.tsv"
        prefix = out / pair / "moments"
        required = [Path(str(prefix) + suffix) for suffix in
                    (".moments.tsv", ".sampling_cov.tsv", ".metadata.json")]
        if not all(p.is_file() for p in required):
            parser.error(f"incomplete moment bundle at {prefix}")
        metadata = json.loads(required[2].read_text())
        if metadata.get("score_sha256") != provenance["score_sha256"] or \
                metadata.get("sumstats_manifest_sha256") != sha256(manifest):
            parser.error(f"bundle input hashes do not match current inputs: {prefix}")
        previous_inputs = {entry["trait"]: entry["sha256"] for entry in
                           metadata.get("qc", {}).get("input_files", [])}
        with manifest.open() as handle:
            for entry in csv.DictReader(handle, delimiter="\t"):
                path = (manifest.parent / entry["path"]).resolve()
                if previous_inputs.get(entry["trait"]) != sha256(path):
                    parser.error(f"bundle GWAS hash mismatch: {entry['trait']}")
        return required

    if args.skip_moments:
        for pair in pairs:
            validate_bundle(pair)
    out.mkdir(parents=True, exist_ok=True)
    def record() -> None:
        temporary = out / "execution.json.tmp"
        temporary.write_text(json.dumps(provenance, indent=2) + "\n")
        temporary.replace(out / "execution.json")

    def run(command: list[str]) -> None:
        print("Running:", " ".join(command), flush=True)
        provenance["commands"].append(command)
        record()
        try:
            subprocess.run(command, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            provenance.update(status="failed", failed_command=command, error=str(exc))
            record()
            raise

    record()
    for pair in pairs:
        manifest = data / "manifests" / f"{pair}.tsv"
        pair_out = out / pair
        prefix = pair_out / "moments"
        pair_out.mkdir(parents=True, exist_ok=True)
        estimate = [sys.executable, "-m", "gcmsem", "estimate", "--sumstats", str(manifest),
                    "--scores", str(score_path), "--out", str(prefix), "--max-order", "4",
                    "--blocks", str(args.blocks)]
        if args.overwrite:
            estimate.append("--overwrite")
        if not args.skip_moments:
            run(estimate)
        required = validate_bundle(pair)
        fit = [args.rscript, str(r_script), str(prefix), str(pair_out / "sem"),
               str(args.starts), str(args.seed)]
        run(fit)
        provenance["pairs"][pair] = {"manifest_sha256": sha256(manifest),
            "moment_bundle_sha256": {p.name: sha256(p) for p in required}}
        record()
    provenance["status"] = "completed"
    provenance["completed_utc"] = datetime.now(timezone.utc).isoformat()
    record()
    print(f"Completed both examples: {out}")


if __name__ == "__main__":
    main()
