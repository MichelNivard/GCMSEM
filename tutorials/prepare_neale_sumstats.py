#!/usr/bin/env python3
"""Stream Neale round-2 GWAS results into GCMSEM's per-allele input schema.

Match GRCh37 position AND alleles to the tags in an existing signed-score file.
The score archive and its reference source-variant counts are never changed.

Neale's own LDSC export assigns A1=ALT and A2=REF to the unchanged signed t
statistic: https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/27.export_ldsc_sumstats.py
The beta/se/tstat columns come from the same regression results in
0.1/23.export_results.py. minor_allele and minor_AF describe frequency, and do
not change which allele beta refers to. The output manifest must therefore
specify effect_scale='allele'.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shlex
import sys
import tempfile
import time

from gcmsem.scores import load_scores


REQUIRED = ("variant", "minor_allele", "minor_AF", "low_confidence_variant",
            "n_complete_samples", "beta", "se")
OUTPUT = ("SNP", "A1", "A2", "BETA", "SE", "N", "EAF", "CHR", "BP")
CHROMOSOMES = frozenset(str(i) for i in range(1, 23)) | {"X", "Y", "XY", "MT", "M"}
BASES = frozenset("ACGT")
ORIENTATION_SOURCE = (
    "https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/27.export_ldsc_sumstats.py"
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def chromosome(value):
    value = str(value).upper()
    if value.startswith("CHR"):
        value = value[3:]
    if value.isdigit():
        value = str(int(value))
    if value not in CHROMOSOMES:
        raise ValueError("Invalid chromosome")
    return value


def variant_key(chrom, bp, a1, a2):
    chrom = chromosome(chrom)
    position = int(bp)
    a1, a2 = str(a1).upper(), str(a2).upper()
    if position <= 0 or a1 not in BASES or a2 not in BASES or a1 == a2:
        raise ValueError("Expected a positive position and two different SNP alleles")
    return chrom, position, "".join(sorted((a1, a2)))


def score_mapping(scores):
    mapping = {}
    for snp, chrom, bp, a1, a2 in zip(
            scores["snp"], scores["chrom"], scores["bp"], scores["a1"], scores["a2"], strict=True):
        key = variant_key(chrom, bp, a1, a2)
        if key in mapping:
            raise ValueError(f"Ambiguous score mapping for {key}: {mapping[key]} and {snp}")
        mapping[key] = str(snp)
    return mapping


def open_input(path):
    # BGZF files are concatenated gzip members; gzip reads and verifies them all.
    if str(path).lower().endswith((".gz", ".bgz", ".bgzf")):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return Path(path).open(encoding="utf-8", newline="")


def convert_rows(source, destination, mapping, log):
    reader = csv.reader(source, delimiter="\t")
    header = next(reader, None)
    if not header or len(header) != len(set(header)):
        raise ValueError("Neale input needs a nonempty header with unique column names")
    missing = set(REQUIRED) - set(header)
    if missing:
        raise ValueError(f"Missing Neale columns: {', '.join(sorted(missing))}")
    indices = {name: header.index(name) for name in REQUIRED}
    counts = dict.fromkeys(("input_rows", "unsupported_variant", "not_score_tag",
                           "matched_score_tags", "low_confidence", "invalid_confidence_flag",
                           "invalid_numeric", "invalid_minor_allele", "invalid_minor_AF",
                           "retained", "eaf_complemented"), 0)
    writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
    writer.writerow(OUTPUT)
    seen = set()
    for line_number, row in enumerate(reader, 2):
        counts["input_rows"] += 1
        if len(row) != len(header):
            raise ValueError(f"Malformed Neale input at line {line_number}: expected {len(header)} fields, got {len(row)}")
        parts = row[indices["variant"]].split(":")
        try:
            if len(parts) != 4:
                raise ValueError("Expected CHR:POS:REF:ALT")
            key = variant_key(*parts)
        except ValueError:
            counts["unsupported_variant"] += 1
            continue
        snp = mapping.get(key)
        if snp is None:
            counts["not_score_tag"] += 1
            continue
        if snp in seen:
            raise ValueError(f"Duplicate matched SNP {snp} at input line {line_number}")
        seen.add(snp)
        counts["matched_score_tags"] += 1
        confidence = row[indices["low_confidence_variant"]].lower()
        if confidence == "true":
            counts["low_confidence"] += 1
            continue
        if confidence != "false":
            counts["invalid_confidence_flag"] += 1
            continue
        try:
            beta, se, n = (float(row[indices[name]]) for name in ("beta", "se", "n_complete_samples"))
            if not all(map(math.isfinite, (beta, se, n))) or se <= 0 or n <= 0 or not n.is_integer():
                raise ValueError("Invalid beta, standard error, or sample size")
        except ValueError:
            counts["invalid_numeric"] += 1
            continue
        ref, alt = parts[2].upper(), parts[3].upper()
        minor = row[indices["minor_allele"]].upper()
        if minor not in (ref, alt):
            counts["invalid_minor_allele"] += 1
            continue
        try:
            maf = float(row[indices["minor_AF"]])
            if not math.isfinite(maf) or not 0 < maf <= 0.5:
                raise ValueError("Invalid minor-allele frequency")
        except ValueError:
            counts["invalid_minor_AF"] += 1
            continue
        eaf = maf if minor == alt else 1 - maf
        counts["eaf_complemented"] += int(minor == ref)
        writer.writerow((snp, alt, ref, row[indices["beta"]], row[indices["se"]],
                         int(n), format(eaf, ".15g"), key[0], key[1]))
        counts["retained"] += 1
        if counts["retained"] % 250000 == 0:
            log(f"Retained {counts['retained']:,} tags after reading {counts['input_rows']:,} rows")
    counts["score_tags_absent"] = len(mapping) - len(seen)
    if not counts["retained"]:
        raise ValueError("No usable score-tag SNPs remained after preparation")
    return counts


def prepare(sumstats, scores, out, *, overwrite=False, command=None):
    """Write a gzip TSV and adjacent .json/.log, returning preparation metadata.

    Complete streaming/CRC validation precedes replacing any existing output.
    File replacements are atomic individually; the TSV is published last.
    """
    started = time.monotonic()
    sumstats, scores, out = (Path(p).resolve() for p in (sumstats, scores, out))
    if not str(out).endswith(".gz"):
        raise ValueError("--out must end in .gz (for example height.tsv.gz)")
    outputs = [out, Path(str(out) + ".json"), Path(str(out) + ".log")]
    if sumstats in outputs or scores in outputs:
        raise ValueError("Output paths must differ from input paths")
    for path in (sumstats, scores):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in outputs:
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output exists: {path}; use --overwrite to replace it")
    out.parent.mkdir(parents=True, exist_ok=True)
    messages = []

    def log(message):
        messages.append(message)
        print(message, file=sys.stderr, flush=True)

    log("GCMSEM Neale GWAS preparation")
    log(f"Input: {sumstats}")
    log(f"Scores: {scores}")
    archive = load_scores(scores)
    mapping = score_mapping(archive)
    source_variants = int(archive["metadata"]["source_variants"])
    del archive
    log(f"Score tags: {len(mapping):,}; reference source variants: {source_variants:,}")
    log("BETA/SE scale: per ALT allele; EAF is the ALT-allele frequency")
    with tempfile.TemporaryDirectory(prefix=".neale-", dir=out.parent) as temporary:
        stage = Path(temporary)
        data_path = stage / "sumstats.tsv.gz"
        with open_input(sumstats) as source, data_path.open("wb") as binary:
            with gzip.GzipFile(fileobj=binary, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as destination:
                    counts = convert_rows(source, destination, mapping, log)
        # Reading to EOF checks the complete gzip stream and its CRC.
        with gzip.open(data_path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        log(f"Input rows: {counts['input_rows']:,}; matched score tags: {counts['matched_score_tags']:,}")
        for name, value in counts.items():
            if name not in ("input_rows", "matched_score_tags", "retained"):
                log(f"{name}: {value:,}")
        log(f"Retained SNPs: {counts['retained']:,}")
        metadata = {
            "schema_version": "1.0", "kind": "neale_round2_gcmsem_sumstats",
            "effect_scale": "allele", "effect_allele": "ALT", "other_allele": "REF",
            "frequency": "ALT frequency derived from minor_allele and minor_AF",
            "orientation_source": ORIENTATION_SOURCE,
            "genome_build": "GRCh37", "columns": list(OUTPUT),
            "input": {"path": str(sumstats), "sha256": sha256(sumstats)},
            "scores": {"path": str(scores), "sha256": sha256(scores),
                       "tag_variants": len(mapping), "source_variants": source_variants},
            "output": {"path": str(out), "sha256": sha256(data_path)},
            "counts": counts, "elapsed_seconds": round(time.monotonic() - started, 3),
            "command": command,
            "next_step": "Set effect_scale=allele in the manifest. estimate_moments performs allele harmonization, frequency QC, palindrome filtering and common-trait SNP selection.",
        }
        log(f"Completed in {metadata['elapsed_seconds']:.1f} seconds")
        log(f"Output: {out}")
        (stage / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
        (stage / "prepare.log").write_text("\n".join(messages) + "\n")
        # Check again before publication. No output is touched on a parse/QC error.
        for path in outputs:
            if path.exists() and not overwrite:
                raise FileExistsError(f"Output appeared during processing: {path}")
        for source, destination in ((stage / "metadata.json", outputs[1]),
                                    (stage / "prepare.log", outputs[2]), (data_path, out)):
            if overwrite:
                os.replace(source, destination)
            else:
                # Hard linking within this directory avoids clobbering a racing writer.
                os.link(source, destination)
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sumstats", required=True, help="Official Neale .tsv.bgz (or uncompressed TSV)")
    parser.add_argument("--scores", required=True, help="Existing GCMSEM signed-score NPZ archive")
    parser.add_argument("--out", required=True, help="Prepared TSV gzip; .json and .log are written alongside")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output after successful validation")
    args = parser.parse_args(argv)
    try:
        prepare(args.sumstats, args.scores, args.out, overwrite=args.overwrite,
                command=shlex.join(sys.argv if argv is None else [parser.prog, *argv]))
    except (OSError, ValueError, EOFError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
