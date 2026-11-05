#!/usr/bin/env python3
"""Rebuild the released 1 Mb rolling-window L1-L4 score arrays from 1000G EUR.

Only derived scores and provenance are output. The BED/BIM/FAM/AFREQ reference
must remain in a local working directory and is never copied into the release.
The released score file supplies tag metadata and deletion groups, which are
not independent LD blocks. All eligible reference SNPs remain source variants.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import shlex
import sys
import time
import urllib.request

import numpy as np

from build_ldsc_weights import (EDGES, decode_normalized, digest, map_reference,
                                read_tags, reference_files)

TARGET_SHA256 = "bd0b589191c78284cb1c98236ab7e14b7d9d7b34ab6e1d2810958b254d93c1c6"
TARGET_COUNTS = np.array([464696, 205404, 314059, 261563, 473718])


def ensure_afreq(directory, download):
    path = directory / "1000G_EUR.afreq"
    url = "https://zenodo.org/api/records/6614170/files/1000G_EUR.afreq/content"
    size, md5 = 56496658, "7746ab38e42c6abd5980c8ce57172327"
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        part = path.with_suffix(".afreq.part")
        with urllib.request.urlopen(url, timeout=120) as src, part.open("wb") as dst:
            while chunk := src.read(4 * 1024 * 1024):
                dst.write(chunk)
        if part.stat().st_size != size or digest(part, "md5") != md5:
            raise ValueError("AFREQ checksum mismatch")
        part.replace(path)
    if path.stat().st_size != size or digest(path, "md5") != md5:
        raise ValueError("AFREQ checksum mismatch")
    return {"source_url": url, "bytes": size, "md5": md5, "sha256": digest(path)}


def read_targets(reference_dir, exclusions):
    with exclusions.open() as handle:
        excluded = {r["SNP"]: (r["CHR"], r["BP"], r["A1"], r["A2"])
                    for r in csv.DictReader(handle, delimiter="\t")}
    if len(excluded) != 286:
        raise ValueError("Expected the 286 archived frequency-outlier exclusions")
    indices, chroms, positions, frequencies, seen = [], [], [], [], set()
    with (reference_dir / "1000G_EUR.bim").open() as bim, (reference_dir / "1000G_EUR.afreq").open() as af:
        next(af)
        for index, (bl, fl) in enumerate(zip(bim, af, strict=True)):
            c, snp, cm, bp, a1, a2 = bl.split()
            f = fl.split()
            if (snp, a1, a2) != (f[1], f[3], f[2]):
                raise ValueError(f"BIM/AFREQ row mismatch: {snp}")
            p = float(f[4])
            if snp in excluded:
                if (c, bp, a1, a2) != excluded[snp]:
                    raise ValueError(f"Frequency-exclusion mismatch: {snp}")
                seen.add(snp)
                continue
            if not (1 <= int(c) <= 22 and len(a1) == len(a2) == 1 and
                    a1 in "ACGT" and a2 in "ACGT" and a1 != a2 and 0 < p < 1):
                continue
            if c == "6" and 25_000_000 <= int(bp) <= 34_000_000:
                continue
            indices.append(index); chroms.append(c); positions.append(int(bp)); frequencies.append(p)
    if seen != set(excluded):
        raise ValueError("Some frequency exclusions are absent from reference")
    indices = np.array(indices, dtype=np.int64)
    target_hash = hashlib.sha256(indices.tobytes()).hexdigest()
    if target_hash != TARGET_SHA256:
        raise ValueError(f"Source universe differs from original reference: {target_hash}")
    p = np.array(frequencies)
    return dict(indices=indices, chrom=np.array(chroms), bp=np.array(positions),
                af=np.minimum(p, 1 - p), signs=np.where(p <= .5, 1., -1.))


def score_rect(x, bp, bins, row_indices, n_reference, window, tile):
    n, nb = len(row_indices), 5
    left = np.searchsorted(bp, bp[row_indices] - window, side="left")
    right = np.searchsorted(bp, bp[row_indices] + window, side="right")
    cumulative = np.zeros((len(x) + 1, nb), dtype=np.int32)
    cumulative[1:] = np.cumsum(bins[:, None] == np.arange(nb), axis=0)
    counts = cumulative[right] - cumulative[left]
    scores = np.empty((4, n, nb), dtype=np.float64)
    oracle = np.zeros(4)
    for start in range(0, n, tile):
        stop = min(n, start + tile)
        lo, hi = int(left[start]), int(right[stop - 1])
        r = x[row_indices[start:stop]] @ x[lo:hi].T
        np.clip(r, -1, 1, out=r)
        target = np.arange(lo, hi)
        keep = ((target[None, :] >= left[start:stop, None]) &
                (target[None, :] < right[start:stop, None]))
        r[~keep] = 0
        r[np.arange(stop - start), row_indices[start:stop] - lo] = 1
        power = r.copy()
        for q in range(4):
            for b in range(nb):
                scores[q, start:stop, b] = np.sum(power[:, bins[lo:hi] == b], axis=1, dtype=np.float64)
            if q < 3:
                power *= r
        if start == 0:
            j = min(7, stop - 1)
            rr = x[left[j]:right[j]].astype(np.float64) @ x[row_indices[j]].astype(np.float64)
            rr[row_indices[j] - left[j]] = 1
            for q in range(4):
                truth = np.array([np.sum(rr[bins[left[j]:right[j]] == b] ** (q + 1)) for b in range(nb)])
                oracle[q] = np.max(np.abs(scores[q, j] - truth))
                np.testing.assert_allclose(scores[q, j], truth, rtol=3e-5, atol=3e-5)
    scores[1] = ((n_reference - 1) * scores[1] - counts) / (n_reference - 2)
    return scores, oracle


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--frequency-exclusions", type=Path, required=True)
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--out", type=Path, required=True, help="Separate rebuilt score archive, never overwrites the released input")
    args = ap.parse_args()
    if args.out.exists() or args.out.resolve() == args.scores.resolve():
        ap.error("Choose a new output path")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(level=logging.INFO, format="%(asctime)sZ %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(str(args.out) + ".log")])
    start = time.monotonic()
    logging.info("COMMAND %s", shlex.join(sys.argv))
    reference = reference_files(args.reference_dir, args.download)
    reference["1000G_EUR.afreq"] = ensure_afreq(args.reference_dir, args.download)
    from gcmsem.scores import load_scores, save_scores
    scores = load_scores(args.scores)
    tags = read_tags(args.scores)
    tag_indices, signs = map_reference(tags, args.reference_dir / "1000G_EUR.bim")
    targets = read_targets(args.reference_dir, args.frequency_exclusions)
    bins = np.searchsorted(EDGES[1:-1], targets["af"], side="right")
    np.testing.assert_array_equal(np.bincount(bins, minlength=5), TARGET_COUNTS)
    np.testing.assert_array_equal(scores["counts"], TARGET_COUNTS)
    n_reference = sum(1 for _ in (args.reference_dir / "1000G_EUR.fam").open())
    bed = np.memmap(args.reference_dir / "1000G_EUR.bed", mode="r", dtype=np.uint8)
    np.testing.assert_array_equal(bed[:3], [108, 27, 1])
    rebuilt = np.empty((4, len(tags["snp"]), 5), dtype=np.float64)
    details = []
    for chrom in sorted(set(tags["chrom"]), key=int):
        rows = np.flatnonzero(tags["chrom"] == chrom)
        target_rows = np.flatnonzero(targets["chrom"] == chrom)
        positions = np.searchsorted(targets["indices"][target_rows], tag_indices[rows])
        np.testing.assert_array_equal(targets["indices"][target_rows][positions], tag_indices[rows])
        x, af = decode_normalized(bed, targets["indices"][target_rows],
                                  targets["signs"][target_rows], n_reference)
        af_error = float(np.max(np.abs(af - targets["af"][target_rows])))
        if af_error > 1e-5:
            raise ValueError(f"Reference AF mismatch on chromosome {chrom}")
        result, oracle = score_rect(x, targets["bp"][target_rows], bins[target_rows],
                                    positions, n_reference, 1_000_000, 512)
        # Full array comparison against the released historical computations.
        differences = []
        for q in range(4):
            old = scores[f"L{q + 1}"][rows]
            np.testing.assert_allclose(result[q], old, rtol=3e-5, atol=3e-5)
            differences.append(float(np.max(np.abs(result[q] - old))))
        rebuilt[:, rows] = result
        details.append(dict(chrom=chrom, tags=len(rows), sources=len(target_rows),
                            oracle_max_abs_error_by_order=oracle.tolist(),
                            release_max_abs_error_by_order=differences, max_af_error=af_error))
        logging.info("CHROM %s tags=%d sources=%d release_max_error=%s", chrom, len(rows), len(target_rows), differences)
    audit = {
        "reference": "https://zenodo.org/records/6614170", "reference_files": reference,
        "n_reference": n_reference, "tag_count": len(tags["snp"]),
        "source_count": len(targets["indices"]), "source_bin_counts": TARGET_COUNTS.tolist(),
        "target_bed_index_sha256": TARGET_SHA256,
        "frequency_exclusions_sha256": digest(args.frequency_exclusions),
        "input_scores_sha256": digest(args.scores),
        "window_bp": 1_000_000, "cross_chromosome_ld": False,
        "finite_reference_corrected_orders": [2], "uncorrected_orders": [1, 3, 4],
        "all_score_elements_compared": int(rebuilt.size),
        "all_score_elements_within_tolerance": True,
        "comparison_atol": 3e-5, "comparison_rtol": 3e-5,
        "max_abs_difference_by_order": np.max([d["release_max_abs_error_by_order"] for d in details], axis=0).tolist(),
        "raw_genotypes_in_output": False, "chromosomes": details,
        "elapsed_seconds": time.monotonic() - start,
    }
    for q in range(4):
        scores[f"L{q + 1}"] = rebuilt[q]
    scores["metadata"]["rebuild_audit"] = audit
    save_scores(scores, args.out)
    audit["rebuilt_scores_sha256"] = digest(args.out)
    args.out.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    logging.info("SUCCESS all_elements=%d max_abs_difference=%s seconds=%.1f", rebuilt.size,
                 audit["max_abs_difference_by_order"], audit["elapsed_seconds"])


if __name__ == "__main__":
    main()
