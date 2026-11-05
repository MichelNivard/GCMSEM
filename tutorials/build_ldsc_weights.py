#!/usr/bin/env python3
"""Recompute tutorial HM3-only LDSC weights; never package the reference genotypes.

The input is the released GCMSEM score archive (only its tag metadata are used).
For exact historical reconstruction, --aligned accepts the original prealigned
height/weight archive instead. Both routes use the same SNP set as both tag and
source population, a same-chromosome +/-1 Mb window, and the LDSC finite-sample
r-squared correction. These are weighting scores, not all-target moment scores.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import platform
import shlex
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

REFERENCE = "https://zenodo.org/records/6614170"
FILES = {
    "1000G_EUR.bed": (231387159, "a163c74ea92676c9ada76fe5eb3b2ec9"),
    "1000G_EUR.bim": (51537337, "81b1ee40a05ef5be98ad3e496c71d708"),
    "1000G_EUR.fam": (9557, "669a4260fda7a9e1dd7df374aff294ea"),
}
EDGES = np.array([0, .05, .1, .2, .3, .500000001])


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def reference_files(directory, download):
    directory.mkdir(parents=True, exist_ok=True)
    provenance = {}
    for name, (size, md5) in FILES.items():
        path = directory / name
        url = f"https://zenodo.org/api/records/6614170/files/{name}/content"
        if not path.exists():
            if not download:
                raise FileNotFoundError(f"Missing {path}; obtain the reference or use --download")
            logging.info("DOWNLOAD %s", name)
            part = path.with_suffix(path.suffix + ".part")
            with urllib.request.urlopen(url, timeout=120) as src, part.open("wb") as dst:
                while chunk := src.read(4 * 1024 * 1024):
                    dst.write(chunk)
            if part.stat().st_size != size or digest(part, "md5") != md5:
                raise ValueError(f"Reference checksum mismatch: {name}")
            part.replace(path)
        if path.stat().st_size != size or digest(path, "md5") != md5:
            raise ValueError(f"Reference checksum mismatch: {name}")
        provenance[name] = {"source_url": url, "bytes": size, "md5": md5,
                            "sha256": digest(path)}
        logging.info("REFERENCE verified %s", name)
    return provenance


def read_tags(path, legacy=False):
    with np.load(path, allow_pickle=False) as data:
        if legacy:
            snp = data["rsid"]
            chrom = data["chr"].astype(str)
            bp = data["bp"]
            p = data["p_a1_ref"]
            a1 = data["minor_allele"]
            a2 = np.where(a1 == data["a1"], data["a2"], data["a1"])
            af = data["minor_af_ref"]
            indices = data["bed_index"]
        else:
            snp, chrom, bp, a1, a2, af = (
                data[key] for key in ("snp", "chrom", "bp", "a1", "a2", "af"))
            chrom = chrom.astype(str)
            indices = None
    if len(set(snp)) != len(snp) or not np.all((af > 0) & (af <= .5)):
        raise ValueError("Input tag identifiers must be unique and frequencies in (0,.5]")
    return dict(snp=snp, chrom=chrom, bp=bp, a1=a1, a2=a2, af=af,
                legacy_indices=indices)


def map_reference(tags, bim):
    wanted = {str(snp): i for i, snp in enumerate(tags["snp"])}
    indices = np.full(len(wanted), -1, dtype=np.int64)
    signs = np.empty(len(wanted), dtype=np.float64)
    with bim.open() as handle:
        for index, line in enumerate(handle):
            chrom, snp, cm, bp, a1, a2 = line.split()
            row = wanted.get(snp)
            if row is None:
                continue
            if indices[row] != -1:
                raise ValueError(f"Duplicate reference tag: {snp}")
            if chrom != str(tags["chrom"][row]) or int(bp) != tags["bp"][row]:
                raise ValueError(f"Reference coordinate mismatch: {snp}")
            pair = (str(tags["a1"][row]), str(tags["a2"][row]))
            if pair == (a1, a2):
                signs[row] = 1
            elif pair == (a2, a1):
                signs[row] = -1
            else:
                raise ValueError(f"Reference allele mismatch: {snp}")
            indices[row] = index
    if np.any(indices < 0):
        raise ValueError(f"{np.sum(indices < 0)} tags absent from the reference")
    if tags["legacy_indices"] is not None:
        np.testing.assert_array_equal(indices, tags["legacy_indices"])
    return indices, signs


def decode_normalized(bed, indices, signs, n_reference):
    stride = (n_reference + 3) // 4
    raw = bed[3:].reshape(-1, stride)[indices]
    codes = ((raw[:, :, None] >> (2 * np.arange(4, dtype=np.uint8))) & 3)
    codes = codes.reshape(len(raw), stride * 4)[:, :n_reference]
    x = np.array([2, np.nan, 1, 0], dtype=np.float32)[codes]
    mean = np.nanmean(x, axis=1, dtype=np.float64)
    x -= mean.astype(np.float32)[:, None]
    np.nan_to_num(x, copy=False, nan=0)
    norm = np.sqrt(np.sum(x * x, axis=1, dtype=np.float64))
    if not np.all(norm > 0):
        raise ValueError("Monomorphic/all-missing tag in the reference")
    x *= (signs / norm).astype(np.float32)[:, None]
    return np.ascontiguousarray(x), np.where(signs == 1, mean / 2, 1 - mean / 2)


def score_chromosome(x, bp, bins, n_reference, window, tile):
    n, nb = len(x), len(EDGES) - 1
    left = np.searchsorted(bp, bp - window, side="left")
    right = np.searchsorted(bp, bp + window, side="right")
    counts = np.zeros((n + 1, nb), dtype=np.int32)
    counts[1:] = np.cumsum(bins[:, None] == np.arange(nb), axis=0)
    window_counts = counts[right] - counts[left]
    corrected = np.empty((n, nb), dtype=np.float64)
    raw_totals = np.zeros(nb)
    oracle_max_error = 0.
    for start in range(0, n, tile):
        stop = min(start + tile, n)
        lo, hi = int(left[start]), int(right[stop - 1])
        r = x[start:stop] @ x[lo:hi].T
        np.clip(r, -1, 1, out=r)
        targets = np.arange(lo, hi)
        keep = ((targets[None, :] >= left[start:stop, None]) &
                (targets[None, :] < right[start:stop, None]))
        r[~keep] = 0
        r[np.arange(stop - start), np.arange(start, stop) - lo] = 1
        square = r * r
        raw = np.stack([np.sum(square[:, bins[lo:hi] == b], axis=1,
                               dtype=np.float64) for b in range(nb)], axis=1)
        corrected[start:stop] = ((n_reference - 1) * raw -
                                 window_counts[start:stop]) / (n_reference - 2)
        raw_totals += raw.sum(axis=0)
        # Independent float64 dot-product oracle for a tag in the first tile.
        if start == 0:
            j = min(7, stop - 1)
            rr = x[left[j]:right[j]].astype(np.float64) @ x[j].astype(np.float64)
            rr[j - left[j]] = 1
            truth = np.array([np.sum(rr[bins[left[j]:right[j]] == b] ** 2)
                              for b in range(nb)])
            oracle_max_error = float(np.max(np.abs(truth - raw[j])))
            np.testing.assert_allclose(raw[j], truth, rtol=2e-5, atol=2e-5)
    return corrected, raw_totals, oracle_max_error


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--scores", type=Path, help="Released GCMSEM score archive")
    inputs.add_argument("--aligned", type=Path, help="Historical aligned height/weight archive")
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    parser.add_argument("--tile", type=int, default=512)
    args = parser.parse_args(argv)
    if args.out.exists():
        parser.error("Output exists; choose a new path")
    if args.window_bp < 0 or args.tile < 1:
        parser.error("window-bp must be nonnegative and tile positive")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    handlers = [logging.StreamHandler(), logging.FileHandler(str(args.out) + ".log")]
    logging.Formatter.converter = time.gmtime
    logging.basicConfig(level=logging.INFO, handlers=handlers,
                        format="%(asctime)sZ %(levelname)s %(message)s")
    start = time.monotonic()
    logging.info("COMMAND %s", shlex.join(sys.argv))
    provenance = reference_files(args.reference_dir, args.download)
    input_path = args.scores or args.aligned
    tags = read_tags(input_path, legacy=args.aligned is not None)
    indices, signs = map_reference(tags, args.reference_dir / "1000G_EUR.bim")
    n_reference = sum(1 for _ in (args.reference_dir / "1000G_EUR.fam").open())
    bed = np.memmap(args.reference_dir / "1000G_EUR.bed", mode="r", dtype=np.uint8)
    np.testing.assert_array_equal(bed[:3], [108, 27, 1])
    bins = np.searchsorted(EDGES[1:-1], tags["af"], side="right")
    result = np.empty((len(indices), 5), dtype=np.float64)
    details = []
    for chrom in sorted(set(tags["chrom"]), key=int):
        loc = np.flatnonzero(tags["chrom"] == chrom)
        loc = loc[np.argsort(tags["bp"][loc], kind="stable")]
        x, af = decode_normalized(bed, indices[loc], signs[loc], n_reference)
        error = float(np.max(np.abs(af - tags["af"][loc])))
        if error > 1e-5:
            raise ValueError(f"Decoded AF mismatch on chromosome {chrom}: {error}")
        q2, raw_sum, oracle_error = score_chromosome(
            x, tags["bp"][loc], bins[loc], n_reference, args.window_bp, args.tile)
        result[loc] = q2
        detail = dict(chrom=chrom, rows=len(loc), max_af_error=error,
                      oracle_max_abs_error=oracle_error, raw_q2_sum_by_bin=raw_sum.tolist())
        details.append(detail)
        logging.info("CHROM %s rows=%s mean_L2=%.10g oracle_error=%.3g",
                     chrom, len(loc), q2.sum(axis=1).mean(), oracle_error)
    metadata = {
        "schema_version": "1.0", "kind": "LDSC_HM3_weighting_scores",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "reference": REFERENCE, "reference_doi": "10.5281/zenodo.6614170",
        "reference_creator": "Aw, Alan", "reference_license": "CC-BY-4.0",
        "reference_files": provenance, "n_reference": n_reference,
        "input_sha256": digest(input_path), "input_type": "scores" if args.scores else "historical_aligned",
        "tag_count": len(indices), "source_count": len(indices),
        "source_definition": "Identical retained primary HM3 tags; no other source variants",
        "window_bp": args.window_bp, "same_chromosome_only": True,
        "maf_edges": [0, .05, .1, .2, .3, .5],
        "maf_counts": np.bincount(bins, minlength=5).tolist(),
        "normalization": "Mean-imputed centered dosages normalized to unit empirical Euclidean norm; float32 correlation products, float64 sums",
        "finite_reference_correction": "((n_reference-1)*sum(r^2)-window_variant_count)/(n_reference-2); self term one",
        "negative_values": "Retained; the downstream unmodified LDSC weighting code applies its own floor",
        "mean_by_bin": result.mean(axis=0).tolist(),
        "mean_total": float(result.sum(axis=1).mean()),
        "chromosomes": details, "raw_genotypes_in_output": False,
        "python_version": platform.python_version(), "numpy_version": np.__version__,
        "elapsed_seconds": time.monotonic() - start,
    }
    with args.out.open("wb") as handle:
        np.savez_compressed(handle, snp=tags["snp"], chrom=tags["chrom"],
                            bp=tags["bp"], L2=result.sum(axis=1),
                            L2_by_maf=result, metadata=np.array(json.dumps(metadata)))
    args.out.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    logging.info("SUCCESS tags=%d mean_L2=%.12g sha256=%s seconds=%.1f",
                 len(indices), metadata["mean_total"], digest(args.out), metadata["elapsed_seconds"])


if __name__ == "__main__":
    main()
