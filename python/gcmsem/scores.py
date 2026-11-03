"""Signed LD-power scores, reference-panel QC, and portable score archives.

For disjoint source strata c, Lq[j,c] = sum(k in c) R[j,k]**q.
R is signed LD between genotype-standardized, minor-allele-coded genotypes.
Blocks are an explicit independence/truncation assumption; this module never
silently invents genomic blocks. A regression/tag subset does not subset the
source universe. Third and fourth powers are conditional on the supplied LD,
not finite-reference unbiased estimates of population LD powers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import platform
import tempfile
from typing import Iterable, Mapping
import warnings

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

SCHEMA_VERSION = "1.0"
_REQUIRED = ("snp", "a1", "a2", "af", "chrom", "bp", "block")


def _log(log, message, level="info"):
    if log is None:
        getattr(logging.getLogger("gcmsem"), level)(message)
    elif hasattr(log, level):
        getattr(log, level)(message)
    else:
        log(message)


def _warn(log, message):
    _log(log, message, "warning")
    if log is None and not logging.getLogger("gcmsem").hasHandlers():
        warnings.warn(message, RuntimeWarning, stacklevel=3)


def _strings(value, n, name):
    a = np.asarray(value).astype(str)
    if a.ndim == 0:
        a = np.repeat(a, n)
    if a.shape != (n,) or np.any(np.char.strip(a) == ""):
        raise ValueError(f"{name} must contain {n} nonempty entries")
    if np.any([any(c.isspace() for c in x) for x in a]):
        raise ValueError(f"{name} must not contain whitespace")
    return a


def _metadata(block, n):
    missing = [k for k in ("snp", "a1", "a2", "chrom", "bp") if k not in block]
    if missing:
        raise ValueError(f"Block is missing SNP metadata: {', '.join(missing)}")
    out = {k: _strings(block[k], n, k) for k in ("snp", "a1", "a2", "chrom")}
    out["a1"] = np.char.upper(out["a1"])
    out["a2"] = np.char.upper(out["a2"])
    if not np.all(np.isin(out["a1"], list("ACGT"))) or not np.all(np.isin(out["a2"], list("ACGT"))):
        raise ValueError("Only biallelic A/C/G/T SNP alleles are supported")
    if np.any(out["a1"] == out["a2"]):
        raise ValueError("Reference A1 and A2 must differ")
    bp = np.asarray(block["bp"], dtype=float)
    if bp.shape != (n,) or not np.isfinite(bp).all() or np.any(bp <= 0) or np.any(bp != np.floor(bp)):
        raise ValueError("bp must contain positive integer base-pair positions")
    out["bp"] = bp.astype(np.int64)
    if len(set(out["snp"])) != n:
        raise ValueError("Duplicate SNP identifiers within reference block")
    if len(set(out["chrom"])) != 1:
        raise ValueError("Each LD block must lie on one chromosome")
    return out


def _prepare_block(block, check_psd, min_maf, max_variants, log):
    if ("ld" in block) == ("genotype" in block):
        raise ValueError("Each block needs exactly one of ld or genotype")
    if "ld" in block:
        r = block["ld"]
        if not hasattr(r, "shape"):
            r = np.asarray(r)
        if len(r.shape) == 2 and max(r.shape) > max_variants:
            raise ValueError(f"LD block exceeds max_variants={max_variants}")
        if sparse.issparse(r):
            r = r.toarray()
        r = np.asarray(r, dtype=float)
        if r.ndim != 2 or r.shape[0] != r.shape[1] or r.shape[0] == 0:
            raise ValueError("LD must be a nonempty square matrix")
        n = r.shape[0]
        if n > max_variants:
            raise ValueError(f"LD block has {n} variants; exceeds max_variants={max_variants}")
        meta = _metadata(block, n)
        if "af" not in block:
            raise ValueError("Precomputed LD requires A1 allele frequencies af")
        af = np.asarray(block["af"], dtype=float)
        if af.shape != (n,) or not np.isfinite(af).all() or np.any((af <= 0) | (af >= 1)):
            raise ValueError("Precomputed LD af must be finite and strictly between zero and one")
        if not np.isfinite(r).all() or not np.allclose(r, r.T, rtol=0, atol=1e-8):
            raise ValueError("LD must be finite and symmetric")
        if not np.allclose(np.diag(r), 1, rtol=0, atol=1e-8) or np.max(np.abs(r)) > 1 + 1e-8:
            raise ValueError("LD must be a correlation matrix with unit diagonal and |r| <= 1")
        if check_psd and eigh(r, subset_by_index=[0, 0], check_finite=False, eigvals_only=True)[0] < -1e-7:
            raise ValueError("LD matrix is not positive semidefinite")
        n_reference = block.get("n_reference")
        qc = {"input_variants": n, "all_missing": 0, "monomorphic": 0, "missing_dosages": 0}
    else:
        g = np.asarray(block["genotype"], dtype=float)
        if g.ndim != 2 or min(g.shape) == 0 or g.shape[0] < 3:
            raise ValueError("genotype must be samples x variants, with at least three samples")
        n = g.shape[1]
        if n > max_variants:
            raise ValueError(f"Genotype block has {n} variants; exceeds max_variants={max_variants}")
        if np.any(np.isinf(g)) or np.any((g[np.isfinite(g)] < 0) | (g[np.isfinite(g)] > 2)):
            raise ValueError("Genotypes must be dosages in [0,2], with NaN for missing")
        meta = _metadata(block, n)
        finite = np.isfinite(g)
        n_observed = finite.sum(axis=0)
        means = np.divide(np.nansum(g, axis=0), n_observed, out=np.zeros(n), where=n_observed > 0)
        af = means / 2
        centered = np.where(finite, g - means, 0.0)
        norm = np.sqrt(np.sum(centered ** 2, axis=0))
        keep = (n_observed > 0) & (norm > np.finfo(float).eps) & (af > 0) & (af < 1)
        qc = {"input_variants": n, "all_missing": int(np.sum(n_observed == 0)),
              "monomorphic": int(np.sum((n_observed > 0) & ~keep)),
              "missing_dosages": int(np.sum(~finite))}
        if not np.any(keep):
            raise ValueError("Reference block has no polymorphic, nonmissing variants")
        if "af" in block:
            supplied = np.asarray(block["af"], dtype=float)
            if supplied.shape != (n,) or not np.isfinite(supplied[keep]).all() or np.max(np.abs(supplied[keep] - af[keep])) > 0.02:
                raise ValueError("Supplied reference AF disagrees with genotype A1 dosage frequency by >0.02")
        centered = centered[:, keep] / norm[keep]
        r = centered.T @ centered
        r = (r + r.T) / 2
        np.fill_diagonal(r, 1.0)
        meta = {k: v[keep] for k, v in meta.items()}
        af = af[keep]
        n_reference = g.shape[0]
        if block.get("n_reference") is not None and int(block["n_reference"]) != n_reference:
            raise ValueError("n_reference must match genotype sample count")
    # A1 is always the minor allele; do this before signed odd-power scores.
    # At exactly 0.5, use lexical allele order to make recoding deterministic.
    flip = (af > 0.5) | ((af == 0.5) & (meta["a1"] > meta["a2"]))
    sign = np.where(flip, -1.0, 1.0)
    r = r * sign[:, None] * sign[None, :]
    a1, a2 = meta["a1"].copy(), meta["a2"].copy()
    meta["a1"] = np.where(flip, a2, a1)
    meta["a2"] = np.where(flip, a1, a2)
    af = np.where(flip, 1.0 - af, af)
    keep = af >= min_maf - 1e-12
    qc["minor_allele_flips"] = int(np.sum(flip))
    qc["below_min_maf"] = int(np.sum(~keep))
    if not np.any(keep):
        raise ValueError("Reference block has no variants after min_maf filtering")
    meta = {k: v[keep] for k, v in meta.items()}
    meta["af"] = af[keep]
    r = r[np.ix_(keep, keep)]
    return r, meta, n_reference, qc


def build_scores(blocks: Iterable[Mapping], maf_edges=(0.0, 0.5), tag_snps=None,
                 min_maf=0.0, finite_reference=False, max_distance_bp=None,
                 check_psd=True, max_variants=10000, log=None):
    """Build signed LD-power scores, consuming one independent block at a time.

    Each block is a mapping with ``id`` (or scalar ``block``), ``snp``, ``a1``,
    ``a2``, ``chrom``, ``bp``, and either ``ld`` plus ``af``, or ``genotype``
    (samples x variants, A1 dosage; NaN missing). ``n_reference`` is mandatory
    with precomputed LD and finite_reference=True. AF always refers to A1.
    All metadata arrays must follow matrix column order; no reordering is guessed.

    ``maf_edges`` defines disjoint [left,right) source strata (last bin closed).
    ``tag_snps`` only filters output rows. Source strata/counts retain all QC-passed
    reference variants. ``finite_reference`` applies r²-(1-r²)/(n-2) off diagonal,
    retaining self-LD=1, with no clipping; it does not debias q=1,3,4 scores.
    Return arrays and a JSON-serializable metadata dictionary; use save_scores.
    """
    edges = np.asarray(maf_edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or edges[0] != 0 or edges[-1] != 0.5 or np.any(np.diff(edges) <= 0):
        raise ValueError("maf_edges must increase strictly from 0 to 0.5")
    if not np.isfinite(min_maf) or not 0 <= min_maf < 0.5:
        raise ValueError("min_maf must be in [0,0.5)")
    if max_distance_bp is not None and (not np.isfinite(max_distance_bp) or max_distance_bp < 0):
        raise ValueError("max_distance_bp must be nonnegative")
    if int(max_variants) != max_variants or max_variants < 1:
        raise ValueError("max_variants must be a positive integer")
    tags = None if tag_snps is None else set(str(x) for x in tag_snps)
    if tags is not None and not tags:
        raise ValueError("tag_snps is empty")
    c = len(edges) - 1
    counts = np.zeros(c, dtype=np.int64)
    pieces = {k: [] for k in (*_REQUIRED, "L1", "L2", "L3", "L4")}
    seen_snp, seen_block, found_tags = set(), set(), set()
    block_metadata = []
    for block in blocks:
        block_id = str(block.get("id", block.get("block", "")))
        if not block_id or any(x.isspace() for x in block_id) or block_id in seen_block:
            raise ValueError("Every LD block must have a nonempty, unique scalar id")
        seen_block.add(block_id)
        r, meta, nref, qc = _prepare_block(block, check_psd, min_maf, max_variants, log)
        duplicated = seen_snp.intersection(meta["snp"])
        if duplicated:
            raise ValueError(f"SNP identifiers occur in multiple blocks: {sorted(duplicated)[:3]}")
        seen_snp.update(meta["snp"])
        bin_af = meta["af"].copy()
        for edge in edges:
            bin_af[np.isclose(bin_af, edge, rtol=0, atol=1e-12)] = edge
        bins = np.minimum(np.searchsorted(edges, bin_af, side="right") - 1, c - 1)
        block_counts = np.bincount(bins, minlength=c)
        counts += block_counts
        use = np.ones(len(r), dtype=bool) if tags is None else np.array([x in tags for x in meta["snp"]])
        found_tags.update(meta["snp"][use])
        ld_origin = "sample_genotypes" if "genotype" in block else str(block.get("ld_origin", "supplied_ld"))
        if nref is not None:
            if not np.isfinite(float(nref)) or float(nref) < 1 or float(nref) != int(nref):
                raise ValueError("n_reference, when supplied, must be a positive integer")
        if finite_reference:
            if ld_origin in ("population", "generating_truth"):
                raise ValueError("Finite-reference correction must not be applied to population/generating LD")
            if nref is None or not np.isfinite(float(nref)) or float(nref) <= 2 or float(nref) != int(nref):
                raise ValueError("Finite-reference L2 correction requires integer n_reference > 2 in every block")
        tag_rows = np.flatnonzero(use)
        if len(tag_rows):
            # Bound score intermediates even when a reference block is large.
            indicator = np.eye(c, dtype=float)[bins]
            block_scores = {q: np.empty((len(tag_rows), c)) for q in range(1, 5)}
            for offset in range(0, len(tag_rows), 512):
                selected = tag_rows[offset:offset + 512]
                rr = r[selected, :]
                allowed = np.ones_like(rr, dtype=bool)
                if max_distance_bp is not None:
                    allowed = np.abs(meta["bp"][selected, None] - meta["bp"][None, :]) <= max_distance_bp
                for q in range(1, 5):
                    power = rr ** q
                    if q == 2 and finite_reference:
                        power = power - (1 - power) / (float(nref) - 2)
                        power[np.arange(len(selected)), selected] = 1.0
                    block_scores[q][offset:offset + len(selected)] = (power * allowed) @ indicator
            for q in range(1, 5):
                pieces[f"L{q}"].append(block_scores[q])
            for k in ("snp", "a1", "a2", "af", "chrom", "bp"):
                pieces[k].append(meta[k][use])
            pieces["block"].append(np.repeat(block_id, len(tag_rows)))
        block_metadata.append({"id": block_id, "chrom": str(meta["chrom"][0]),
                               "source_variants": int(len(r)), "tag_variants": int(use.sum()),
                               "source_counts": block_counts.tolist(),
                               "n_reference": None if nref is None else int(nref),
                               "ld_origin": ld_origin, "qc": qc})
        _log(log, f"LD block {block_id}: {len(r)} sources, {use.sum()} tags; strata {block_counts.tolist()}; QC {qc}")
    if not block_metadata:
        raise ValueError("No LD blocks supplied")
    if not pieces["snp"]:
        raise ValueError("No reference variants match the tag set")
    if np.any(counts <= 0):
        raise ValueError(f"Empty MAF source strata {np.flatnonzero(counts <= 0).tolist()}; merge bins or provide more reference variants")
    out = {k: np.concatenate(v, axis=0) for k, v in pieces.items()}
    out["counts"] = counts
    unknown_tags = 0 if tags is None else len(tags - found_tags)
    if unknown_tags:
        _warn(log, f"{unknown_tags} requested tag SNPs are absent or failed reference QC")
    warning = ("L3 and L4 use plug-in signed reference-LD powers; finite-reference bias is not corrected. "
               "Downstream uncertainty conditions on this reference. Cross-block LD is assumed zero.")
    _warn(log, warning)
    out["metadata"] = {"schema_version": SCHEMA_VERSION, "kind": "signed_ld_power_scores",
                       "effect_scale": "genotype_standardized", "allele_orientation": "minor_allele",
                       "maf_edges": edges.tolist(), "min_maf": float(min_maf),
                       "l2_finite_reference_correction": bool(finite_reference),
                       "l2_formula": "off-diagonal r^2-(1-r^2)/(n_reference-2); diagonal 1" if finite_reference else "r^2",
                       "higher_ld_power_bias_corrected": False, "ld_psd_checked": bool(check_psd),
                       "max_distance_bp": None if max_distance_bp is None else int(max_distance_bp), "source_variants": int(counts.sum()),
                       "tag_variants": int(len(out["snp"])), "requested_tags_missing": unknown_tags,
                       "blocks": block_metadata, "warnings": [warning],
                       "python_version": platform.python_version(), "numpy_version": np.__version__}
    _validate_scores(out)
    return out


def _validate_scores(scores):
    missing = [k for k in (*_REQUIRED, "L1", "L2", "L3", "L4", "counts", "metadata") if k not in scores]
    if missing:
        raise ValueError(f"Score archive is missing {missing}")
    n = len(scores["snp"])
    if not n:
        raise ValueError("Score archive has no tag SNPs")
    for k in ("snp", "a1", "a2", "chrom", "block"):
        scores[k] = _strings(scores[k], n, k)
    if len(set(scores["snp"])) != n:
        raise ValueError("Score SNP identifiers must be unique")
    if not np.all(np.isin(scores["a1"], list("ACGT"))) or not np.all(np.isin(scores["a2"], list("ACGT"))) or np.any(scores["a1"] == scores["a2"]):
        raise ValueError("Invalid score alleles")
    scores["af"] = np.asarray(scores["af"], dtype=float)
    if scores["af"].shape != (n,) or not np.isfinite(scores["af"]).all() or np.any((scores["af"] <= 0) | (scores["af"] > 0.5 + 1e-12)):
        raise ValueError("Score AF must be minor-allele frequency in (0,0.5]")
    bp = np.asarray(scores["bp"], dtype=float)
    if bp.shape != (n,) or not np.isfinite(bp).all() or np.any(bp <= 0) or np.any(bp != np.floor(bp)):
        raise ValueError("Score bp must be positive integers")
    scores["bp"] = bp.astype(np.int64)
    counts = np.asarray(scores["counts"], dtype=float)
    if counts.ndim != 1 or not len(counts) or not np.isfinite(counts).all() or np.any(counts <= 0) or np.any(counts != np.floor(counts)):
        raise ValueError("Score source counts must be positive integers")
    scores["counts"] = counts
    if counts.sum() < n:
        raise ValueError("Tag count cannot exceed the reference source count")
    for q in range(1, 5):
        a = np.asarray(scores[f"L{q}"], dtype=float)
        if a.shape != (n, len(counts)) or not np.isfinite(a).all():
            raise ValueError(f"L{q} must be a finite tags x strata matrix")
        scores[f"L{q}"] = a
    if np.any(scores["L4"] < -1e-12) or np.any(scores["L4"].sum(axis=1) < 1 - 1e-8):
        raise ValueError("L4 must be nonnegative and include the unit self-LD term")
    meta = scores["metadata"]
    if not isinstance(meta, dict) or meta.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unknown score schema_version; expected 1.0")
    if meta.get("effect_scale") != "genotype_standardized" or meta.get("allele_orientation") != "minor_allele":
        raise ValueError("Unsupported score scale or allele orientation")
    if meta.get("source_variants") != int(counts.sum()) or meta.get("tag_variants") != n:
        raise ValueError("Score counts disagree with metadata")
    if not meta.get("l2_finite_reference_correction", False):
        if np.any(scores["L2"] < -1e-12) or np.any(scores["L2"].sum(axis=1) < 1 - 1e-8):
            raise ValueError("Uncorrected L2 must be nonnegative and include unit self-LD")
    edges = np.asarray(meta.get("maf_edges", []), dtype=float)
    if len(edges) != len(counts) + 1 or edges[0] != 0 or edges[-1] != 0.5 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError("Score MAF edges disagree with source strata")


def save_scores(scores, path):
    """Atomically write a validated, pickle-free compressed NPZ score archive."""
    _validate_scores(scores)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {k: v for k, v in scores.items() if k != "metadata"}
    arrays["metadata"] = np.array(json.dumps(scores["metadata"], sort_keys=True, allow_nan=False))
    fd, temp = tempfile.mkstemp(prefix=".scores-", suffix=".npz", dir=path.parent)
    os.close(fd)
    try:
        np.savez_compressed(temp, **arrays)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)
    return path


def load_scores(path):
    """Read a validated score archive; metadata is returned as a dictionary."""
    with np.load(path, allow_pickle=False) as data:
        out = {k: data[k] for k in data.files}
    try:
        out["metadata"] = json.loads(str(out["metadata"].item()))
    except (KeyError, ValueError, AttributeError, json.JSONDecodeError) as exc:
        raise ValueError("Score metadata must be a JSON scalar string") from exc
    _validate_scores(out)
    return out


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _chr(x):
    x = str(x)
    return x[3:] if x.lower().startswith("chr") else x


def _read_bim(prefix):
    rows = []
    with open(str(prefix) + ".bim", encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) != 6:
                raise ValueError("BIM must contain six whitespace-separated columns")
            chrom, snp, cm, bp, a1, a2 = fields
            rows.append((chrom, snp, int(bp), a1, a2))
    if not rows:
        raise ValueError("BIM contains no variants")
    return rows


def plink_blocks(prefix, block_map, max_variants=10000):
    """Stream SNP-major PLINK BED into explicit genomic LD blocks.

    block_map TSV columns: chrom, start, end, block (1-based inclusive).
    Every BIM variant must lie in exactly one interval; input is not silently
    partitioned into arbitrary chunks. BED A1 dosage codes are 00=2, 01=missing,
    10=1, 11=0. FAM rows determine the reference sample count.
    """
    prefix, block_map = Path(prefix), Path(block_map)
    bim = _read_bim(prefix)
    with open(str(prefix) + ".fam", encoding="utf-8") as f:
        n_sample = sum(bool(line.strip()) for line in f)
    if n_sample < 3:
        raise ValueError("PLINK reference needs at least three FAM samples")
    with open(block_map, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not {"chrom", "start", "end", "block"}.issubset(reader.fieldnames or []):
            raise ValueError("Block map needs chrom,start,end,block tab-separated columns")
        intervals = list(reader)
    assignments = np.full(len(bim), -1, dtype=int)
    by_chrom = {}
    interval_ids = set()
    for b, row in enumerate(intervals):
        try:
            start, end = int(row["start"]), int(row["end"])
        except ValueError as exc:
            raise ValueError("Block map start/end must be integers") from exc
        if start <= 0 or end < start:
            raise ValueError("Block map intervals must be valid 1-based inclusive ranges")
        if not row["block"] or row["block"] in interval_ids:
            raise ValueError("Block map IDs must be nonempty and unique")
        interval_ids.add(row["block"])
        by_chrom.setdefault(_chr(row["chrom"]), []).append((start, end, b))
    # O(M log B) assignment rather than a Python scan of all variants per block.
    bim_chrom = np.asarray([_chr(row[0]) for row in bim])
    bim_bp = np.asarray([row[2] for row in bim], dtype=np.int64)
    for chrom, entries in by_chrom.items():
        entries.sort()
        starts = np.asarray([row[0] for row in entries], dtype=np.int64)
        ends = np.asarray([row[1] for row in entries], dtype=np.int64)
        original_ids = np.asarray([row[2] for row in entries], dtype=int)
        if np.any(starts[1:] <= ends[:-1]):
            raise ValueError("Block map intervals overlap")
        variants = np.flatnonzero(bim_chrom == chrom)
        positions = bim_bp[variants]
        where = np.searchsorted(starts, positions, side="right") - 1
        safe = np.maximum(where, 0)
        matched = (where >= 0) & (positions <= ends[safe])
        assignments[variants[matched]] = original_ids[where[matched]]
    if np.any(assignments < 0):
        raise ValueError(f"Block map leaves {np.sum(assignments < 0)} BIM variants unassigned")
    bytes_per_snp = (n_sample + 3) // 4
    expected = 3 + bytes_per_snp * len(bim)
    bed = Path(str(prefix) + ".bed")
    if bed.stat().st_size != expected:
        raise ValueError("BED byte size is inconsistent with FAM and BIM")
    lookup = np.array([2.0, np.nan, 1.0, 0.0])
    with open(bed, "rb") as stream:
        if stream.read(3) != b"\x6c\x1b\x01":
            raise ValueError("BED must use PLINK SNP-major mode (magic 6c1b01)")
        for b, row in enumerate(intervals):
            ix = np.flatnonzero(assignments == b)
            if not len(ix):
                continue
            if len(ix) > max_variants:
                raise ValueError(f"PLINK block {row['block']} exceeds max_variants={max_variants}")
            g = np.empty((n_sample, len(ix)))
            for j, source_index in enumerate(ix):
                stream.seek(3 + int(source_index) * bytes_per_snp)
                packed = np.frombuffer(stream.read(bytes_per_snp), dtype=np.uint8)
                code = ((packed[:, None] >> np.array([0, 2, 4, 6])) & 3).reshape(-1)[:n_sample]
                g[:, j] = lookup[code]
            rows = [bim[i] for i in ix]
            yield {"id": row["block"], "genotype": g, "snp": [x[1] for x in rows],
                   "chrom": [x[0] for x in rows], "bp": [x[2] for x in rows],
                   "a1": [x[3] for x in rows], "a2": [x[4] for x in rows]}


def manifest_blocks(manifest, max_variants=10000):
    """Read TSV block,path,format[,n_reference,block_map] relative to manifest.

    format=npz: archive has ld or genotype plus SNP metadata; block column sets id.
    format=bed: path is a PLINK prefix and block_map is its explicit interval TSV.
    """
    manifest = Path(manifest).resolve()
    with open(manifest, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"block", "path", "format"}.issubset(reader.fieldnames or []):
            raise ValueError("Reference manifest needs block,path,format tab-separated columns")
        for row in reader:
            path = Path(row["path"])
            if not path.is_absolute():
                path = manifest.parent / path
            if row["format"].lower() == "npz":
                with np.load(path, allow_pickle=False) as data:
                    block = {k: data[k] for k in data.files}
                block["id"] = row["block"]
                if row.get("n_reference"):
                    block["n_reference"] = int(row["n_reference"])
                yield block
            elif row["format"].lower() == "bed":
                if not row.get("block_map"):
                    raise ValueError("BED manifest rows require block_map")
                block_map = Path(row["block_map"])
                if not block_map.is_absolute():
                    block_map = manifest.parent / block_map
                yield from plink_blocks(path, block_map, max_variants=max_variants)
            else:
                raise ValueError("Reference format must be npz or bed")


def _manifest_provenance(manifest):
    manifest = Path(manifest).resolve()
    files = [manifest]
    with open(manifest, newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = Path(row["path"])
            path = path if path.is_absolute() else manifest.parent / path
            if row["format"].lower() == "bed":
                files.extend(Path(str(path) + suffix) for suffix in (".bed", ".bim", ".fam"))
                bm = Path(row["block_map"])
                files.append(bm if bm.is_absolute() else manifest.parent / bm)
            else:
                files.append(path)
    return [{"path": str(p), "sha256": _sha256(p), "bytes": p.stat().st_size} for p in dict.fromkeys(files)]


def add_score_parser(subparsers):
    """Register ``gcmsem score`` with an argparse subparser collection."""
    p = subparsers.add_parser("score", help="Build signed LD scores of orders 1–4")
    p.add_argument("--manifest", required=True, help="Reference block TSV (block,path,format)")
    p.add_argument("--out", required=True, help="Output score NPZ filename")
    p.add_argument("--maf-edges", default="0,0.5", help="Comma-separated stratum edges, from 0 to 0.5")
    p.add_argument("--min-maf", type=float, default=0.0)
    p.add_argument("--tags", help="One tag/HM3 SNP identifier per line; no source pruning")
    p.add_argument("--finite-reference", action="store_true", help="Apply off-diagonal finite-reference L2 correction")
    p.add_argument("--max-distance-bp", type=int, help="Optional additional within-block LD window")
    p.add_argument("--skip-psd-check", action="store_true", help="Skip expensive minimum LD eigenvalue checks")
    p.add_argument("--max-block-variants", type=int, default=10000)
    p.add_argument("--overwrite", action="store_true", help="Replace an existing score archive")
    p.set_defaults(func=score_cli)
    return p


def score_cli(args):
    if Path(args.out).exists() and not getattr(args, "overwrite", False):
        raise FileExistsError("Score archive already exists; choose another --out or use --overwrite")
    tags = None
    if args.tags:
        with open(args.tags, encoding="utf-8") as stream:
            tags = [line.strip() for line in stream if line.strip()]
    try:
        edges = tuple(float(x) for x in args.maf_edges.split(","))
    except ValueError as exc:
        raise ValueError("--maf-edges must be comma-separated numbers") from exc
    scores = build_scores(manifest_blocks(args.manifest, args.max_block_variants), maf_edges=edges,
                          tag_snps=tags, min_maf=args.min_maf, finite_reference=args.finite_reference,
                          max_distance_bp=args.max_distance_bp, check_psd=not args.skip_psd_check,
                          max_variants=args.max_block_variants)
    from . import __version__
    scores["metadata"]["gcmsem_version"] = __version__
    scores["metadata"]["command"] = getattr(args, "command", None)
    scores["metadata"]["input_files"] = _manifest_provenance(args.manifest)
    if args.tags:
        scores["metadata"]["tag_file"] = {"path": str(Path(args.tags).resolve()), "sha256": _sha256(args.tags)}
    save_scores(scores, args.out)
    _log(None, f"Wrote {args.out}: {scores['metadata']['tag_variants']} tags, {scores['metadata']['source_variants']} source variants, {len(scores['counts'])} strata")
    return 0
