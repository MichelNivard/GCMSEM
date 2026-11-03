"""Strict GWAS harmonization and versioned, labeled SEM interchange files."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import tempfile

import numpy as np

from .moments import jackknife_cov, validate_traits


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def _rows(path):
    with _open_text(path) as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"empty table: {path}")
        delimiter = "\t" if "\t" in first else "," if "," in first else None
        if delimiter:
            names = [x.strip().lower() for x in next(csv.reader([first], delimiter=delimiter))]
            reader = csv.DictReader(handle, fieldnames=names, delimiter=delimiter)
            if len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError(f"duplicate column names: {path}")
            for number, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f"wrong number of fields in {path}, row {number}")
                yield number, {str(k).strip().lower(): str(v).strip() for k, v in row.items()}
        else:
            columns = [x.lower() for x in first.split()]
            if len(columns) != len(set(columns)):
                raise ValueError(f"duplicate column names: {path}")
            for number, line in enumerate(handle, 2):
                if not line.strip():
                    continue
                values = line.split()
                if len(values) != len(columns):
                    raise ValueError(f"wrong number of fields in {path}, row {number}")
                yield number, dict(zip((c.lower() for c in columns), values))


def read_manifest(path):
    """Read trait/path/effect_scale manifest; paths relative to manifest directory."""
    path = Path(path).resolve()
    entries = []
    for number, row in _rows(path):
        required = {"trait", "path", "effect_scale"}
        if not required <= set(row):
            raise ValueError("sumstats manifest requires trait, path, effect_scale columns")
        if row["effect_scale"] not in ("allele", "genotype_standardized"):
            raise ValueError(f"manifest row {number}: effect_scale must be allele or genotype_standardized")
        f = Path(row["path"])
        if not f.is_absolute():
            f = path.parent / f
        if not f.is_file():
            raise ValueError(f"summary-statistics file does not exist: {f}")
        entries.append(dict(row, path=str(f.resolve())))
    if not entries:
        raise ValueError("sumstats manifest is empty")
    validate_traits([row["trait"] for row in entries], len(entries))
    return entries


_ALIASES = {
    "snp": ("snp", "rsid", "rs_id", "variant_id"),
    "a1": ("a1", "effect_allele", "ea", "alt"),
    "a2": ("a2", "other_allele", "nea", "ref"),
    "beta": ("beta", "b", "effect"),
    "se": ("se", "standard_error", "stderr"),
    "eaf": ("eaf", "effect_allele_frequency", "a1freq", "a1_freq"),
}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def harmonize_sumstats(manifest, scores, *, min_info=None, min_n=None,
                       af_tolerance=0.1, allow_strand_flips=False,
                       keep_palindromic=False, logger=None):
    """Align trait tables to score tags and retain a common complete-case SNP set.

    Required GWAS fields SNP,A1,A2,BETA,SE accept documented aliases. EAF is
    optional for nonpalindromic SNPs and refers to the reported effect allele.
    A/T and C/G SNPs are dropped by default. Keeping them requires EAF and a
    unique frequency-based orientation; midpoint ambiguity is always dropped.
    Duplicated IDs cause an error, including duplicates outside the tag set.
    No p-value filter is performed. Source annotation counts are never reduced
    when tags are removed. Inputs are assumed on the same genome build.
    """
    log = logger or logging.getLogger("gcmsem")
    entries = read_manifest(manifest) if isinstance(manifest, (str, Path)) else list(manifest)
    traits = validate_traits([r["trait"] for r in entries], len(entries))
    if not 0 <= af_tolerance < 0.5:
        raise ValueError("af_tolerance must be in [0,0.5)")
    reference = np.asarray(scores["snp"]).astype(str)
    if len(set(reference)) != len(reference):
        raise ValueError("score SNP IDs are duplicated")
    lookup = {s: i for i, s in enumerate(reference)}
    a1, a2 = np.asarray(scores["a1"]).astype(str), np.asarray(scores["a2"]).astype(str)
    af = np.asarray(scores["af"], dtype=float)
    if np.any((af <= 0) | (af >= 1)) or not np.isfinite(af).all():
        raise ValueError("score allele frequencies must lie strictly between zero and one")
    n, p = len(reference), len(entries)
    beta, se = np.full((n, p), np.nan), np.full((n, p), np.nan)
    stats, input_files = [], []
    for col, entry in enumerate(entries):
        scale = entry.get("effect_scale")
        if scale not in ("allele", "genotype_standardized"):
            raise ValueError("every trait manifest entry requires an explicit effect_scale")
        path = Path(entry["path"])
        qc = {"trait": traits[col], "rows": 0, "not_in_reference": 0, "invalid_numeric": 0,
              "invalid_alleles": 0, "allele_mismatch": 0, "palindromic": 0,
              "frequency_mismatch": 0, "low_info": 0, "low_n": 0, "retained": 0,
              "sign_flipped": 0, "strand_flipped": 0}
        seen, columns = set(), None
        for number, row in _rows(path):
            qc["rows"] += 1
            if columns is None:
                columns = {}
                for canonical, aliases in _ALIASES.items():
                    found = [a for a in aliases if a in row]
                    if len(found) > 1:
                        raise ValueError(f"{path}: ambiguous {canonical} columns {found}; keep one")
                    if not found and canonical != "eaf":
                        raise ValueError(f"{path}: missing required {canonical} column")
                    columns[canonical] = found[0] if found else None
                if min_info is not None and "info" not in row:
                    raise ValueError(f"{path}: min_info requested but INFO column missing")
                if min_n is not None and "n" not in row:
                    raise ValueError(f"{path}: min_n requested but N column missing")
            snp = row[columns["snp"]]
            if not snp:
                raise ValueError(f"{path}, row {number}: empty SNP ID")
            if snp in seen:
                raise ValueError(f"{path}, row {number}: duplicate SNP ID {snp}")
            seen.add(snp)
            if snp not in lookup:
                qc["not_in_reference"] += 1
                continue
            i = lookup[snp]
            try:
                b, s = float(row[columns["beta"]]), float(row[columns["se"]])
                freq = None if columns["eaf"] is None else float(row[columns["eaf"]])
                if not np.isfinite(b) or not np.isfinite(s) or s <= 0 \
                        or (freq is not None and (not np.isfinite(freq) or not 0 <= freq <= 1)):
                    raise ValueError()
                if min_info is not None:
                    info = float(row["info"])
                    if not np.isfinite(info) or info < min_info:
                        qc["low_info"] += 1
                        continue
                if min_n is not None:
                    sample_n = float(row["n"])
                    if not np.isfinite(sample_n) or sample_n < min_n:
                        qc["low_n"] += 1
                        continue
            except (ValueError, TypeError):
                qc["invalid_numeric"] += 1
                continue
            ea, oa = row[columns["a1"]].upper(), row[columns["a2"]].upper()
            if ea not in "ACGT" or oa not in "ACGT" or len(ea) != 1 or len(oa) != 1 or ea == oa:
                qc["invalid_alleles"] += 1
                continue
            palindromic = {ea, oa} in ({"A", "T"}, {"C", "G"})
            candidates = []
            if (ea, oa) == (a1[i], a2[i]):
                candidates.append((1, False))
            if (ea, oa) == (a2[i], a1[i]):
                candidates.append((-1, False))
            if allow_strand_flips or palindromic:
                ca, co = ea.translate(_COMPLEMENT), oa.translate(_COMPLEMENT)
                if (ca, co) == (a1[i], a2[i]):
                    candidates.append((1, True))
                if (ca, co) == (a2[i], a1[i]):
                    candidates.append((-1, True))
            if palindromic:
                if not keep_palindromic or freq is None:
                    qc["palindromic"] += 1
                    continue
                candidates = [(sign, strand) for sign, strand in candidates
                              if abs((freq if sign == 1 else 1 - freq) - af[i]) <= af_tolerance]
                # Both orientations can match near AF=0.5; never resolve by luck.
                if len({sign for sign, _ in candidates}) != 1:
                    qc["palindromic"] += 1
                    continue
            if not candidates:
                qc["allele_mismatch"] += 1
                continue
            sign, strand = candidates[0]
            if freq is not None and abs((freq if sign == 1 else 1 - freq) - af[i]) > af_tolerance:
                qc["frequency_mismatch"] += 1
                continue
            factor = np.sqrt(2 * af[i] * (1 - af[i])) if scale == "allele" else 1.0
            beta[i, col], se[i, col] = sign * b * factor, s * factor
            qc["retained"] += 1
            qc["sign_flipped"] += int(sign == -1)
            qc["strand_flipped"] += int(strand)
        if columns is None:
            raise ValueError(f"summary-statistics table has no data rows: {path}")
        log.info("QC %s", json.dumps(qc, sort_keys=True))
        stats.append(qc)
        input_files.append({"trait": traits[col], "path": str(path.resolve()),
                            "sha256": sha256_file(path), "input_effect_scale": scale,
                            "phenotype_unit": entry.get("phenotype_unit", "unspecified")})
    keep = np.isfinite(beta).all(axis=1) & np.isfinite(se).all(axis=1)
    if keep.sum() < 3:
        raise ValueError("fewer than three common valid tag SNPs after harmonization")
    aligned = dict(scores)
    for key in ("snp", "a1", "a2", "af", "chrom", "bp", "block", "L1", "L2", "L3", "L4"):
        aligned[key] = np.asarray(scores[key])[keep]
    log.info("HARMONIZATION common_tags=%d original_tags=%d source_variants=%g", keep.sum(), n, np.sum(scores["counts"]))
    return beta[keep], se[keep], aligned, traits, {
        "trait_qc": stats, "input_files": input_files, "n_tags_before_qc": n,
        "n_tags_after_qc": int(keep.sum()), "af_tolerance": af_tolerance,
        "allow_strand_flips": allow_strand_flips, "keep_palindromic": keep_palindromic,
        "source_counts_unchanged": True, "common_tag_complete_case": True,
        "genome_build_check": "User must ensure shared genome build; matching IDs alone does not verify build"}


def _atomic_write(path, write):
    path = Path(path)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            write(handle)
        os.replace(name, path)
    except BaseException:
        try:
            os.unlink(name)
        except FileNotFoundError:
            pass
        raise


def _table(path, header, rows):
    def write(handle):
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    _atomic_write(path, write)


def _json_clean(value):
    if isinstance(value, dict):
        return {str(k): _json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_clean(v) for v in value]
    if isinstance(value, np.generic):
        return _json_clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_results(result, prefix, *, overwrite=False, conditional=True):
    """Write labeled pooled/conditional moments, full V, deletions and metadata.

    All names carry explicit trait tuples. Metadata is installed last and is the
    completion marker. Existing results require overwrite=True. Annotation
    bundles use suffix .stratum-N and can be read directly by the R package.
    """
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    traits, indices = result["traits"], result["indices"]
    keys = ["|".join(traits[i] for i in t) for t in indices]
    values, V = np.asarray(result["estimates"]), np.asarray(result["sampling_cov"])
    if values.shape != (len(keys),) or V.shape != (len(keys), len(keys)):
        raise ValueError("result estimates/covariance do not match tuple keys")
    suffixes = (".moments.tsv", ".sampling_cov.tsv", ".jackknife.tsv", ".covariance.tsv",
                ".correlation.tsv", ".metadata.json")
    paths = [Path(str(prefix) + suffix) for suffix in suffixes]
    if not overwrite and any(path.exists() for path in paths):
        raise FileExistsError(f"output exists for {prefix}; use overwrite=True or a new prefix")
    # Invalidate the completion marker before replacing any bundle member.
    # An interrupted overwrite must never masquerade as a complete older run.
    if overwrite:
        paths[5].unlink(missing_ok=True)
    sd = np.sqrt(np.maximum(np.diag(V), 0))
    _table(paths[0], ["key", "order", "estimate", "se"],
           ((key, len(t), value, s) for key, t, value, s in zip(keys, indices, values, sd)))
    _table(paths[1], ["key", *keys], ((key, *row) for key, row in zip(keys, V)))
    _table(paths[2], ["block", *keys], ((i, *row) for i, row in enumerate(result["jackknife"])))
    lookup = {t: value for t, value in zip(indices, values)}
    covariance = np.array([[lookup[tuple(sorted((a, b)))] for b in range(len(traits))]
                           for a in range(len(traits))])
    variance = np.diag(covariance)
    denom = np.sqrt(np.maximum(variance, 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = covariance / np.outer(denom, denom)
    correlation[(variance <= 0), :] = np.nan
    correlation[:, (variance <= 0)] = np.nan
    _table(paths[3], ["trait", *traits], ((name, *row) for name, row in zip(traits, covariance)))
    _table(paths[4], ["trait", *traits], ((name, *row) for name, row in zip(traits, correlation)))
    metadata = dict(result["metadata"], output_prefix=str(prefix.resolve()),
                    moment_keys=keys, output_files=[path.name for path in paths])
    extra = []
    conditional_values = np.asarray(result.get("conditional", []))
    if conditional and conditional_values.ndim == 2 and len(conditional_values) > 1:
        fractions = np.asarray(result["metadata"]["jackknife_fractions"])
        for c, point in enumerate(conditional_values):
            jk = result["conditional_jackknife"][:, c, :]
            sub_covariance = jackknife_cov(jk, point, fractions)
            active = np.diag(sub_covariance) > 0
            rank = 0
            if active.any():
                pseudo = point + (1 / fractions[:, None] - 1) * (point - jk)
                factor = (pseudo - fractions @ pseudo)[:, active] \
                         * np.sqrt(fractions / (1 - fractions) / len(jk))[:, None] \
                         / np.sqrt(np.diag(sub_covariance)[active])
                rank = int(np.linalg.matrix_rank(factor, tol=1e-4))
            sub_lookup = dict(zip(indices, point))
            sub_genetic = np.array([[sub_lookup[tuple(sorted((a,b))) ] for b in range(len(traits))]
                                    for a in range(len(traits))])
            sub_metadata = dict(result["metadata"], target="conditional_annotation", annotation_index=c,
                                sampling_covariance_rank=rank,
                                genetic_covariance_min_eigenvalue=float(np.linalg.eigvalsh(sub_genetic)[0]))
            sub = dict(result, estimates=point, sampling_cov=sub_covariance,
                       jackknife=jk, metadata=sub_metadata)
            extra += write_results(sub, str(prefix) + f".stratum-{c + 1}", overwrite=overwrite, conditional=False)
        metadata["conditional_prefixes"] = [prefix.name + f".stratum-{c + 1}" for c in range(len(conditional_values))]
    _atomic_write(paths[5], lambda handle: json.dump(_json_clean(metadata), handle, indent=2, allow_nan=False))
    return [str(p) for p in paths] + extra
