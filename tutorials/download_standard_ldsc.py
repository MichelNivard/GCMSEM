#!/usr/bin/env python3
"""Download unchanged, checksum-pinned standard LDSC reference files.

The original Broad URLs in the LDSC wiki return 404 (checked 2026-11-05).
These DOI archives contain Oliver Pain's gzipped copies of the original files.
They are downloaded to the user's local reference directory, not redistributed
as part of the GCMSEM data deposit. No scores are recalculated by this script.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import logging
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import urllib.request


FILES = [
    dict(filename="eur_w_ld_chr.tar.gz", doi="10.5281/zenodo.8182036",
         url="https://zenodo.org/api/records/8182036/files/eur_w_ld_chr.tar.gz/content",
         original_url="https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2",
         bytes=33357890, md5="e2f16343c4cfaa76caa7d0c03d26b489",
         sha256="9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069"),
    dict(filename="w_hm3.snplist.gz", doi="10.5281/zenodo.7773502",
         url="https://zenodo.org/api/records/7773502/files/w_hm3.snplist.gz/content",
         original_url="https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2",
         bytes=5221397, md5="153ecc2bcfa740afafe656e6a384d769",
         sha256="3c876903ea996851fefabbb0460190ddbc8366fb7bc1ab4dea02335714d45a48"),
]
HM3_SHA256 = "ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def checked_download(item, destination):
    path = destination / item["filename"]
    if not path.exists():
        part = path.with_suffix(path.suffix + ".part")
        logging.info("DOWNLOAD %s", item["url"])
        with urllib.request.urlopen(item["url"], timeout=120) as src, part.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        if part.stat().st_size != item["bytes"] or sha256(part) != item["sha256"]:
            raise ValueError(f"Downloaded file failed checksum: {part}")
        part.replace(path)
    if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
        raise ValueError(f"Existing reference archive differs from the pinned file: {path}")
    logging.info("VERIFIED %s bytes=%d sha256=%s", path.name, item["bytes"], item["sha256"])
    return path


def extract_reference(archive, destination):
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        for member in members:
            name = PurePosixPath(member.name)
            if (name.is_absolute() or ".." in name.parts or
                    not name.parts or name.parts[0] != "eur_w_ld_chr" or
                    not (member.isfile() or member.isdir())):
                raise ValueError(f"Unsafe archive member: {member.name}")
        for member in members:
            path = destination.joinpath(*PurePosixPath(member.name).parts)
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            with handle.extractfile(member) as source:
                if path.exists():
                    expected = hashlib.sha256(source.read()).hexdigest()
                    if sha256(path) != expected:
                        raise ValueError(f"Existing extracted reference differs: {path}")
                else:
                    with path.open("wb") as target:
                        shutil.copyfileobj(source, target, length=1024 * 1024)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True, help="Local reference cache, outside the GCMSEM data deposit")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    archives = [checked_download(item, args.out) for item in FILES]
    extract_reference(archives[0], args.out)
    hm3 = args.out / "w_hm3.snplist"
    if not hm3.exists():
        with gzip.open(archives[1], "rb") as source, hm3.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    if sha256(hm3) != HM3_SHA256:
        raise ValueError("Extracted HM3 allele list differs from the pinned file")
    chromosomes = []
    for chrom in range(1, 23):
        prefix = args.out / "eur_w_ld_chr" / str(chrom)
        score = Path(str(prefix) + ".l2.ldscore.gz")
        count = Path(str(prefix) + ".l2.M_5_50")
        with gzip.open(score, "rt") as handle:
            header = next(handle).split()
            rows = sum(1 for _ in handle)
        m = int(count.read_text().strip())
        if header != ["CHR", "SNP", "BP", "CM", "MAF", "L2"] or m <= 0:
            raise ValueError(f"Unexpected standard reference format on chromosome {chrom}")
        chromosomes.append(dict(chromosome=chrom, rows=rows, M_5_50=m,
                                score_sha256=sha256(score), M_5_50_sha256=sha256(count)))
    if sum(c["rows"] for c in chromosomes) != 1290028 or sum(c["M_5_50"] for c in chromosomes) != 1173569:
        raise ValueError("Unexpected standard reference dimensions")
    provenance = {
        "files": FILES, "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_description": "Unchanged Alkes Group standard LDSC files, archived as gzip by Oliver Pain",
        "archive_license": "CC-BY-4.0; retain original Alkes Group and archive creator attribution",
        "documentation": "https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation",
        "both_ref_ld_chr_and_w_ld_chr": str(args.out / "eur_w_ld_chr") + "/",
        "M_5_50_total": 1173569, "score_rows_total": 1290028,
        "hm3_rows": 1217311, "hm3_sha256": HM3_SHA256,
        "chromosomes": chromosomes,
        "hla": "Use the published 6.* files; 6_old.* are included upstream but unused",
        "score_recalculation": False,
    }
    (args.out / "standard_reference.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    logging.info("SUCCESS standard reference ready: %s", args.out)


if __name__ == "__main__":
    main()
