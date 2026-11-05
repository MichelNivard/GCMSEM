#!/usr/bin/env python3
"""Fetch pinned tutorial assets from Zenodo; verify bytes before extraction.

The checked-in manifest pins content, while --record can select its permanent
Zenodo record after deposition. A record number is never inferred or invented.
No credentials are needed to download a published open record.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile
import urllib.request
import zipfile


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def valid_asset(asset):
    name = asset["filename"]
    if not name or Path(name).name != name or name in (".", ".."):
        raise ValueError("Asset filename must be a plain basename")
    checksum = asset["sha256"]
    if len(checksum) != 64 or any(c not in "0123456789abcdef" for c in checksum):
        raise ValueError("Each asset requires an exact SHA-256 checksum")
    if not isinstance(asset["size_bytes"], int) or asset["size_bytes"] <= 0:
        raise ValueError("Each asset requires its positive byte size")


def download(asset, url, directory):
    valid_asset(asset)
    if not url.startswith("https://"):
        raise ValueError("Downloads must use HTTPS")
    destination = directory / asset["filename"]
    if destination.exists():
        if destination.stat().st_size == asset["size_bytes"] and digest(destination) == asset["sha256"]:
            print(f"VERIFIED cached {destination}", flush=True)
            return destination
        raise ValueError(f"Existing file failed its checksum: {destination}; move it aside before retrying")
    fd, temporary = tempfile.mkstemp(prefix="download-", suffix=".part", dir=directory)
    try:
        with os.fdopen(fd, "wb") as output, urllib.request.urlopen(url, timeout=120) as response:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        if Path(temporary).stat().st_size != asset["size_bytes"] or digest(temporary) != asset["sha256"]:
            raise ValueError(f"Downloaded bytes failed verification: {asset['filename']}")
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    print(f"VERIFIED downloaded {destination}", flush=True)
    return destination


def extract(archive, destination):
    """Allow ordinary files only, never symlinks or path traversal."""
    with zipfile.ZipFile(archive) as source:
        seen = set()
        for item in source.infolist():
            path = PurePosixPath(item.filename)
            mode = item.external_attr >> 16
            if (path.is_absolute() or ".." in path.parts or "\\" in item.filename
                    or not path.parts or ":" in path.parts[0]
                    or (mode & 0o170000) not in (0, 0o040000, 0o100000)):
                raise ValueError(f"Unsafe archive entry: {item.filename}")
            if item.filename in seen:
                raise ValueError(f"Duplicate archive entry: {item.filename}")
            seen.add(item.filename)
            target = (destination / item.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise ValueError(f"Archive path escapes destination: {item.filename}")
            if target.exists() and not item.is_dir():
                with source.open(item) as handle:
                    checksum = hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(handle.read()).hexdigest()
                if not target.is_file() or digest(target) != checksum:
                    raise ValueError(f"Refusing to replace modified file: {target}")
        source.extractall(destination)
    print(f"EXTRACTED {archive.name} to {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("data_manifest.json"))
    parser.add_argument("--record", help="Published Zenodo record number, overriding the manifest location")
    parser.add_argument("--out", type=Path, default=Path("tutorial-data"))
    parser.add_argument("--asset", action="append", help="Asset ID to download; repeat as needed; default all")
    parser.add_argument("--no-extract", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema_version") != "1.0":
        raise ValueError("Unsupported tutorial manifest version")
    record = args.record or manifest.get("zenodo_record_id")
    if not record or not str(record).isdigit():
        raise ValueError("The data release has not yet been published. Supply --record with its published Zenodo record number; no download URL is fabricated.")
    request = urllib.request.Request(f"https://zenodo.org/api/records/{record}", headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        record_data = json.load(response)
    files = {item["key"]: item for item in record_data["files"]}
    assets = manifest["assets"]
    if args.asset:
        missing = set(args.asset) - {a["id"] for a in assets}
        if missing:
            raise ValueError(f"Unknown asset IDs: {sorted(missing)}")
        assets = [a for a in assets if a["id"] in args.asset]
    args.out.mkdir(parents=True, exist_ok=True)
    cache = args.out / ".downloads"
    cache.mkdir(exist_ok=True)
    for asset in assets:
        valid_asset(asset)
        if asset["filename"] not in files:
            raise ValueError(f"Record {record} does not contain {asset['filename']}")
        entry = files[asset["filename"]]
        if entry["size"] != asset["size_bytes"]:
            raise ValueError("Record file size differs from the pinned release")
        url = entry["links"].get("self") or entry["links"].get("content")
        if not url:
            raise ValueError("Record has no download link")
        archive = download(asset, url, cache)
        if not args.no_extract and asset.get("format") == "zip":
            extract(archive, args.out)
    print("SUCCESS all requested assets match the pinned release", flush=True)


if __name__ == "__main__":
    main()
