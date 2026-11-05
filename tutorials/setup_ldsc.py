#!/usr/bin/env python3
"""Fetch the pinned official Python 3 LDSC used by the empirical benchmark.

Install NumPy, pandas, SciPy and bitarray first; this script never changes the
Python environment. It applies only the documented pandas CLI syntax fix and
verifies that the regression, IRWLS and jackknife code is unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

REPOSITORY = "https://github.com/CBIIT/ldsc.git"
COMMIT = "1f09cf0c5e61a38a281016fefcc1e0583cc78402"
PATCH_PATH = "ldscore/sumstats.py"
OLD = ".drop('SNP', 1)"
NEW = ".drop('SNP', axis=1)"
STATISTICS = ("ldscore/regressions.py", "ldscore/irwls.py", "ldscore/jackknife.py")
DEPENDENCIES = {"numpy": "2.3.5", "pandas": "2.2.3", "scipy": "1.17.0", "bitarray": "3.12.3"}


def git(path, *arguments):
    return subprocess.check_output(["git", "-C", str(path), *arguments], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("external/ldsc"))
    args = parser.parse_args()
    destination = args.out.resolve()
    versions = {}
    for package in DEPENDENCIES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            parser.error("Missing dependency. In your tutorial environment run: python -m pip install " +
                         " ".join(f"{name}=={version}" for name, version in DEPENDENCIES.items()))
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="ldsc-", dir=destination.parent) as temporary:
            checkout = Path(temporary) / "checkout"
            subprocess.run(["git", "clone", "--no-checkout", "--depth", "1", REPOSITORY, str(checkout)], check=True)
            subprocess.run(["git", "-C", str(checkout), "fetch", "--depth", "1", "origin", COMMIT], check=True)
            subprocess.run(["git", "-C", str(checkout), "checkout", "--detach", COMMIT], check=True)
            os.replace(checkout, destination)
    if not (destination / ".git").is_dir() or git(destination, "rev-parse", "HEAD") != COMMIT:
        parser.error(f"Existing --out must be the pinned LDSC checkout {COMMIT}; choose an empty new path")
    changed = set(git(destination, "diff", "--name-only", "HEAD").splitlines())
    if changed - {PATCH_PATH}:
        parser.error(f"Existing checkout has unrelated modifications: {sorted(changed - {PATCH_PATH})}")
    original = subprocess.check_output(["git", "-C", str(destination), "show", f"HEAD:{PATCH_PATH}"], text=True)
    if original.count(OLD) != 1:
        raise ValueError("Pinned pandas compatibility expression is not unique")
    patched = original.replace(OLD, NEW)
    source = destination / PATCH_PATH
    if source.read_text() not in (original, patched):
        raise ValueError("Refusing to replace additional local edits in ldscore/sumstats.py")
    source.write_text(patched)
    hashes = {}
    for relative in STATISTICS:
        current = (destination / relative).read_bytes()
        expected = subprocess.check_output(["git", "-C", str(destination), "show", f"HEAD:{relative}"])
        if current != expected:
            raise ValueError(f"Statistical code differs from the pinned commit: {relative}")
        hashes[relative] = hashlib.sha256(current).hexdigest()
    subprocess.run([sys.executable, str(destination / "ldsc.py"), "--help"], check=True, stdout=subprocess.DEVNULL)
    provenance = {"repository": REPOSITORY, "commit": COMMIT,
                  "package_version": "3.0.2", "cli_banner_version": "3.0.1",
                  "compatibility_patch": {"file": PATCH_PATH, "old": OLD, "new": NEW},
                  "statistical_source_sha256": hashes, "python": platform.python_version(),
                  "dependencies": versions, "benchmark_dependency_versions": DEPENDENCIES,
                  "cli_help_check": "passed"}
    target = destination.parent / (destination.name + "-setup.json")
    target.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"READY pinned official LDSC: {destination}")
    print(f"PROVENANCE {target}")
    if versions != DEPENDENCIES:
        print("NOTE dependency versions differ from the executed tutorial benchmark; inspect provenance.")


if __name__ == "__main__":
    main()
