#!/usr/bin/env python3
"""Fetch and verify untouched bulik/ldsc 1.0.1 for the standard LDSC tutorial.

The original LDSC commands run under a separate Python 2.7 environment. This
Python 3 setup helper does not translate or patch either LDSC command. Use
--fetch-only to obtain the source before creating its environment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

REPOSITORY = "https://github.com/bulik/ldsc.git"
COMMIT = "2fdeeb3b44379408794154993dbd6101b8946b7e"
BENCHMARK = {"numpy": "1.16.0", "scipy": "0.18.0", "pandas": "0.20.3",
             "bitarray": "0.8.3", "nose": "1.3.7", "python-dateutil": "2.8.2",
             "pytz": "2020.5", "six": "1.16.0"}


def git(path, *arguments):
    return subprocess.check_output(["git", "-C", str(path), *arguments], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("external/standard-ldsc"))
    parser.add_argument("--python2", default="python2.7", help="Python executable in the separate stock LDSC environment")
    parser.add_argument("--fetch-only", action="store_true")
    args = parser.parse_args()
    destination = args.out.resolve()
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="stock-ldsc-", dir=destination.parent) as temporary:
            checkout = Path(temporary) / "checkout"
            subprocess.run(["git", "clone", "--no-checkout", "--depth", "1", REPOSITORY, str(checkout)], check=True)
            subprocess.run(["git", "-C", str(checkout), "fetch", "--depth", "1", "origin", COMMIT], check=True)
            subprocess.run(["git", "-C", str(checkout), "checkout", "--detach", COMMIT], check=True)
            os.replace(checkout, destination)
    if not (destination / ".git").is_dir() or git(destination, "rev-parse", "HEAD") != COMMIT:
        parser.error(f"Existing --out must be the pinned original LDSC checkout {COMMIT}; choose a new path")
    if git(destination, "status", "--porcelain", "--untracked-files=no"):
        parser.error("Stock LDSC must have no tracked modifications; do not patch the commands or regression code")
    sources = [destination / name for name in ("ldsc.py", "munge_sumstats.py")]
    sources += sorted((destination / "ldscore").glob("*.py"))
    hashes = {str(path.relative_to(destination)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    provenance = {"repository": REPOSITORY, "commit": COMMIT, "version": "1.0.1",
                  "tracked_source_modified": False, "source_sha256": hashes,
                  "reference_procedure": "https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation",
                  "benchmark_dependencies": BENCHMARK}
    if not args.fetch_only:
        check = (
            "import sys,json,platform,pkg_resources; "
            "assert sys.version_info[:2]==(2,7), 'Original LDSC requires Python 2.7'; "
            "import numpy,scipy,pandas,bitarray; "
            "from scipy.stats import chi2; "
            "assert chi2.isf(0.05,1)>3; "
            "print(json.dumps({'python':platform.python_version(),'dependencies':"
            "{name:pkg_resources.get_distribution(name).version for name in " + repr(list(BENCHMARK)) + "}}))"
        )
        runtime = json.loads(subprocess.check_output([args.python2, "-c", check], text=True))
        provenance.update(runtime)
        for script in ("ldsc.py", "munge_sumstats.py"):
            subprocess.run([args.python2, str(destination / script), "--help"], check=True,
                           stdout=subprocess.DEVNULL)
        provenance["cli_help_checks"] = ["ldsc.py passed", "munge_sumstats.py passed"]
        if runtime["dependencies"] != BENCHMARK:
            print("NOTE installed dependency versions differ from the executed benchmark; inspect provenance.")
    target = destination.parent / (destination.name + "-setup.json")
    target.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"READY untouched bulik/ldsc: {destination}")
    print(f"PROVENANCE {target}")
    if args.fetch_only:
        print("Source fetched only. Create the separate Python 2.7 environment, then rerun with --python2.")


if __name__ == "__main__":
    main()
