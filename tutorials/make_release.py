#!/usr/bin/env python3
"""Build deterministic, data-only tutorial assets for external archival hosting.

This release builder uses explicit top-level selections, excludes hidden/temp
files, and rejects raw reference formats. It does not upload or publish files.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

CAUSAL = ("50_irnt", "21002_irnt", "21001_raw", "4080_raw")
PANEL = ("GSCAN_DrinksPerWeek", "30890_raw", "4080_irnt", "1200", "20544_11",
         "30780_raw", "3063_irnt", "6150_1", "6150_2", "50_irnt", "21002_irnt",
         "21001_irnt", "4079_irnt", "3062_irnt", "48_irnt")
FORBIDDEN = (".bed", ".bim", ".fam", ".bgen", ".vcf", ".vcf.gz", ".pgen", ".pvar", ".psam")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def files_under(root, relative):
    target = root / relative
    if not target.exists():
        raise FileNotFoundError(target)
    return [p for p in ([target] if target.is_file() else target.rglob("*"))
            if p.is_file() and not any(part.startswith(".") for part in p.relative_to(root).parts)
            and not p.name.endswith((".tmp", ".part"))]


def make_archive(root, out, name, selected):
    selected = sorted(set(selected))
    inventory = []
    for path in selected:
        if path.is_symlink() or path.name.lower().endswith(FORBIDDEN):
            raise ValueError(f"Raw reference or link cannot enter release: {path}")
        relative = path.relative_to(root).as_posix()
        if not path.stat().st_size:
            raise ValueError(f"Empty release file: {relative}")
        inventory.append({"path": relative, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    target = out / name
    temporary = target.with_suffix(".zip.part")
    with zipfile.ZipFile(temporary, "w", allowZip64=True) as archive:
        for path, item in zip(selected, inventory):
            info = zipfile.ZipInfo(item["path"], date_time=(2026, 11, 5, 0, 0, 0))
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_STORED if path.suffix in (".gz", ".npz", ".zip") else zipfile.ZIP_DEFLATED
            with path.open("rb") as source, archive.open(info, "w", force_zip64=True) as destination:
                h = hashlib.sha256()
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    destination.write(block)
                    h.update(block)
                if h.hexdigest() != item["sha256"]:
                    raise ValueError(f"Input changed during packaging: {path}")
        info = zipfile.ZipInfo(f"RELEASE_FILES_{name[:-4]}.json", date_time=(2026, 11, 5, 0, 0, 0))
        info.external_attr = 0o100644 << 16
        archive.writestr(info, json.dumps(inventory, indent=2) + "\n")
    with zipfile.ZipFile(temporary) as check:
        failure = check.testzip()
        if failure:
            raise ValueError(f"Archive CRC failed: {failure}")
    temporary.replace(target)
    return {"filename": name, "format": "zip", "size_bytes": target.stat().st_size,
            "sha256": sha256(target), "file_count": len(inventory)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("data_manifest.json"))
    parser.add_argument("--record", help="Actual published Zenodo record ID, if already known")
    parser.add_argument("--include-results", action="store_true")
    args = parser.parse_args()
    root = args.data_dir.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    common = files_under(root, "provenance") + files_under(root, "RIGHTS.md")
    selections = {
        "scores": files_under(root, "scores") + common,
        "causal": [root / "sumstats" / f"{code}.tsv.gz" for code in CAUSAL]
                  + [root / "manifests" / f"{name}.tsv" for name in ("height_weight", "bmi_sbp")] + common,
        "panel": [root / "sumstats" / f"{code}.tsv.gz" for code in PANEL]
                 + [root / "panel.npz", root / "manifests/panel.tsv", root / "manifests/bmi_sbp_irnt.tsv"] + common,
    }
    if args.include_results:
        # Only the completed primary analyses enter this release. Development
        # attempts and custom-reference LDSC sensitivities remain separate.
        gcmsem = [p for p in (root / "results/correlations").glob("gcmsem*") if p.is_file()]
        if not gcmsem or not (root / "results/standard_ldsc_verified/status.json").is_file():
            raise ValueError("Complete GCMSEM and verified standard LDSC outputs are required")
        standard = [p for p in files_under(root, "results/standard_ldsc_verified")
                    if "work" not in p.relative_to(root / "results/standard_ldsc_verified").parts]
        selections["results"] = files_under(root, "results/causal") + gcmsem + standard + common
    assets = []
    for key, selection in selections.items():
        print(f"PACKAGING {key} files={len(set(selection))}", flush=True)
        asset = make_archive(root, args.out, f"gcmsem-{key}-v1.zip", selection)
        assets.append({"id": key, **asset})
    manifest = {"schema_version": "1.0", "release": "empirical-v1",
                "zenodo_record_id": args.record,
                "publication_status": "published" if args.record else "prepared; awaiting Zenodo account authorization",
                "raw_reference_included": False, "assets": assets}
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    (args.out / "SHA256SUMS").write_text("".join(f"{a['sha256']}  {a['filename']}\n" for a in assets))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
