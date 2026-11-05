"""Allele orientation and failure safety for the public Neale input tutorial."""
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
import pytest

from gcmsem.scores import build_scores, save_scores


SCRIPT = Path(__file__).resolve().parents[2] / "tutorials" / "prepare_neale_sumstats.py"
spec = importlib.util.spec_from_file_location("prepare_neale_sumstats", SCRIPT)
neale = importlib.util.module_from_spec(spec)
spec.loader.exec_module(neale)
HEADER = "variant\tminor_allele\tminor_AF\tlow_confidence_variant\tn_complete_samples\tbeta\tse\n"


def fixture_scores(tmp_path, *, ambiguous=False):
    # Two different SNP allele pairs at the same position must map separately.
    score = build_scores([dict(id="1:a", ld=np.eye(3), snp=["rsAG", "rsAC", "rsCT"],
                               a1=["A", "A", "C"], a2=["G", "G" if ambiguous else "C", "T"],
                               chrom=["1"] * 3, bp=[100, 100, 200], af=[.2, .3, .4])],
                         log=lambda _: None)
    path = tmp_path / "scores.npz"
    save_scores(score, path)
    return path


def raw_file(tmp_path, rows, *, header=HEADER):
    path = tmp_path / "raw.tsv.bgz"
    with gzip.open(path, "wt") as handle:
        handle.write(header + "\n".join(rows) + "\n")
    return path


def test_alt_beta_and_effect_frequency_survive_preparation(tmp_path):
    scores = fixture_scores(tmp_path)
    raw = raw_file(tmp_path, ["1:100:A:G\tA\t0.2\tfalse\t1000\t0.125\t0.02",
                              "chr01:100:A:C\tC\t0.3\tfalse\t900\t-0.25\t0.04",
                              "1:200:C:T\tT\t0.4\ttrue\t800\t0.5\t0.1",
                              "1:300:A:G\tG\t0.1\tfalse\t800\t0.1\t0.1"])
    out = tmp_path / "prepared.tsv.gz"
    metadata = neale.prepare(raw, scores, out)
    with gzip.open(out, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0] == dict(SNP="rsAG", A1="G", A2="A", BETA="0.125", SE="0.02",
                           N="1000", EAF="0.8", CHR="1", BP="100")
    assert rows[1]["SNP"] == "rsAC"
    assert float(rows[1]["EAF"]) == .3 and float(rows[1]["BETA"]) == -.25
    assert metadata["effect_scale"] == "allele"
    assert metadata["counts"]["retained"] == 2
    assert metadata["counts"]["eaf_complemented"] == 1
    assert metadata["counts"]["low_confidence"] == 1
    assert metadata["counts"]["not_score_tag"] == 1
    assert metadata["scores"]["source_variants"] == 3
    assert json.loads(Path(str(out) + ".json").read_text())["output"]["sha256"] == neale.sha256(out)
    assert "Retained SNPs: 2" in Path(str(out) + ".log").read_text()


@pytest.mark.parametrize("row,count", [
    ("1:100:A:G\tG\t0.2\tfalse\t100\tnan\t0.1", "invalid_numeric"),
    ("1:100:A:G\tG\t0.2\tfalse\t100\t0.2\t0", "invalid_numeric"),
    ("1:100:A:G\tG\t0.2\tfalse\t1.5\t0.2\t0.1", "invalid_numeric"),
    ("1:100:A:G\tC\t0.2\tfalse\t100\t0.2\t0.1", "invalid_minor_allele"),
    ("1:100:A:G\tG\t0.9\tfalse\t100\t0.2\t0.1", "invalid_minor_AF"),
    ("1:100:A:G\tG\tnan\tfalse\t100\t0.2\t0.1", "invalid_minor_AF"),
    ("1:100:A:G\tG\t0.2\tunknown\t100\t0.2\t0.1", "invalid_confidence_flag"),
    ("1:100:A:AG\tAG\t0.2\tfalse\t100\t0.2\t0.1", "unsupported_variant"),
])
def test_invalid_matched_rows_have_explicit_counts(row, count):
    source = io.StringIO(HEADER + row + "\n1:200:C:T\tT\t0.4\tfalse\t100\t0.2\t0.1\n")
    counts = neale.convert_rows(source, io.StringIO(),
                               {("1", 100, "AG"): "rsAG", ("1", 200, "CT"): "rsCT"}, lambda _: None)
    assert counts[count] == 1 and counts["retained"] == 1


def test_ambiguous_score_mapping_rejected(tmp_path):
    scores = fixture_scores(tmp_path, ambiguous=True)
    with pytest.raises(ValueError, match="Ambiguous score mapping"):
        neale.score_mapping(neale.load_scores(scores))


def test_existing_outputs_survive_bad_input_and_require_overwrite(tmp_path):
    scores = fixture_scores(tmp_path)
    raw = raw_file(tmp_path, ["1:100:A:G\tG\t0.2\tfalse\t100\t0.2\t0.1"])
    out = tmp_path / "prepared.tsv.gz"
    neale.prepare(raw, scores, out)
    paths = [out, Path(str(out) + ".json"), Path(str(out) + ".log")]
    original = [path.read_bytes() for path in paths]
    with pytest.raises(FileExistsError):
        neale.prepare(raw, scores, out)
    raw.write_bytes(b"not a gzip file")
    with pytest.raises(OSError):
        neale.prepare(raw, scores, out, overwrite=True)
    assert [path.read_bytes() for path in paths] == original
    assert not list(tmp_path.glob(".neale-*"))


@pytest.mark.parametrize("header,rows,message", [
    ("variant\tbeta\n", ["1:100:A:G\t0.1"], "Missing Neale columns"),
    (HEADER, ["1:100:A:G\tG"], "Malformed Neale input"),
    (HEADER, ["1:100:A:G\tG\t0.2\tfalse\t100\t0.2\t0.1"] * 2, "Duplicate matched SNP"),
])
def test_bad_schema_or_duplicate_never_publishes(tmp_path, header, rows, message):
    raw = raw_file(tmp_path, rows, header=header)
    scores = fixture_scores(tmp_path)
    out = tmp_path / "prepared.tsv.gz"
    with pytest.raises(ValueError, match=message):
        neale.prepare(raw, scores, out)
    assert not any(tmp_path.glob("prepared.tsv.gz*"))


def test_truncated_gzip_is_rejected(tmp_path):
    raw = raw_file(tmp_path, ["1:100:A:G\tG\t0.2\tfalse\t100\t0.2\t0.1"])
    raw.write_bytes(raw.read_bytes()[:-8])
    with pytest.raises(EOFError):
        neale.prepare(raw, fixture_scores(tmp_path), tmp_path / "prepared.tsv.gz")
    assert not any(tmp_path.glob("prepared.tsv.gz*"))
