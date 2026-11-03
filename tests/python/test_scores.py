"""Numerical and ingestion contracts for signed LD-power score construction."""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse

from gcmsem.scores import (build_scores, load_scores, save_scores, plink_blocks,
                           manifest_blocks, score_cli, add_score_parser)


def block(r=None, *, ident="1:a", snp=None, af=None, **extra):
    if r is None:
        r = np.array([[1, -.4, .2], [-.4, 1, .1], [.2, .1, 1.]])
    n = r.shape[0]
    return dict(id=ident, ld=r, snp=snp or [f"{ident}:rs{i}" for i in range(n)],
                a1=["A"]*n, a2=["G"]*n, chrom=[ident.split(":")[0]]*n,
                bp=np.arange(1, n+1)*100, af=np.array(af if af is not None else [.1]*n), **extra)


def test_signed_power_identity_and_source_tag_distinction():
    b = block(af=[.05, .1, .4])
    s = build_scores([b], maf_edges=(0,.2,.5), tag_snps=[b["snp"][0]], log=lambda m: None)
    assert s["counts"].tolist() == [2,1]
    assert len(s["snp"]) == 1
    for q in range(1,5):
        assert s[f"L{q}"][0,0] == pytest.approx(1+(-.4)**q)
        assert s[f"L{q}"][0,1] == pytest.approx(.2**q)
    assert s["L3"][0,0] < 1


def test_allele_recode_invariance_and_sparse_input():
    b = block(af=[.1,.3,.2])
    base = build_scores([b], log=lambda m: None)
    recoded = dict(b)
    recoded["a1"] = ["G","A","G"]
    recoded["a2"] = ["A","G","A"]
    recoded["af"] = [0.9, .3, .8]
    signs = np.array([-1,1,-1])
    recoded["ld"] = sparse.csr_matrix(signs[:,None]*b["ld"]*signs[None,:])
    result = build_scores([recoded], log=lambda m: None)
    for key in ("L1","L2","L3","L4","af"):
        np.testing.assert_allclose(base[key], result[key], atol=1e-15)
    for key in ("a1", "a2"):
        np.testing.assert_equal(base[key], result[key])


def test_variable_blocks_and_finite_reference_diagonal():
    a = block(r=np.array([[1,.4],[.4,1]]), ident="1:a", n_reference=12)
    b = block(r=np.eye(4), ident="2:b", n_reference=12)
    s = build_scores([a,b], finite_reference=True, log=lambda m: None)
    np.testing.assert_allclose(s["L2"][:2,0], 1 + .16 - .84/10)
    # Independent population variants still have noisy sample-LD correction below zero.
    np.testing.assert_allclose(s["L2"][2:,0], 1-3/10)
    assert s["counts"].tolist() == [6]
    assert [x["source_variants"] for x in s["metadata"]["blocks"]] == [2,4]
    assert s["metadata"]["higher_ld_power_bias_corrected"] is False


def test_window_masks_off_diagonal_but_preserves_sources():
    b = block()
    s = build_scores([b], max_distance_bp=0, log=lambda m: None)
    for q in range(1,5):
        np.testing.assert_array_equal(s[f"L{q}"],np.ones((3,1)))
    assert s["counts"].tolist() == [3]


def test_genotype_missing_monomorphic_and_minor_orientation():
    g=np.array([[0,2,0,np.nan],[1,2,2,np.nan],[2,2,np.nan,np.nan],[0,2,2,np.nan]])
    b=block(r=np.eye(4)); b.pop("ld"); b.pop("af"); b["genotype"]=g
    s=build_scores([b],log=lambda m:None)
    assert s["counts"].tolist()==[2]
    np.testing.assert_allclose(s["af"],[3/8,1/3])
    assert s["a1"].tolist()==["A","G"]
    qc=s["metadata"]["blocks"][0]["qc"]
    assert qc["monomorphic"]==1 and qc["all_missing"]==1 and qc["missing_dosages"]==5
    centered=np.array([[0,2/3],[1,-4/3],[2,0],[0,-4/3]])
    centered[:,0]-=centered[:,0].mean()
    # Direct mean-imputed Pearson correlation uses observed dosage mean.
    raw=g[:,[0,2]].copy(); means=np.nanmean(raw,axis=0)
    raw=np.where(np.isfinite(raw),raw,means); raw[:,1]=2-raw[:,1]
    rr=np.corrcoef(raw,rowvar=False)
    np.testing.assert_allclose(s["L1"][:,0],rr.sum(axis=1))


@pytest.mark.parametrize("change,match",[
    ({"ld":np.array([[1,2,0],[2,1,0],[0,0,1]])},"correlation"),
    ({"ld":np.array([[1,.9,.9],[.9,1,-.9],[.9,-.9,1]])},"positive semidefinite"),
    ({"af":[.1,np.nan,.2]},"af must"),
    ({"snp":["same","same","other"]},"Duplicate"),
    ({"chrom":["1","1","2"]},"one chromosome"),
    ({"bp":[10,20.5,30]},"integer"),
    ({"a2":["A","G","G"]},"differ"),
])
def test_malformed_blocks_rejected(change,match):
    b=block(); b.update(change)
    with pytest.raises(ValueError,match=match): build_scores([b],log=lambda m:None)


def test_cross_block_duplicates_empty_bins_and_reference_size():
    with pytest.raises(ValueError,match="unique scalar"):
        build_scores([block(),block()],log=lambda m:None)
    a=block(); b=block(ident="1:b"); b["snp"]=a["snp"]
    with pytest.raises(ValueError,match="multiple blocks"):
        build_scores([a,b],log=lambda m:None)
    with pytest.raises(ValueError,match="Empty MAF"):
        build_scores([a],maf_edges=(0,.2,.5),log=lambda m:None)
    with pytest.raises(ValueError,match="n_reference"):
        build_scores([a],finite_reference=True,log=lambda m:None)


def test_archive_round_trip_and_tamper_detection(tmp_path):
    s=build_scores([block()],log=lambda m:None)
    path=save_scores(s,tmp_path/"scores.npz")
    loaded=load_scores(path)
    assert loaded["metadata"]==s["metadata"]
    for q in range(1,5): np.testing.assert_array_equal(loaded[f"L{q}"],s[f"L{q}"])
    with np.load(path,allow_pickle=False) as f: arrays={k:f[k] for k in f.files}
    arrays["counts"]=np.array([99]); np.savez(path,**arrays)
    with pytest.raises(ValueError,match="counts disagree"): load_scores(path)


def _write_bed(tmp_path):
    prefix=tmp_path/"ref"
    # Five samples checks trailing padded two-bit codes. Columns A1 dosages.
    g=np.array([[2,0,1],[1,0,2],[0,1,0],[np.nan,2,1],[2,2,2]],dtype=float)
    (tmp_path/"ref.fam").write_text("\n".join(f"F I{i} 0 0 0 -9" for i in range(len(g)))+"\n")
    (tmp_path/"ref.bim").write_text("1 rs1 0 100 A G\n1 rs2 0 200 C T\n2 rs3 0 500 A C\n")
    packed=bytearray(b"\x6c\x1b\x01")
    for col in g.T:
        code=[1 if np.isnan(x) else {2:0,1:2,0:3}[x] for x in col]
        code+= [0]*((-len(code))%4)
        packed.extend(sum(code[i+j]<<(2*j) for j in range(4)) for i in range(0,len(code),4))
    (tmp_path/"ref.bed").write_bytes(packed)
    bm=tmp_path/"blocks.tsv"; bm.write_text("chrom\tstart\tend\tblock\n1\t1\t300\t1:a\n2\t1\t1000\t2:b\n")
    return prefix,bm,g


def test_plink_bed_decoding_and_explicit_blocks(tmp_path):
    prefix,bm,g=_write_bed(tmp_path)
    blocks=list(plink_blocks(prefix,bm))
    assert [b["id"] for b in blocks]==["1:a","2:b"]
    np.testing.assert_equal(blocks[0]["genotype"],g[:,:2])
    np.testing.assert_equal(blocks[1]["genotype"],g[:,2:])
    score=build_scores(blocks,log=lambda m:None)
    assert score["counts"].tolist()==[3]
    assert score["metadata"]["blocks"][0]["n_reference"]==5
    bm.write_text("chrom\tstart\tend\tblock\n1\t1\t300\t1:a\n")
    with pytest.raises(ValueError,match="unassigned"): list(plink_blocks(prefix,bm))


def test_manifest_cli_and_provenance(tmp_path):
    b=block(); b.pop("id")
    np.savez(tmp_path/"block.npz",**b)
    manifest=tmp_path/"blocks.tsv"; manifest.write_text("block\tpath\tformat\n1:a\tblock.npz\tnpz\n")
    import argparse
    parser=argparse.ArgumentParser(); add_score_parser(parser.add_subparsers())
    args=parser.parse_args(["score","--manifest",str(manifest),"--out",str(tmp_path/"out.npz")])
    assert score_cli(args)==0
    s=load_scores(tmp_path/"out.npz")
    assert len(s["metadata"]["input_files"])==2
    assert all(len(x["sha256"])==64 for x in s["metadata"]["input_files"])
