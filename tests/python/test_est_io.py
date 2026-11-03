import csv
import json
from pathlib import Path

import numpy as np
import pytest

from gcmsem.cli import main
from gcmsem.io import _rows, harmonize_sumstats, read_manifest, write_results
from gcmsem.moments import estimate_moments
from gcmsem.scores import build_scores, save_scores


def table(path, header, rows):
    with Path(path).open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)


def make_fixture(tmp_path):
    blocks=[]
    for b in range(12):
        m=10
        blocks.append(dict(id=f'b{b}',ld=np.eye(m),snp=np.array([f'rs{b*m+i}' for i in range(m)]),
                           a1=np.repeat('A',m),a2=np.repeat('C',m),af=np.full(m,.2),
                           chrom=np.repeat('1',m),bp=np.arange(b*1000+1,b*1000+m+1)))
    scores=build_scores(blocks)
    rng=np.random.default_rng(991)
    y=rng.normal(size=(120,2)); y[:,1]+=.3*y[:,0]
    for c,name in enumerate(['A','B']):
        rows=[]
        for i,snp in enumerate(scores['snp']):
            flip=i%3==0
            rows.append([snp,'C' if flip else 'A','A' if flip else 'C',
                         (-1 if flip else 1)*y[i,c]/np.sqrt(.32),.1/np.sqrt(.32),.8 if flip else .2])
        table(tmp_path/f'{name}.tsv',['SNP','A1','A2','BETA','SE','EAF'],rows)
    table(tmp_path/'manifest.tsv',['trait','path','effect_scale'],[['A','A.tsv','allele'],['B','B.tsv','allele']])
    return scores,y


def test_harmonization_reorders_sign_and_scale(tmp_path):
    score,y=make_fixture(tmp_path)
    beta,se,aligned,traits,qc=harmonize_sumstats(tmp_path/'manifest.tsv',score)
    np.testing.assert_allclose(beta,y,atol=1e-14)
    np.testing.assert_allclose(se,.1)
    assert traits==['A','B']
    assert qc['trait_qc'][0]['sign_flipped']==40
    assert aligned['counts'].tolist()==score['counts'].tolist()
    assert all(len(f['sha256'])==64 for f in qc['input_files'])


def test_complete_case_qc_preserves_source_counts(tmp_path):
    score,y=make_fixture(tmp_path)
    file=tmp_path/'B.tsv'; lines=file.read_text().splitlines()
    row=lines[3].split('\t');row[4]='-1';lines[3]='\t'.join(row)
    file.write_text('\n'.join(lines)+'\n')
    beta,se,aligned,traits,qc=harmonize_sumstats(tmp_path/'manifest.tsv',score)
    assert len(beta)==119
    assert qc['trait_qc'][1]['invalid_numeric']==1
    assert aligned['counts'].sum()==120


def test_duplicate_snp_and_case_colliding_headers_fail(tmp_path):
    score,y=make_fixture(tmp_path)
    f=tmp_path/'A.tsv'; lines=f.read_text().splitlines();f.write_text('\n'.join([*lines,lines[1]])+'\n')
    with pytest.raises(ValueError,match='duplicate SNP'):
        harmonize_sumstats(tmp_path/'manifest.tsv',score)
    table(tmp_path/'bad.tsv',['SNP','BETA','beta'],[['rs1',1,2]])
    with pytest.raises(ValueError,match='duplicate column'):
        list(_rows(tmp_path/'bad.tsv'))


def test_palindrome_requires_unique_frequency_orientation(tmp_path):
    score,y=make_fixture(tmp_path)
    score['a2'][0]='T'
    for name in ['A','B']:
        f=tmp_path/f'{name}.tsv';lines=f.read_text().splitlines();row=lines[1].split('\t')
        # Original first row was swapped: C/A, -.8 per-allele direction.
        row[1]='T';lines[1]='\t'.join(row);f.write_text('\n'.join(lines)+'\n')
    assert len(harmonize_sumstats(tmp_path/'manifest.tsv',score)[0])==119
    beta,*_=harmonize_sumstats(tmp_path/'manifest.tsv',score,keep_palindromic=True)
    assert len(beta)==120
    np.testing.assert_allclose(beta,y)
    score['af'][0]=.5
    assert len(harmonize_sumstats(tmp_path/'manifest.tsv',score,keep_palindromic=True)[0])==119


def test_cli_real_file_roundtrip_logs_and_overwrite_guard(tmp_path):
    score,y=make_fixture(tmp_path)
    save_scores(score,tmp_path/'score.npz')
    table(tmp_path/'noise.tsv',['trait','A','B'],[['A',0,0],['B',0,0]])
    args=['estimate','--sumstats',str(tmp_path/'manifest.tsv'),'--scores',str(tmp_path/'score.npz'),
          '--out',str(tmp_path/'out'),'--blocks','12','--noise-intercept',str(tmp_path/'noise.tsv')]
    assert main(args)==0
    meta=json.loads((tmp_path/'out.metadata.json').read_text())
    assert meta['traits']==['A','B'] and meta['n_moments']==14
    assert meta['noise_mode']=='supplied_fixed'
    assert meta['score_sha256'] and meta['sumstats_manifest_sha256']
    log=(tmp_path/'out.log').read_text()
    for text in ['GCMSEM','COMMAND','QC','JACKKNIFE','OUTPUT','SUCCESS']:
        assert text in log
    rows=list(_rows(tmp_path/'out.moments.tsv'))
    assert len(rows)==14 and rows[-1][1]['key']=='B|B|B|B'
    assert main(args)==2
    assert main([*args,'--overwrite'])==0


def test_schema_requires_explicit_input_units(tmp_path):
    (tmp_path/'a.tsv').write_text('SNP\tA1\tA2\tBETA\tSE\nrs1\tA\tC\t1\t.1\n')
    table(tmp_path/'manifest.tsv',['trait','path'],[['A','a.tsv']])
    with pytest.raises(ValueError,match='effect_scale'):
        read_manifest(tmp_path/'manifest.tsv')


def test_conditional_bundle_covariance_and_completion_marker(tmp_path,monkeypatch):
    score,y=make_fixture(tmp_path)
    # Partition tags/targets into two independent disjoint source strata.
    annotation=np.tile([0,1],60)
    score['counts']=np.array([60,60])
    for q in range(1,5): score[f'L{q}']=np.eye(2)[annotation]
    result=estimate_moments(y,np.ones_like(y),score,noise_intercept=np.zeros((2,2)),n_blocks=12)
    paths=write_results(result,tmp_path/'x')
    assert len(paths)==18
    for stratum in (1,2):
        meta=json.loads((tmp_path/f'x.stratum-{stratum}.metadata.json').read_text())
        assert meta['target']=='conditional_annotation'
        rows=list(_rows(tmp_path/f'x.stratum-{stratum}.moments.tsv'))
        got=np.array([float(r['estimate']) for _,r in rows])
        np.testing.assert_allclose(got,result['conditional'][stratum-1])
    # An interrupted overwrite invalidates the old completed-run marker.
    import gcmsem.io as io
    monkeypatch.setattr(io,'_table',lambda *a,**k: (_ for _ in ()).throw(OSError('interrupted')))
    with pytest.raises(OSError,match='interrupted'):
        write_results(result,tmp_path/'x',overwrite=True)
    assert not (tmp_path/'x.metadata.json').exists()
