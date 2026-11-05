#!/usr/bin/env python3
"""Download/check original public GWAS and create raw chi-square caches."""
import argparse
import hashlib
import json
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd


def verify_or_download(row, raw_dir, download):
    dest = raw_dir/row['filename']
    if not dest.exists():
        if not download: raise FileNotFoundError(f'{dest}: use --download')
        print('Downloading', row['code'], flush=True)
        temporary = dest.with_suffix(dest.suffix+'.partial')
        with urllib.request.urlopen(row['url'], timeout=90) as response, temporary.open('wb') as f:
            for block in iter(lambda: response.read(4*1024*1024), b''): f.write(block)
        temporary.rename(dest)
    h = hashlib.sha256()
    with dest.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''): h.update(block)
    if h.hexdigest() != row['sha256']: raise ValueError(f'SHA-256 mismatch: {dest}')
    print('Verified', row['code'], flush=True)
    return dest


def order_update(h, values): h.update(('\n'.join(values)+'\n').encode())


def mapping(path, out):
    h = hashlib.sha256(); parts = []
    for chunk in pd.read_csv(path, sep='\t', compression='gzip', usecols=['variant','rsid'],
            dtype=str, keep_default_na=False, chunksize=400000):
        order_update(h,chunk.variant)
        parts.append(pd.to_numeric(chunk.rsid.str.extract(r'^rs([0-9]+)$',expand=False),
            errors='coerce').fillna(0).to_numpy(dtype=np.int64))
    ids = np.concatenate(parts); np.save(out/'neale_row_rsids.npy',ids)
    meta = dict(n_rows=len(ids),variant_order_sha256=h.hexdigest(),missing_rsid=int((ids==0).sum()))
    (out/'neale_row_order.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Variant map:',len(ids),'rows',flush=True)


def neale(task):
    code, path, out = task; h = hashlib.sha256(); stats = []; qc = []
    for chunk in pd.read_csv(path, sep='\t', compression='gzip',
            usecols=['variant','tstat','low_confidence_variant'],
            dtype={'variant':str,'tstat':float,'low_confidence_variant':str},chunksize=400000):
        order_update(h,chunk.variant); z = chunk.tstat.to_numpy(); stats.append(z*z)
        qc.append((chunk.low_confidence_variant=='false').to_numpy() & np.isfinite(z))
    chisq = np.concatenate(stats); valid = np.concatenate(qc)
    meta = json.loads((out/'neale_row_order.json').read_text())
    if len(chisq) != meta['n_rows'] or h.hexdigest() != meta['variant_order_sha256']:
        raise ValueError(f'{code}: GWAS/annotation row order differs; use an explicit variant join')
    np.savez_compressed(out/(code+'.raw.npz'),chisq=chisq,valid=valid)
    print(code, 'parsed',len(chisq),'eligible',int(valid.sum()),flush=True)


def gscan(path, out):
    parts = []; ids = []
    for chunk in pd.read_csv(path,sep='\t',compression='gzip',usecols=['RSID','STAT'],
            dtype={'RSID':str,'STAT':float},chunksize=400000):
        ids.append(pd.to_numeric(chunk.RSID.str.extract(r'^rs([0-9]+)$',expand=False),
            errors='coerce').fillna(0).to_numpy(dtype=np.int64))
        parts.append(chunk.STAT.to_numpy())
    ids = np.concatenate(ids); chisq = np.concatenate(parts)
    valid = np.isfinite(chisq) & (chisq>=0) & (ids>0)
    np.savez_compressed(out/'GSCAN_DrinksPerWeek.raw.npz',rsids=ids,chisq=chisq,valid=valid)
    print('GSCAN parsed',len(ids),'eligible',int(valid.sum()),flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('sld4m-source-gwas.tsv'))
    p.add_argument('--raw-dir',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--download',action='store_true')
    p.add_argument('--workers',type=int,default=2)
    a = p.parse_args(); a.raw_dir.mkdir(parents=True,exist_ok=True); a.out.mkdir(parents=True,exist_ok=True)
    rows = pd.read_csv(a.manifest,sep='\t',dtype=str).to_dict('records')
    paths = {r['code']:verify_or_download(r,a.raw_dir,a.download) for r in rows}
    mapping(paths.pop('variants'),a.out)
    gscan_path = paths.pop('GSCAN_DrinksPerWeek')
    with ProcessPoolExecutor(a.workers) as pool:
        list(pool.map(neale,[(code,path,a.out) for code,path in paths.items()]))
    gscan(gscan_path,a.out)
    (a.out/'source-manifest.tsv').write_bytes(a.manifest.read_bytes())


if __name__ == '__main__': main()
