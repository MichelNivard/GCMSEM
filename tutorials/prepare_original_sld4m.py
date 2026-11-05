#!/usr/bin/env python3
"""Repack the author's full reference losslessly and match raw GWAS by rsID."""
import argparse
import hashlib
import io
import json
import struct
import zipfile
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

TRAITS = dict(GSCAN_DrinksPerWeek='DrinksPerWeek', **{
    '30890_raw': 'VitaminD', '4080_irnt': 'SBP', '1200': 'Insomnia',
    '20544_11': 'Depression', '30780_raw': 'LDL', '3063_irnt': 'FEV1',
    '6150_1': 'HeartAttack', '6150_2': 'Angina', '50_irnt': 'Height',
    '21002_irnt': 'Weight', '21001_irnt': 'BMI', '4079_irnt': 'DBP',
    '3062_irnt': 'FVC', '48_irnt': 'Waist'})


class StoredZipMember(io.RawIOBase):
    """Seek within an uncompressed ZIP member without making a second copy."""
    def __init__(self, path, info):
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError('Expected the original ZIP_STORED reference archive')
        self.handle = open(path, 'rb')
        self.handle.seek(info.header_offset)
        header = struct.unpack('<4s5H3I2H', self.handle.read(30))
        self.base = info.header_offset + 30 + header[-2] + header[-1]
        self.length, self.position = info.file_size, 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else offset + (self.position if whence == 1 else self.length)
        return self.position

    def read(self, size=-1):
        size = min(size if size >= 0 else self.length-self.position, self.length-self.position)
        self.handle.seek(self.base+self.position)
        data = self.handle.read(size); self.position += len(data)
        return data

    def readinto(self, buffer):
        data = self.read(len(buffer)); buffer[:len(data)] = data; return len(data)

    def close(self):
        self.handle.close(); super().close()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()


def convert(archive, out):
    import h5py
    with zipfile.ZipFile(archive) as z:
        info = z.getinfo('baselineLD.1kg.l2l4.mat')
        example = z.read('example.m').decode()
        (out/'upstream-example.m').write_text(example)
        with z.open('height_example.mat') as h:
            height = loadmat(h, variable_names=['chisq','sumstat_RSIDs_as_ints'])
    with StoredZipMember(archive, info) as raw, h5py.File(raw, 'r') as f:
        ncols, nrows = f['annot'].shape
        metadata = dict(nrows=nrows, ncols=ncols,
            RSIDs=f['RSIDs'][0,:].astype(np.int64), EW2=f['EW2'][0,:], EW4=f['EW4'][0,:])
        headers = [''.join(chr(x) for x in f[r][:].ravel()) for r in f['header'][:,0]]
        for key in ('ell2', 'ell4', 'annot'):
            dest = out/(key+'.f64'); sums = np.zeros(ncols)
            print('Converting', key, (nrows,ncols), flush=True)
            with dest.open('wb') as target:
                for start in range(0,nrows,100000):
                    block = np.ascontiguousarray(f[key][:,start:start+100000].T, dtype='<f8')
                    block.tofile(target); sums += block.sum(axis=0)
            if dest.stat().st_size != nrows*ncols*8: raise RuntimeError('Incomplete conversion')
            metadata[key+'_sums'] = sums
        savemat(out/'reference_metadata.mat', metadata, do_compression=True, oned_as='column')
        provenance = dict(reference_member=info.filename, reference_rows=nrows, annotations=headers,
            archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
            storage='Lossless original float64 values, original row/column order',
            no_blocks=100, upstream_example=example,
            reference_download='https://www.dropbox.com/sh/iiyftw01gdpt6un/AACU7AmWK45RxTmDJvRkdKhIa?dl=1')
        (out/'reference-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    match(metadata['RSIDs'], height['sumstat_RSIDs_as_ints'].ravel(), height['chisq'].ravel(),
          np.ones(height['chisq'].size,dtype=bool), out, 'author_height_example')


def match(reference_ids, ids, chisq, valid, out, trait):
    eligible = np.flatnonzero(valid & (ids>0) & np.isfinite(chisq) & (chisq>=0))
    # MATLAB intersect(...,'stable'): unique IDs, first occurrence, reference order.
    _, ir, ig = np.intersect1d(reference_ids, ids[eligible], return_indices=True)
    order = np.argsort(ir); ir = ir[order]; ig = eligible[ig[order]]
    savemat(out/(trait+'.input.mat'), dict(idx_ref=ir+1, chisq=chisq[ig]), do_compression=True, oned_as='column')
    meta = dict(trait=trait, input_rows=len(ids), eligible_rows=len(eligible), matched_snps=len(ir),
        blocks=100, native_rows_used=(len(ir)//100)*100, native_tail_rows_unused=len(ir)%100,
        selection='finite nonnegative chi-square, publisher low-confidence filter (Neale), own reference rsID intersection',
        chisq_sha256=hashlib.sha256(np.asarray(chisq[ig],dtype='<f8').tobytes()).hexdigest())
    (out/(trait+'.input.json')).write_text(json.dumps(meta,indent=2)+'\n')
    print(trait, len(ir), 'matched SNPs', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--raw-cache', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--skip-conversion', action='store_true')
    a = p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    if not a.skip_conversion: convert(a.archive,a.out)
    reference_ids = loadmat(a.out/'reference_metadata.mat',variable_names=['RSIDs'])['RSIDs'].ravel()
    neale_ids = np.load(a.raw_cache/'neale_row_rsids.npy')
    for code, trait in TRAITS.items():
        with np.load(a.raw_cache/(code+'.raw.npz')) as data:
            match(reference_ids, data['rsids'] if 'rsids' in data else neale_ids,
                  data['chisq'], data['valid'], a.out, trait)


if __name__ == '__main__': main()
