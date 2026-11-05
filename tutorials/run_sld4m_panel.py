#!/usr/bin/env python3
"""Run the unchanged upstream S-LD4M function on each full-GWAS input."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
from scipy.io import loadmat
from prepare_original_sld4m import TRAITS, sha256

PIN = 'd8ca72212d1618abfa0a9695d27fa49f45f851cd9641674809fe6a3ed4b5917a'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--upstream',type=Path,required=True)
    p.add_argument('--octave',default='octave')
    p.add_argument('--jobs',type=int,default=1)
    p.add_argument('--traits',nargs='+',choices=list(TRAITS.values()))
    a = p.parse_args()
    if sha256(a.upstream/'SLD4M.m') != PIN: raise ValueError('Use the pinned, unmodified SLD4M.m')
    def quote(x): return "'"+str(Path(x).resolve()).replace("'","''")+"'"
    def run(trait):
        output = a.reference/(trait+'.result.mat')
        pending = a.reference/(trait+'.result.pending.mat')
        pending.unlink(missing_ok=True)
        command = 'addpath('+quote(Path(__file__).parent)+'); run_original_sld4m('+','.join(map(quote,
            [a.reference,a.reference/(trait+'.input.mat'),a.upstream,pending]))+');'
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        log_path = a.reference/(trait+'.log')
        print('Starting',trait,flush=True)
        with log_path.open('w') as log:
            subprocess.run([a.octave,'--quiet','--eval',command],stdout=log,stderr=subprocess.STDOUT,
                check=True,env=env)
        result = loadmat(pending,variable_names=['ma_jk'])
        if result['ma_jk'].dtype.kind != 'f' or result['ma_jk'].shape != (75,100):
            raise RuntimeError(f'{trait}: invalid native result storage')
        pending.replace(output)
        print(trait,log_path.read_text().splitlines()[-1],flush=True)
    with ThreadPoolExecutor(a.jobs) as pool: list(pool.map(run,a.traits or TRAITS.values()))


if __name__ == '__main__': main()
