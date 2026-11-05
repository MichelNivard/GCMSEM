#!/usr/bin/env python3
"""Known-LD check using independent Gaussian source mixtures and GWAS noise.

Vary block lengths and LD strengths. All variants are sources and tags, so no
representative-tag approximation is needed. This checks the moment mapping;
it is not an interval-coverage study or a realistic population simulation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import savemat
from scipy.sparse import csc_matrix
from gcmsem.io import write_results
from gcmsem.moments import estimate_moments


def simulate(out, seed=20261105, n_blocks=1200):
    out = Path(out);out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    lengths = rng.integers(20, 121, n_blocks)
    M = int(lengths.sum()); N = 200000; target_H = .4
    probabilities = [.02, .2]; names = ['Sparse', 'Polygenic']
    beta, se, powers, blocks, truth_v, truth_D = [], [], [], [], [], []
    for i, length in enumerate(lengths):
        rho = rng.uniform(.15,.9) * rng.choice([-1,1])
        R = rho**np.abs(np.arange(length)[:,None]-np.arange(length)[None,:])
        v = np.column_stack([(rng.random(length)<pi)*target_H/(M*pi) for pi in probabilities])
        direct = rng.normal(size=v.shape)*np.sqrt(v)
        observed = R@direct + np.linalg.cholesky(R)@rng.normal(size=v.shape)/np.sqrt(N)
        beta.append(observed);se.append(np.full_like(observed,1/np.sqrt(N)))
        powers.append(np.column_stack([(R**q).sum(1) for q in range(1,5)]))
        blocks.extend([f'block{i:04d}']*length)
        truth_v.append(v.sum(0));truth_D.append((v*((R**2)@v)).sum(0))
    beta,se,P=np.vstack(beta),np.vstack(se),np.vstack(powers)
    counts=np.array([M]);bins=np.zeros(M,dtype=int)
    sc={f'L{q}':P[:,q-1:q] for q in range(1,5)}
    sc.update(counts=counts,block=np.array(blocks),chrom=np.repeat('1',M))
    # Use a divisor of M so the unmodified upstream drops no tail. In GCMSEM,
    # original generating LD blocks remain whole, even if deletion sizes vary.
    B=max(b for b in range(30,201) if M % b == 0) if any(M%b==0 for b in range(30,201)) else 100
    # The fixed seed gives M=84,111 and B=159. No simulation source is dropped.
    if M % B: raise ValueError('Choose a simulation seed whose M has a divisor in 30:200')
    savemat(out/'sld4m_input.mat',dict(beta=beta,se=se,ell2=sc['L2'],ell4=sc['L4'],
        annot=csc_matrix(np.ones((M,1))),weights=(1/np.maximum(P[:,1],1))[:,None],
        idxldsc=np.arange(1,M+1,dtype=float)[:,None],no_blocks=B,counts=counts[:,None]))
    np.savez(out/'design.npz',L2=sc['L2'],bins=bins,counts=counts,chrom=sc['chrom'],block=sc['block'])
    for i,name in enumerate(names):
        r=estimate_moments(beta[:,i:i+1],se[:,i:i+1],sc,traits=[name],n_blocks=B)
        write_results(r,out/name)
        np.savez(out/(name+'.estimates.npz'),conditional=r['conditional'],conditional_jackknife=r['conditional_jackknife'])
    H=np.sum(truth_v,axis=0);D=np.sum(truth_D,axis=0)
    a=dict(traits=names,n_blocks=B,n_tags=M,source_counts=[M],n_samples=N,seed=seed,
        generating_ld_blocks=n_blocks,block_length_min=int(lengths.min()),block_length_max=int(lengths.max()),
        probabilities=probabilities,true_effective_counts=(H**2/D).tolist(),
        true_genetic_variances=H.tolist(),ld_reference='exact generating LD; all source SNPs are tags')
    (out/'analysis.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',default='tutorial-results/polygenicity-simulation')
    a=p.parse_args();simulate(a.out)
