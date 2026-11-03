#!/usr/bin/env python3
"""Estimate all 34 cumulants in the actual GWASBrewer Monte Carlo experiment.

The target is the generating population cumulant, not the empirical cumulant
of the particular set of sampled joint effects. Both comparisons are exported.
"""
from __future__ import annotations
import argparse
import itertools
import json
import logging
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.io import mmread

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"python"))
from gcmsem.scores import build_scores, save_scores, load_scores


def reference_blocks(work):
    info=pd.read_csv(work/"ld/blocks.tsv",sep="\t")
    af=pd.read_csv(work/"ld/af.tsv",sep="\t").af.to_numpy()
    offset=0
    for rep in range(2):
        for row in info.itertuples():
            n=int(row.size); ix=np.arange(offset,offset+n); offset+=n
            yield dict(id=f"rep{rep+1}_b{row.block}",
                       ld=mmread(work/f"ld/block_{row.block:02d}.mtx"),
                       snp=np.array([f"sim{i+1}" for i in ix]),a1=np.repeat("A",n),
                       a2=np.repeat("C",n),af=af[int(row.start)-1:int(row.end)],
                       chrom=np.repeat("19",n),bp=ix*1000+1)


def target_cumulants(T, source_variance, pi, indices):
    source={1:np.zeros_like(source_variance),2:source_variance,
            3:np.sqrt(2)*source_variance**1.5/np.sqrt(pi),
            4:source_variance**2*(6/pi-3)}
    return np.array([np.sum(np.prod(T[list(t)],axis=0)*source[len(t)]) for t in indices])


def empirical_cumulants(joint,indices):
    center=joint-joint.mean(axis=0)
    covariance=center.T@center/len(center)
    vals=[]
    for t in indices:
        if len(t)==1: val=joint[:,t[0]].mean()
        else: val=np.prod(center[:,t],axis=1).mean()
        if len(t)==4:
            a,b,c,d=t
            val-=covariance[a,b]*covariance[c,d]+covariance[a,c]*covariance[b,d]+covariance[a,d]*covariance[b,c]
        vals.append(val)
    return np.array(vals)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--work",type=Path,default=Path("validation/work"))
    parser.add_argument("--results",type=Path,default=Path("validation/results"))
    parser.add_argument("--limit",type=int)
    parser.add_argument("--resume",action="store_true",help="reuse completed replication rows from this exact experiment")
    args=parser.parse_args()
    work=args.work;out=args.results;out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    logger=logging.getLogger("validation")
    from gcmsem.moments import estimate_moments
    from gcmsem.io import write_results
    score_file=work/"scores.npz"
    if score_file.exists(): scores=load_scores(score_file)
    else:
        scores=build_scores(reference_blocks(work),maf_edges=(0.,.5),finite_reference=False,
                            check_psd=False,log=logger)
        save_scores(scores,score_file)
    # The literal scenario name "null" is not a missing value.
    manifest=pd.read_csv(work/"manifest.tsv",sep="\t",keep_default_na=False)
    if args.limit: manifest=manifest.iloc[:args.limit]
    detail=[]; summary=[]
    if args.resume and (out/"moment_replicates.tsv").exists() and (out/"run_status.tsv").exists():
        detail=pd.read_csv(out/"moment_replicates.tsv",sep="\t",keep_default_na=False).to_dict("records")
        summary=pd.read_csv(out/"run_status.tsv",sep="\t",keep_default_na=False).to_dict("records")
    completed={x["run"] for x in summary if x["status"]=="ok"}
    for rec in manifest.itertuples():
        if rec.name in completed and (work/"fits"/(rec.name+".metadata.json")).exists():
            logger.info("Reusing completed replication %s",rec.name)
            continue
        t0=time.monotonic(); data=pd.read_csv(work/(rec.name+".tsv.gz"),sep="\t")
        sx=np.sqrt(2*data.AF.to_numpy()*(1-data.AF.to_numpy()))
        beta=data[["beta_A","beta_B","beta_C"]].to_numpy()*sx[:,None]
        se=data[["se_A","se_B","se_C"]].to_numpy()*sx[:,None]
        truth_joint=data[["joint_A","joint_B","joint_C"]].to_numpy()
        T=np.loadtxt(work/(rec.name+".T.tsv"),delimiter="\t",ndmin=2)
        source_variance=pd.read_csv(work/(rec.name+".variance.tsv"),sep="\t").dir_h2.to_numpy()/rec.J
        try:
            res=estimate_moments(beta,se,scores,traits=["A","B","C"],max_order=4,n_blocks=int(rec.blocks),logger=logger)
            indices=res["indices"]
            expected=target_cumulants(T,source_variance,.15,indices)
            sampled=empirical_cumulants(truth_joint,indices)
            err=np.sqrt(np.diag(res["sampling_cov"]))
            for j,t in enumerate(indices):
                key="|".join("ABC"[i] for i in t)
                detail.append(dict(run=rec.name,scenario=rec.scenario,replicate=rec.replicate,key=key,order=len(t),
                    target=expected[j],sample_truth=sampled[j],estimate=res["estimates"][j],se=err[j],
                    covered=bool(abs(res["estimates"][j]-expected[j])<=1.96*err[j])))
            # Uniform, deterministic scaling makes SEM optimizers well conditioned.
            # Since all observed expected variances are .2/J, path coefficients
            # are unchanged. Sampling covariance is transformed by the same factors.
            scale=np.array([(rec.J/.2)**(len(t)/2) for t in indices])
            scaled=dict(res)
            scaled["estimates"]=res["estimates"]*scale
            scaled["sampling_cov"]=res["sampling_cov"]*scale[:,None]*scale[None,:]
            scaled["jackknife"]=res["jackknife"]*scale[None,:]
            scaled["conditional"]=res["conditional"]*scale[None,:]
            scaled["conditional_jackknife"]=res["conditional_jackknife"]*scale[None,None,:]
            scaled["metadata"]=dict(res["metadata"],validation_uniform_effect_multiplier=float(np.sqrt(rec.J/.2)))
            prefix=work/"fits"/rec.name;prefix.parent.mkdir(exist_ok=True)
            write_results(scaled,prefix,overwrite=True)
            noise=np.loadtxt(work/(rec.name+".noise.tsv"),ndmin=2)
            summary.append(dict(run=rec.name,scenario=rec.scenario,status="ok",seconds=time.monotonic()-t0,
                                n_moments=len(indices),n_blocks=len(res["jackknife"]),
                                max_noise_intercept_error=float(np.max(abs(res["noise_intercept"]-noise)))))
            np.savetxt(work/"fits"/(rec.name+".truth.tsv"),expected*scale)
        except Exception as e:
            logger.exception("Failed %s",rec.name)
            summary.append(dict(run=rec.name,scenario=rec.scenario,status="failed",error=str(e),seconds=time.monotonic()-t0))
        pd.DataFrame(detail).to_csv(out/"moment_replicates.tsv",sep="\t",index=False)
        pd.DataFrame(summary).to_csv(out/"run_status.tsv",sep="\t",index=False)
    details=pd.DataFrame(detail)
    aggregate=[]
    grouped=[] if details.empty else details.groupby(["scenario","key"],sort=False)
    for (scenario,key),d in grouped:
        aggregate.append(dict(scenario=scenario,key=key,order=int(d.order.iloc[0]),replicates=len(d),
            target=float(d.target.iloc[0]),mean_estimate=d.estimate.mean(),bias=(d.estimate-d.target).mean(),
            empirical_sd=d.estimate.std(),mean_se=d.se.mean(),coverage=d.covered.mean(),
            rms_z=np.sqrt(np.mean(((d.estimate-d.target)/d.se)**2))))
    pd.DataFrame(aggregate).to_csv(out/"moment_summary.tsv",sep="\t",index=False)
    (out/"simulation_metadata.json").write_text(json.dumps(dict(
        software="GWASBrewer 0.3.0.0259",seed_base=20261103,
        reference="GWASBrewer bundled HapMap3 chromosome 19; minor-allele recoding",
        ld_patterns=39,independent_blocks=78,min_block=140,max_block=863,
        source_distribution="15% nonzero: centered standardized gamma(shape=2); 85% zero",
        genotype_scale="standardized",sampling_noise="GWASBrewer correlated normal errors",
        sample_size=int(manifest.N.iloc[0]),sample_overlap_fraction=.5,
        interval="normal 1.96 * full-refit delete-block jackknife SE",
        caveat="Small Monte Carlo study; does not establish general interval calibration. LD is known generating LD; reference-estimation uncertainty is excluded.",
        runs=summary),indent=2)+"\n")


if __name__=="__main__": main()
