#!/usr/bin/env python3
"""Summarize executed validation, retaining failed and withheld-inference runs."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
r=root/"validation/results"
mom=pd.read_csv(r/"moment_replicates.tsv",sep="\t",keep_default_na=False)
run=pd.read_csv(r/"run_status.tsv",sep="\t",keep_default_na=False)
sem=pd.read_csv(r/"sem_status.tsv",sep="\t",keep_default_na=False)
par=pd.read_csv(r/"sem_parameters.tsv",sep="\t",na_values=["NA"],keep_default_na=False)
exact=pd.read_csv(r/"exact_model_recovery.tsv",sep="\t",keep_default_na=False)
toggles=pd.read_csv(r/"moment_toggles.tsv",sep="\t",keep_default_na=False)
meta=json.loads((r/"simulation_metadata.json").read_text())

def md(headers,rows):
    return "| "+" | ".join(headers)+" |\n|"+"|".join(["---"]*len(headers))+"|\n"+"\n".join("| "+" | ".join(map(str,row))+" |" for row in rows)
rows=[]
for scenario,d in sem.groupby("scenario",sort=False,dropna=False):
    rows.append([scenario,len(d),int(d.converged.sum()),int(d.formal_se_available.sum()),int(d.source_admissible.sum()),int(d.inference_valid.sum())])
sem_table=md(["Model","Fits attempted","Converged","Formal sandwich SEs available","Admissible source moments","Regular local inference"],rows)
coverage=[]
for (scenario,order),d in mom.groupby(["scenario","order"],sort=False,dropna=False):
    each=d.groupby("key").covered.mean()
    coverage.append([scenario,int(order),d.run.nunique(),f"{d.covered.mean():.1%}",f"{each.min():.0%}–{each.max():.0%}"])
cov_table=md(["Model","Cumulant order","Replicates","Pooled entry coverage","Per-entry range"],coverage)
paths=[]
summary_paths=[]
for (scenario,parameter),d in par[par.parameter.isin(["ab","bc","a"])].groupby(["scenario","parameter"],sort=False):
    truth=float(d.truth.iloc[0]); est=d.estimate.to_numpy(float);se=d.se.to_numpy(float)
    # Factor sign is arbitrary: align to the positive-loading generating convention.
    if parameter=="a":est=abs(est)
    good=np.isfinite(se)
    cov=float(np.mean(abs(est[good]-truth)<=1.96*se[good])) if good.any() else None
    paths.append([scenario,parameter,f"{truth:.3f}",f"{est.mean():.3f}",f"{est.std(ddof=1):.3f}",int(good.sum()),"—" if cov is None else f"{cov:.0%}"])
    summary_paths.append(dict(scenario=scenario,parameter=parameter,truth=truth,mean=float(est.mean()),
        empirical_sd=float(est.std(ddof=1)),n=len(est),n_reported_se=int(good.sum()),conditional_coverage=cov))
path_table=md(["Model","Path","Truth","Mean estimate","Across-run SD","Formal SEs","Coverage among available formal SEs"],paths)
oracle=[]
for scenario,d in exact.groupby("scenario",sort=False):
    p=d[d.parameter.isin(["ab","bc","a"])]
    vals=", ".join(f"{x.parameter}={x.estimate:.10f}" for x in p.itertuples()) or "No structural path"
    oracle.append([scenario,f"{d.Q.iloc[0]:.3g}",vals])
oracle_table=md(["Model","Exact-target discrepancy","Recovered path(s)"],oracle)
toggle_table=md(["Model","Free higher-cumulant Q","Gaussian-source Q","Covariance-only rank","Covariance-only SEs"],
    [[x.scenario,f"{x.nonGaussian_Q:.3f}",f"{x.Gaussian_Q:.3f}",f"{x.covariance_only_rank}/{x.covariance_only_parameters}",x.covariance_only_inference_valid] for x in toggles.itertuples()])
text=f"""# Executed validation

The software recovers the known models from their exact population moments,
but this Monte Carlo finds **substantial undercoverage of higher-order moment
intervals**. In the 20 chain runs, pooled-entry coverage was 95.8% for order two,
89.5% for order three and **78.0% for order four** (nominal 95%). The fourth-order
per-entry range was 40–90%. These results do not support claiming that the
higher-order intervals are generally calibrated.

These are execution results for the checked-out implementation. The complete
replication-level estimates, failure status, formal sandwich standard errors and
source-moment warnings are in `validation/results/`.

## Design

The experiment calls **GWASBrewer 0.3.0.0259** directly. Its bundled HapMap3
chromosome-19 reference supplies 39 real LD patterns with **140–863 SNPs per
block**. Each pattern is repeated twice with independent simulated effects and
sampling errors: **78 blocks and 38,980 SNPs**. Both signed LD and allele
frequencies are recoded to the minor allele before generating asymmetric effects.
Scoring uses the exact generating matrices; finite-reference correction is off.
Reference-panel uncertainty is therefore outside this experiment.

There are three observed traits in every run. The primary experiment has
**N={meta['sample_size']:,} per trait, 50% pairwise participant overlap**, and
observational correlations 0.2. Source effects are zero with probability 0.85;
otherwise they are centered, standardized gamma(shape=2), with the scale set by
the direct genetic variance. True source means are zero. All zero-effect SNPs
remain in the analysis. The seed for run k is 20261103+k.

- **Chain:** 20 independent runs, A→B=0.4 and B→C=0.25.
- **Null:** 5 independent runs, no structural paths.
- **Factor:** 5 independent runs, an unobserved non-Gaussian F with equal direct
  loadings 0.4 onto A, B and C. All observed disturbances are also non-Gaussian.

The Python stage estimates all **34 unique cumulants through order four** and
the full 34×34 sampling covariance. Overlap intercepts and every lower-order
quantity are refitted in all 78 whole-LD-block deletions. Unequal deletion sizes
use the documented weighted-pseudovalue covariance. The SEM stage fits the
31 second- through fourth-order cumulants with DWLS and full-V sandwich SEs.
All moments are rescaled by one deterministic common factor before fitting, so
the structural paths retain their generating values.

## Population-target check

True targets are computed analytically from the generating source distributions
and propagation matrix. They are distinct from the empirical cumulants of each
finite set of generated joint SNP effects; both appear in the replication table.
Fitting exact population targets tests the SEM implementation independently of
the noisy LD regressions. Identity weighting in this numerical check has no
sampling interpretation, so its arbitrary SEs are not reported.

{oracle_table}

## Estimated-moment and SEM results

{sem_table}

{path_table}

The primary fit is **unconstrained DWLS**. It selects the minimum converged
objective without preferring a solution because its source moments are
admissible. Negative source variance or impossible higher moments are retained
and prominently flagged as classical or higher-order Heywood cases. All returned
point estimates contribute to the mean and across-run SD above; there is no
admissibility-based filtering or truncation.

Formal sandwich SEs are retained when the local numerical calculation is regular,
including in flagged Heywood cases. They are marked separately from regular
model-based inference. The last coverage column conditions only on numerical SE
availability; it is **not unconditional 95% coverage** or a claim that flagged
fits define a valid source distribution. Nonconvergence and rank failures remain
separate from source admissibility.

## Moment interval coverage

{cov_table}

Pooled-entry coverage averages many correlated moment entries; it is descriptive,
not additional independent Monte Carlo evidence. With only 20 runs for the chain
and five for each other model, these rates have substantial Monte Carlo error.
Higher moments can be unstable even when second-order estimates are satisfactory.

## Observed-order and Gaussian-source controls

The first replicate of each scenario was refitted with source third/fourth
cumulants fixed to zero, and separately with only observed second-order moments
in the fitting objective. The covariance-only fit deliberately retains the
unused source higher-cumulant parameters; rank diagnostics must flag them.

{toggle_table}

Q is an uncalibrated discrepancy. These differences are not chi-square tests.

## Stress test and remaining limits

An initial high-signal design used N=1,000,000 with the same 38,980 SNPs. Its
estimated sampling-error intercept matrix was indefinite (minimum eigenvalue
−6.75886 in the first replicate), and the pipeline stopped with an explicit
error rather than projecting it silently. The same data ran when supplied the
known generating overlap covariance. This illustrates the difficulty of
separating an intercept when genetic signal dominates a short reference genome;
it does not justify replacing a failed empirical intercept by an arbitrary one.

The real-LD Monte Carlo also exposed impossible unconstrained source moments.
The retained `legacy_feasibility_selection/` outputs are a historical diagnostic
from 29 runs with the earlier candidate-selection policy; they are not the
primary results. The main 30-run tables use the current unconstrained minimum-Q
selection and retain the Heywood estimates. Optional feasibility penalties are
available as a distinct sensitivity estimator, not as the default. Diagnostic
moment-feasibility checks are necessary conditions, not a proof that every
multivariate distribution exists. Local numerical identification does not identify a causal
direction or remove equivalent source assignments. Heavy-tailed finite-sample
behavior, reference mismatch, binary traits, much larger models, and broad
coverage calibration remain separate validation tasks.

## Reproduce

Install the Python project and the R package, plus the validation-only Python
dependency `pandas` and R package `GWASBrewer`. Run from the repository root:

```bash
Rscript validation/generate_gwasbrewer.R validation/work 20 50000
python validation/run_moments.py
Rscript validation/fit_models.R
python validation/summarize.py
```

The generator writes several hundred MB of reproducible intermediate inputs to
the ignored `validation/work/` directory. Source scripts and small tables are
tracked; LD matrices, raw simulated GWAS and runtime installations are not.
GWASBrewer reference and software details are available at
<https://jean997.github.io/GWASBrewer/>.
"""
fourth_oracle_file=r/"fourth_order_oracle.json"
if fourth_oracle_file.exists():
    fourth_oracle=json.loads(fourth_oracle_file.read_text())
    fourth_section=f"""
## Independent fourth-order population identity check

The chain's 20-run mean B fourth-cumulant estimate was only 0.199 times its
population target. This prompted an independent experiment with
{fourth_oracle['replicates']:,} fresh source draws through the same real LD
blocks. It used no GWAS sampling noise and supplied the exact population mean
and variance. It therefore isolates the fourth-order LD-regression identity
from overlap estimation and nonlinear nuisance fitting.

The oracle estimate averaged **{fourth_oracle['mean_to_target']:.3f} times the
target**, with Monte Carlo SE **{fourth_oracle['monte_carlo_se_to_target']:.3f}
times the target**, consistent with the stated population identity. The
distribution was strongly right-skewed: its median was
{fourth_oracle['median_to_target']:.3f} times the target, SD was
{fourth_oracle['sd_to_target']:.3f} times the target, and its 99th percentile
was {fourth_oracle['quantile_to_target'][-1]:.3f} times the target. The directly
sampled source fourth-moment check averaged
{fourth_oracle['mean_realized_source_oracle_to_target']:.3f} times the target.

Thus the short experiment's apparent attenuation is not evidence of a
systematic 0.2 multiplicative bias in the LD-power identity. Rare, large
LD-mixed fourth-product fluctuations can make a small Monte Carlo mean
misleading. This result does **not** validate the complete fitted-nuisance
estimator or repair its observed interval undercoverage; those remain separate
calibration questions. No estimate or replicate was excluded on sign or
source-admissibility grounds.

Reproduce this supplementary check after generating the LD inputs and scores:

```bash
python validation/check_fourth_order_oracle.py --workdir validation/work --replicates {fourth_oracle['replicates']} --seed {fourth_oracle['seed']}
```

The compact JSON summary and every oracle estimate are retained in
`validation/results/fourth_order_oracle.*`; the command above names its source
script.
"""
    text=text.replace("\n## Observed-order and Gaussian-source controls\n",fourth_section+"\n## Observed-order and Gaussian-source controls\n")
(root/"docs/validation.md").write_text(text)
(r/"summary.json").write_text(json.dumps(dict(paths=summary_paths,n_moment_runs=len(run),
    moment_failures=int((run.status!="ok").sum()),sem_runs=len(sem),
    sem_converged=int(sem.converged.sum()),sem_formal_se_available=int(sem.formal_se_available.sum()),
    sem_source_admissible=int(sem.source_admissible.sum()),sem_regular_local_inference=int(sem.inference_valid.sum())),indent=2)+"\n")
print(text)
