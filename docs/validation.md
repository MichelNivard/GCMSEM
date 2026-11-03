# Executed validation

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
**N=50,000 per trait, 50% pairwise participant overlap**, and
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

| Model | Exact-target discrepancy | Recovered path(s) |
|---|---|---|
| chain | 7.55e-16 | ab=0.4000000002, bc=0.2500000001 |
| null | 0 | No structural path |
| factor | 2.24e-21 | a=0.4000000000 |

## Estimated-moment and SEM results

| Model | Fits attempted | Converged | Formal sandwich SEs available | Admissible source moments | Regular local inference |
|---|---|---|---|---|---|
| chain | 20 | 18 | 18 | 4 | 3 |
| null | 5 | 5 | 5 | 2 | 2 |
| factor | 5 | 4 | 4 | 0 | 0 |

| Model | Path | Truth | Mean estimate | Across-run SD | Formal SEs | Coverage among available formal SEs |
|---|---|---|---|---|---|---|
| chain | ab | 0.400 | 0.406 | 0.139 | 18 | 83% |
| chain | bc | 0.250 | 0.262 | 0.150 | 18 | 89% |
| factor | a | 0.400 | 0.400 | 0.065 | 4 | 100% |

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

| Model | Cumulant order | Replicates | Pooled entry coverage | Per-entry range |
|---|---|---|---|---|
| chain | 1 | 20 | 98.3% | 95%–100% |
| chain | 2 | 20 | 95.8% | 90%–100% |
| chain | 3 | 20 | 89.5% | 80%–95% |
| chain | 4 | 20 | 78.0% | 40%–90% |
| null | 1 | 5 | 93.3% | 80%–100% |
| null | 2 | 5 | 93.3% | 80%–100% |
| null | 3 | 5 | 98.0% | 80%–100% |
| null | 4 | 5 | 93.3% | 80%–100% |
| factor | 1 | 5 | 100.0% | 100%–100% |
| factor | 2 | 5 | 86.7% | 60%–100% |
| factor | 3 | 5 | 92.0% | 80%–100% |
| factor | 4 | 5 | 89.3% | 60%–100% |

Pooled-entry coverage averages many correlated moment entries; it is descriptive,
not additional independent Monte Carlo evidence. With only 20 runs for the chain
and five for each other model, these rates have substantial Monte Carlo error.
Higher moments can be unstable even when second-order estimates are satisfactory.

## Independent fourth-order population identity check

The chain's 20-run mean B fourth-cumulant estimate was only 0.199 times its
population target. This prompted an independent experiment with
1,000 fresh source draws through the same real LD
blocks. It used no GWAS sampling noise and supplied the exact population mean
and variance. It therefore isolates the fourth-order LD-regression identity
from overlap estimation and nonlinear nuisance fitting.

The oracle estimate averaged **0.930 times the
target**, with Monte Carlo SE **0.075
times the target**, consistent with the stated population identity. The
distribution was strongly right-skewed: its median was
0.396 times the target, SD was
2.372 times the target, and its 99th percentile
was 9.579 times the target. The directly
sampled source fourth-moment check averaged
1.001 times the target.

Thus the short experiment's apparent attenuation is not evidence of a
systematic 0.2 multiplicative bias in the LD-power identity. Rare, large
LD-mixed fourth-product fluctuations can make a small Monte Carlo mean
misleading. This result does **not** validate the complete fitted-nuisance
estimator or repair its observed interval undercoverage; those remain separate
calibration questions. No estimate or replicate was excluded on sign or
source-admissibility grounds.

Reproduce this supplementary check after generating the LD inputs and scores:

```bash
python validation/check_fourth_order_oracle.py --workdir validation/work --replicates 1000 --seed 276351
```

The compact JSON summary and every oracle estimate are retained in
`validation/results/fourth_order_oracle.*`; the command above names its source
script.

## Observed-order and Gaussian-source controls

The first replicate of each scenario was refitted with source third/fourth
cumulants fixed to zero, and separately with only observed second-order moments
in the fitting objective. The covariance-only fit deliberately retains the
unused source higher-cumulant parameters; rank diagnostics must flag them.

| Model | Free higher-cumulant Q | Gaussian-source Q | Covariance-only rank | Covariance-only SEs |
|---|---|---|---|---|
| chain | 19.274 | 25.999 | 5/11 | False |
| null | 24.566 | 39.313 | 3/9 | False |
| factor | 6.709 | 29.175 | 4/12 | False |

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
