# GCMSEM 0.3.0

* Added `fit_measures()` with an approximate SE for the minimized DWLS discrepancy Q, a first-order component, and a local Gaussian moment-matched null reference. Calculations use the full moment sampling covariance and project out fitted parameters. Approximate p-values are explicitly labeled; penalized or numerically irregular fits return unavailable uncertainty. Unconstrained Heywood fits retain formal diagnostics with warnings.
* Fit summaries, in-memory logs and optional TSV exports include the new Q diagnostics. The fitting objective and parameter estimates are unchanged.
* Reworked the README and Wiki for first-time users, added model diagrams beside all six causal specifications, and moved detailed R syntax and result handling to dedicated pages. Model-fit assessment precedes effect interpretation.
* The correlation tutorial now estimates the full 15-trait matrix from GWAS files, runs original standard LDSC, reads its results into R, and renders the comparison. The standard LDSC runner can complete munging and fitting without requiring a precomputed GCMSEM result bundle.
* Published the numerical data at Zenodo DOI 10.5281/zenodo.23977299 and activated checksum-verified downloads.

# GCMSEM 0.2.0

- `estimate_moments()` runs the Python summary-statistics engine from R and
  returns a `gcmsem_moments` object. It accepts manifests, named GWAS paths,
  or named lists of GWAS data frames.
- Moment objects include named covariance and correlation matrices, their
  standard errors, the full sampling covariance, genomic deletion estimates,
  conditional moment bundles, provenance, and the estimation log.
- Fit and moment objects have readable print and summary methods. Assigning
  a fit or summary keeps the result in R; `write.to.tsv = TRUE` and
  `export_gcmsem()` provide explicit table export.
- `standardize_moments()` transforms cumulants, their joint covariance, and
  deletion estimates using fixed trait scales and records the scale factors.
- Worked tutorials now include complete R scripts, file descriptions,
  executed console output, and a scatterplot of 105 GCMSEM/standard-LDSC
  correlation pairs. An original Neale Lab GWAS preparation script provides
  a source-download route for the causal examples.
- Unconstrained estimates and computable formal standard errors remain
  available for Heywood cases, with prominent warnings and separate
  convergence, identification, source-admissibility, and inference flags.

The accompanying Python numerical engine remains version 0.1.0.

# GCMSEM 0.1.0

Initial standalone implementation of Genomic Comoment SEM.

* Explicit higher-cumulant declarations separate variable tuples from quoted parameter labels: `k3(x, x, y) = free("s")` and `k4(x, x, y, y) = 0`. Directional punctuation aliases are retired; older product notation emits migration warnings.
* GitHub-compatible math delimiters render equations in the README and mathematical documentation.
* Python signed LD scores, summary-statistic harmonization, Gaussian sampling-noise correction, cumulant regression, multivariate block jackknife, and versioned R interchange.
* Standalone R SEM syntax, source cumulant propagation, exact derivatives, DWLS estimation, and sandwich standard errors with numerical and identification diagnostics.
* Reproducible GWASBrewer validation and separate software tests.

This is the first research release. Passing implementation tests does not establish universal inferential calibration, robustness to misspecified LD, or causal identification.

