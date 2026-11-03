# GCMSEM 0.1.0

Initial standalone implementation of Genomic Comoment SEM.

* Explicit higher-cumulant declarations separate variable tuples from quoted parameter labels: `k3(x, x, y) = free("s")` and `k4(x, x, y, y) = 0`. Directional punctuation aliases are retired; older product notation emits migration warnings.
* GitHub-compatible math delimiters render equations in the README and mathematical documentation.
* Python signed LD scores, summary-statistic harmonization, Gaussian sampling-noise correction, cumulant regression, multivariate block jackknife, and versioned R interchange.
* Standalone R SEM syntax, source cumulant propagation, exact derivatives, DWLS estimation, and sandwich standard errors with numerical and identification diagnostics.
* Reproducible GWASBrewer validation and separate software tests.

This is the first research release. Passing implementation tests does not establish universal inferential calibration, robustness to misspecified LD, or causal identification.

