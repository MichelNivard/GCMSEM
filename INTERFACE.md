# Development interface contract

R package at repository root (`DESCRIPTION`, `R/`, `man/`, `tests/`). Python package `python/gcmsem` with root `pyproject.toml`; CLI `gcmsem`.

Canonical moment key: trait names sorted in declared trait order, separated by `|` (e.g. `Height|Height|Weight`). Orders 1–4; fourth order is excess cumulant, third central cumulant. Input scale is genotype-standardized effects; per-allele inputs converted using reference 2p(1-p), tracked in metadata. Counts include zero-effect variants. LD-reference universe and regression/tag variants must be distinguished.

Interchange output PREFIX.moments.tsv: columns `key`, `order`, `estimate`, `se`. PREFIX.sampling_cov.tsv is a square table with first column `key` and header moment keys. PREFIX.metadata.json includes schema_version `1.0`, traits ordered array, effect_scale, moment_type `cumulant`, and all provenance. PREFIX.jackknife.tsv has block ID and moment columns. PREFIX.covariance.tsv and PREFIX.correlation.tsv are trait square tables. Optional tensor arrays only on explicit request. R read_moments(prefix) reads this bundle.

Score NPZ contract to coordinate between agents: `snp` string array n; `a1`, `a2` string arrays n (LD coding alleles); `af` float n; `chrom`, `bp`, `block` arrays n; `L1`, `L2`, `L3`, `L4` float n x C; `counts` float C (reference/source variants per annotation); metadata JSON string `metadata`. Signed LD powers, not abs; strata disjoint. Scoring agent provides `load_scores(path)` validated dictionary. Estimator takes this format, filters tag rows without changing source counts. Full pipeline refits all orders and noise calibration in genomic jackknife; no significance pruning.

R API: gcmsem_model(syntax, observed, ...), implied_moments(model, par=NULL, orders=2:4, jacobian=FALSE), gcmsem(model, moments, V=NULL, ...), read_moments(prefix), parameter_estimates(fit), summary/coef/vcov/fitted/print. If need deviations notify other agents before coding.

R syntax: `y ~ b*x`; `f =~ a*x + a*y`; `x ~~ v*x`; repeated label ties parameters. Canonical higher cumulant syntax `x ~~~ x*y` means kappa(x,x,y), `x ~~~~ x*y*y` means kappa(x,x,y,y). Add explicit aliases `x ~-- y`=xxy and `x ~~- y`=xyy if practical; canonical tuple notation must cover all N-trait combinations. Gaussian default sets source cumulants of orders 3 and 4 to zero, explicit syntax can free them. Controls distinguish selecting observed orders for fitting from fixing source distributions. No MCMSEM runtime dependency.

Default estimator unconstrained DWLS, full V robust sandwich without N division, multistart nlminb with exact structured forward-mode derivatives. Retain classical/higher-order Heywood estimates and computable formal sandwich SEs with prominent warnings and separate source-admissibility/inference flags. Do not preferentially select admissible source moments or clip free variances by default. Identification or numerical failures must never yield fabricated SEs. Optional feasibility penalties are an explicitly different estimator. No WLSMV or calibrated chi-square p-values claimed. Preserve reciprocal stability and source permutation limitations.

Do not commit large GWAS, licensed reference data, runtime binaries, or old empirical archives. Commit small reproducible synthetic fixtures/results and scripts; validation reports must state what actually ran.
