# A complete small example

Follow the [installation instructions](README.md#install), then run these commands from the repository directory. The example creates its own inputs and needs no external data download.

This example runs immediately after installation and generates its own inputs. From the repository directory, create three synthetic GWAS tables and their signed LD-power scores:

```bash
python examples/make_toy_inputs.py --out toy
gcmsem score --manifest toy/ld_manifest.tsv --out toy/scores.npz
```

`toy/A.tsv`, `toy/B.tsv` and `toy/C.tsv` contain one SNP per row, with `SNP`, `A1`, `A2`, `BETA`, `SE` and `EAF` columns. `toy/sumstats.tsv` names those files and declares per-allele effects. `toy/scores.npz` holds their LD-power scores; `toy/noise.tsv` contains the known cross-trait sampling-noise intercepts for this simulation.

Now start R in the same directory:

```r
library(GCMSEM)

moments <- estimate_moments(
  "toy/sumstats.tsv",
  scores = "toy/scores.npz",
  noise_intercept = "toy/noise.tsv",
  n_blocks = 80
)
moments

scaled <- standardize_moments(moments)
model <- gcmsem_model('
  B ~ ab*A
  C ~ bc*B
', observed = c("A", "B", "C"), gaussian_observed = FALSE)

fit <- gcmsem(model, scaled, n_starts = 8, seed = 20261103)
summary(fit)
parameter_estimates(fit)

moments$R                 # Complete 3 x 3 genetic-correlation matrix
coef(fit)[c("ab", "bc")]  # Named path estimates
vcov(fit)                 # Parameter sampling covariance
tail(fit$log, 3)          # Retained estimation log
```

The moment object prints:

```text
GCMSEM moments
Traits: A, B, C
SNPs: 3,769
Cumulants: 34 (order 1: 3, order 2: 6, order 3: 10, order 4: 15)
Sampling covariance V: 34 x 34   Jackknife groups: 80

Genetic correlations ($R)
       A      B      C
A 1.0000 0.4095 0.0421
B 0.4095 1.0000 0.3293
C 0.0421 0.3293 1.0000
```

The fitted paths are approximately `ab = 0.3367` and `bc = 0.3001`; `summary(fit)` displays their formal SEs alongside all source-moment parameters and model diagnostics.

`standardize_moments()` divides each trait's SNP effects by its estimated marginal SD and transforms every cumulant and its sampling covariance consistently. The scaling constants are held fixed for fitting and SE calculation; `scaled$scaling$scale` records them. Path estimates from this fit are on that SNP-effect-SD scale.

The toy input has three traits, unequal synthetic LD block sizes, non-Gaussian effects and correlated sampling error. For an ordinary GWAS analysis, omit `noise_intercept` to estimate it from LD regressions. The [GWASBrewer validation](../validation/) separately uses the real HapMap3 chromosome 19 LD pattern and varying block lengths.



Continue with [Model specification](Model-specification.md), [Working with R results](R-results.md), or the [real-data causal examples](Causal-examples.md).
