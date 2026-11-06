\name{estimate_polygenicity}
\alias{estimate_polygenicity}
\alias{print.gcmsem_polygenicity}
\alias{summary.gcmsem_polygenicity}
\alias{print.summary_gcmsem_polygenicity}
\alias{as.data.frame.gcmsem_polygenicity}
\alias{coef.gcmsem_polygenicity}
\alias{plot.gcmsem_polygenicity}
\title{Estimate Effective Polygenicity from GWAS Summary Statistics}
\description{Estimate second and fourth source cumulants and a Gaussian-mixture
effective variant count. Each trait is harmonized and analysed independently.
Returns an R object with uncertainty, intermediate estimates, source strata,
genomic deletions, influence diagnostics and logs. Requires the accompanying
Python GCMSEM engine and jsonlite. No comparator software is required.}
\usage{
estimate_polygenicity(sumstats, scores, effect_scale = NULL,
  traits = NULL, n_blocks = 200L, noise_intercept = NULL,
  min_info = NULL, min_n = NULL, af_tolerance = 0.1,
  allow_strand_flips = FALSE, keep_palindromic = FALSE,
  max_memory_mb = 2048, max_condition = 1e6, python = Sys.which("python"),
  write.to.tsv = FALSE, out = NULL, overwrite = FALSE, verbose = FALSE)
\method{print}{gcmsem_polygenicity}(x, digits = 5L, ...)
\method{summary}{gcmsem_polygenicity}(object, trait = NULL, digits = 5L, ...)
\method{print}{summary_gcmsem_polygenicity}(x, ...)
\method{as.data.frame}{gcmsem_polygenicity}(x, row.names = NULL, optional = FALSE, ...)
\method{coef}{gcmsem_polygenicity}(object, ...)
\method{plot}{gcmsem_polygenicity}(x, trait = x$traits[1L], ...)
}
\arguments{
\item{sumstats}{A manifest path or data frame with trait, path and effect_scale
columns; a named character vector of GWAS paths; or a named list of GWAS data
frames. Tables require SNP, A1, A2, BETA and SE. EAF, INFO and N are optional.
Manifest paths are resolved relative to the manifest file; data-frame paths
are relative to the R working directory.}
\item{scores}{Path to a GCMSEM signed LD-power score NPZ file with orders one
through four, source counts and disjoint MAF-stratum metadata.}
\item{effect_scale}{Required for named files or named data frames: either
\code{"allele"} or \code{"genotype_standardized"}, or a named vector for all
input traits. Omit when the manifest supplies the scales. Allele-scale effects
and SEs are multiplied by the reference genotype SD. Phenotype scaling is unchanged.}
\item{traits}{Optional character vector selecting and ordering input traits.
NULL selects all. Each trait uses its own eligible SNPs.}
\item{n_blocks}{Maximum number of chromosome-respecting genomic deletion groups.
Default 200. Whole score blocks are retained and no final rows are trimmed.}
\item{noise_intercept}{NULL estimates the sampling-noise calibration from the LD
regression. Alternatively supply one nonnegative number for all selected traits,
or a named numeric vector with one value for each selected trait. A value of
one assumes sampling variance equals the squared supplied GWAS SE.}
\item{min_info,min_n}{Optional minimum INFO and sample-size thresholds.
When requested, the corresponding GWAS column is required. NULL adds no threshold.}
\item{af_tolerance}{Maximum absolute study/reference allele-frequency difference
when EAF is supplied. Default 0.1; must be in [0,0.5).}
\item{allow_strand_flips}{Allow strand-complement matching. Default FALSE.}
\item{keep_palindromic}{Retain A/T or C/G SNPs only when EAF uniquely resolves
their orientation. Default FALSE drops them.}
\item{max_memory_mb}{Maximum estimated working memory per trait, in MiB.}
\item{max_condition}{Maximum permitted regression-design condition number.}
\item{python}{Executable from an environment with the current Python GCMSEM
engine installed. Virtual-environment symlinks are preserved.}
\item{write.to.tsv}{FALSE returns an R object without persistent output files.
TRUE also calls \code{export_gcmsem()}. Temporary engine files are cleaned up.}
\item{out}{Output filename prefix, required only for \code{write.to.tsv=TRUE}.}
\item{overwrite}{Allow replacement of files managed by the export. Default FALSE.}
\item{verbose}{Print the captured harmonization, estimation and diagnostic log.}
\item{x,object}{A \code{gcmsem_polygenicity} object; for the summary print method,
a \code{summary_gcmsem_polygenicity} object.}
\item{trait}{Traits to display in summary; NULL displays all. The plot accepts
one trait and defaults to the first.}
\item{digits}{Significant digits for printed output.}
\item{row.names,optional}{Standard data-frame conversion arguments; optional is unused.}
\item{...}{Additional graphics arguments for plot; unused by the other methods.}
}
\details{
For source counts \eqn{M_c}, variances \eqn{v_c} and fourth cumulants \eqn{k_c},
the count is \eqn{3H^2/D}, where \eqn{H=\sum_c M_c v_c} and
\eqn{D=\sum_c M_c k_c+3\sum_c M_c v_c a_c}.
\eqn{a_c} is mean second-order tag LD to each source stratum multiplied by
the estimated source variances. Tag means approximate source means.
The main SE recomputes all terms, including these means, in every deletion,
then applies the same unequal-group jackknife as the moment estimator.
The secondary delta SE uses the full joint variance/fourth-cumulant covariance
with mean LD held fixed. It is a sensitivity diagnostic, not a replacement SE.

The count without LD uses \eqn{\sum_c M_c(k_c+3v_c^2)} as denominator.
All counts describe concentration of genetic variance under a zero-mean
independent-source Gaussian variance-mixture model; they are not literal counts
of nonzero effects without additional assumptions. The source SNP universe
determines the target. No sample splitting or p-value threshold is applied.

Warnings flag nonpositive source variances, negative fourth cumulants (which
conflict with the Gaussian-mixture interpretation), violations of
\eqn{k_4\geq k_3^2/v-2v^2}, a denominator within 1.96 jackknife SEs of zero,
nonpositive deletion denominators, undefined ratios, counts outside the source
range and a deletion contributing over half the estimated count variance.
No clipping or deletion exclusion occurs. Undefined ratios have NA uncertainty;
the corresponding denominator and all deletion records remain available.
Symmetric confidence intervals are not provided because ratio uncertainty can
be strongly nonnormal. Reference-panel uncertainty is not included and
higher-order coverage remains under validation.
}
\value{A \code{gcmsem_polygenicity} list. \code{estimates} is a per-trait data
frame containing effective_count, se, delta_se, inverse_count and its SE,
no_ld_count and its SE, tag/source/block counts, and flags.
\code{details[[trait]]} contains stages (estimates and SEs), strata (all four
conditional cumulants and SEs), covariance (joint stage covariance),
conditional_covariance (joint variance/fourth-cumulant covariance across bins),
ld_mean, jackknife (coordinates, fractions, all stage estimates, and each
deletion's share of count variance), warnings, and assumptions.
\code{moments[[trait]]} retains the pooled and conditional moment objects,
deletions, sampling covariance, QC metadata and noise intercept.
\code{log}, \code{metadata}, \code{traits}, and \code{call} retain provenance.
Use \code{saveRDS()} to preserve the complete object. Export writes readable
tables and moment bundles. No between-trait covariance is estimated.}
\seealso{\code{\link{estimate_moments}}, \code{\link{export_gcmsem}}}
\examples{
\dontrun{
p <- estimate_polygenicity(
  c(Depression = "tutorial-data/sumstats/20544_11.tsv.gz"),
  "tutorial-data/scores/eur_hm3.scores.npz",
  effect_scale = "genotype_standardized")
p
summary(p)
p$details$Depression$stages
plot(p, trait = "Depression")
saveRDS(p, "depression-polygenicity.rds")
}
}
