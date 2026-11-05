\name{estimate_moments}
\alias{estimate_moments}
\title{Estimate genomic moments from GWAS summary statistics}
\description{
Harmonize GWAS summary statistics with signed LD-power scores, estimate
cumulants and their joint sampling covariance, and return an R object.
The Python GCMSEM engine performs the numerical estimation.
}
\usage{
estimate_moments(sumstats, scores, effect_scale = NULL,
  max_order = 4L, n_blocks = 200L, moment_keys = NULL,
  noise_intercept = NULL, min_info = NULL, min_n = NULL,
  af_tolerance = 0.1, allow_strand_flips = FALSE,
  keep_palindromic = FALSE, max_moments = 5000L,
  max_memory_mb = 2048, max_condition = 1e6,
  conditional = TRUE, python = Sys.which("python"),
  write.to.tsv = FALSE, out = NULL, overwrite = FALSE,
  verbose = FALSE)
}
\arguments{
\item{sumstats}{A manifest filename; a manifest data frame with columns
\code{trait}, \code{path}, and \code{effect_scale}; a named character vector of
GWAS filenames; or a named list of GWAS data frames. Trait names start with a
letter and contain letters, digits, underscores or periods.}
\item{scores}{Path to a GCMSEM signed LD-power \code{.npz} score archive.}
\item{effect_scale}{For named files or data frames, either one value for all
traits or a named value per trait. Use \code{"allele"} for per-effect-allele
coefficients, or \code{"genotype_standardized"} for effects per genotype SD.
A manifest supplies this column itself.}
\item{max_order}{Largest cumulant order, 2, 3 or 4. Order 4 denotes excess
cumulants. Means and second-order moments are always estimated.}
\item{n_blocks}{Maximum number of genomic jackknife deletion groups, at least 2.}
\item{moment_keys}{Optional character vector of selected tuples, for example
\code{"Height|Height|Weight"}. The engine also includes required lower-order
tuples and all first- and second-order moments.}
\item{noise_intercept}{Optional labeled TSV filename or numeric square matrix
of sampling-error calibration coefficients. Named matrices are aligned to the
manifest traits; unnamed matrices follow manifest order. The default estimates
these coefficients from LD regressions.}
\item{min_info, min_n}{Optional minimum INFO and per-SNP sample-size thresholds.
Requesting either requires the corresponding GWAS column.}
\item{af_tolerance}{Maximum absolute difference between study and reference
effect-allele frequencies after alignment, in the interval [0, 0.5).}
\item{allow_strand_flips}{Allow complementary allele alignment.}
\item{keep_palindromic}{Retain A/T or C/G variants only when study EAF uniquely
resolves their orientation. The default drops them.}
\item{max_moments}{Maximum number of estimated cumulants.}
\item{max_memory_mb}{Numerical-engine memory budget in MiB.}
\item{max_condition}{Maximum permitted standardized regression-design condition.}
\item{conditional}{Also return annotation-specific moment objects when the
score archive contains multiple source strata.}
\item{python}{Python executable with the accompanying GCMSEM engine installed,
for example the executable inside a virtual environment.}
\item{write.to.tsv}{Write the completed object to disk when \code{TRUE}.
The default returns the result in R and removes temporary interchange files.}
\item{out}{Output filename prefix for explicit export, such as
\code{"results/height_weight"}. Required when \code{write.to.tsv = TRUE}.}
\item{overwrite}{Allow explicit export to replace existing output files.}
\item{verbose}{Print the engine log after estimation. The complete log is
available in the returned object regardless of this setting.}
}
\details{
GWAS tables require SNP identifiers, effect and other alleles, BETA and SE.
The standard column names are \code{SNP}, \code{A1}, \code{A2}, \code{BETA},
and \code{SE}; the Python engine also recognizes its documented aliases.
\code{EAF}, \code{INFO} and \code{N} are optional quality-control columns.
Reference and GWAS variants must use the same genome build.

Paths inside a manifest file are relative to that file's directory. Paths in
an R manifest data frame or named vector are relative to the current working
directory. Named GWAS data frames are serialized only inside a temporary
directory; numeric values are written with 17 significant digits.

The \pkg{jsonlite} package and the Python GCMSEM installation are required
for this entry point. Install the Python engine in the chosen environment
with \code{python -m pip install .} from a GCMSEM checkout. Ordinary R model
fitting with \code{gcmsem()} uses the R package directly.

All temporary engine outputs are read before cleanup. Explicit export occurs
only after successful estimation. A failed engine call raises a
\code{gcmsem_engine_error} condition containing \code{status} and \code{log};
the error message includes the final engine-log lines.
}
\value{
A \code{gcmsem_moments} object. Its \code{moments} data frame contains tuple
keys, orders, estimates and SEs. \code{V} is the joint sampling-covariance
matrix of those estimates. \code{S} and \code{R} contain the SNP-effect
covariance and correlation matrices; \code{S_SE} and \code{R_SE} contain
their standard errors. \code{jackknife} holds the deletion estimates with
moments in columns. \code{conditional} is a named list of source-stratum
moment objects. \code{metadata} records harmonization counts, source hashes,
model assumptions and diagnostics. \code{log}, \code{engine} and \code{call}
retain the numerical log, invocation and R call.

Use \code{print()}, \code{summary()}, \code{as.data.frame()} and direct
component extraction to inspect the result; pass it to \code{gcmsem()} to
fit a model. \code{export_gcmsem()} provides a separate explicit export step.
}
\seealso{\code{\link{gcmsem}}, \code{\link{read_moments}},
\code{\link{export_gcmsem}}}
\examples{
\dontrun{
library(GCMSEM)
files <- c(Height = "height.tsv.gz", Weight = "weight.tsv.gz")
moments <- estimate_moments(
  sumstats = files, scores = "eur_hm3.scores.npz",
  effect_scale = "allele", python = ".venv/bin/python"
)
moments
summary(moments)
moments$S
moments$R
head(moments$moments)

model <- gcmsem_model("Weight ~ b*Height", c("Height", "Weight"))
fit <- gcmsem(model, moments, orders = 2)
summary(fit)

# Save R objects or explicitly export tables after inspecting the result.
saveRDS(moments, "height_weight_moments.rds")
export_gcmsem(moments, "results/height_weight")
}
}
