\name{export_gcmsem}
\alias{export_gcmsem}
\title{Export a GCMSEM R object to files}
\usage{export_gcmsem(x, prefix, overwrite = FALSE)}
\arguments{
\item{x}{A \code{gcmsem_fit}, \code{gcmsem_moments}, or \code{gcmsem_polygenicity} object.}
\item{prefix}{Output path prefix, for example \code{"results/height_weight"}. Parent directories are created as needed.}
\item{overwrite}{Allow replacement of existing export files. Obsolete optional
jackknife or log files at the exported prefixes are removed. The default refuses
every collision before writing.}
}
\description{Writes explicitly requested tab-separated tables and supporting metadata while leaving the supplied R object unchanged.}
\details{Moment objects export the cumulant table, complete sampling covariance,
SNP-effect covariance and correlation matrices with SEs, JSON metadata, and any
stored jackknife matrix or log. Named conditional bundles are exported with
additional prefix suffixes. The primary bundle can be read again using
\code{read_moments(prefix)}; conditional bundles can be read by their respective
prefixes. Structured metadata export, including updated trait scaling metadata,
uses \pkg{jsonlite}.
Polygenicity objects export the summary, intermediate estimates, source strata,
full stage and conditional covariance, LD averages, deletion influences,
logs, metadata and complete per-trait moment bundles.
Fit objects export parameter estimates, observed and fitted moments, residuals,
sampling covariance, parameter covariance, restart results, diagnostics,
model syntax and any stored log. Numerical tables retain seventeen significant
digits. Use \code{saveRDS()} and \code{readRDS()} to save and restore the complete R object.}
\value{Invisibly, a character vector of the paths written.}
\examples{
m <- gcmsem_model('y ~ b*x', c('x','y'))
v <- implied_moments(m)
V <- diag(0.01, length(v)); dimnames(V) <- list(names(v), names(v))
fit <- gcmsem(m, v, V, orders=2, n_starts=1)
files <- export_gcmsem(fit, tempfile('gcmsem-'))
unlink(files)
}
