\name{read_moments}
\alias{read_moments}
\title{Read a versioned Python moment bundle}
\usage{read_moments(prefix)}
\arguments{\item{prefix}{Path prefix of .moments.tsv, .sampling_cov.tsv and .metadata.json.}}
\description{Validates ordered keys, SE/covariance agreement, cumulant convention
and schema version, then returns an R object ready for \code{gcmsem()}.}
\details{When optional jsonlite is installed, complete JSON metadata are parsed.
Otherwise complete metadata are preserved as raw JSON and the required convention
and version declarations are validated. Bundles containing R scaling metadata
require \pkg{jsonlite} so that the trait scales can be restored and validated.
No arbitrary JSON evaluation is performed. Conditional annotation bundles are
read separately by their own prefix; use \code{saveRDS()} and \code{readRDS()}
to save and restore an entire R object including all conditional objects.}
\value{A \code{gcmsem_moments} list with \code{moments} (a data frame),
\code{V} (the full sampling covariance of the estimated cumulants),
\code{S} and \code{S_SE} (SNP-effect covariance and SE matrices),
\code{R} and \code{R_SE} (genetic correlation and delta-method SE matrices),
metadata, trait names, diagnostics, the input prefix, any stored trait
\code{scaling}, and any stored jackknife
replicates and log. Matrix rows and columns use the trait names. Correlation
SEs use the full covariance among the relevant second-order moments.
Correlations involving nonpositive marginal variances are \code{NA}; all original
cumulant estimates remain available. Printing gives a readable overview;
\code{summary()} and \code{as.data.frame()} support further interactive inspection.}
