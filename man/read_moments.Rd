\name{read_moments}
\alias{read_moments}
\title{Read a versioned Python moment bundle}
\usage{read_moments(prefix)}
\arguments{\item{prefix}{Path prefix of .moments.tsv, .sampling_cov.tsv and .metadata.json.}}
\description{Validates ordered keys, SE/covariance agreement, cumulant convention
and schema version, then prepares a bundle for gcmsem.}
\details{When optional jsonlite is installed, complete JSON metadata are parsed.
Otherwise complete metadata are preserved as raw JSON and the required convention
and version declarations are validated. No arbitrary JSON evaluation is performed.}
\value{A \code{gcmsem_moments} list with moment data, V, metadata, traits and prefix.}
