\name{standardize_moments}
\alias{standardize_moments}
\title{Rescale an R Moment Object}
\usage{standardize_moments(x, scale = NULL)}
\arguments{
\item{x}{A \code{gcmsem_moments} object returned by \code{estimate_moments()} or \code{read_moments()}.}
\item{scale}{One positive number per trait. Named values are matched to trait names. The default uses the square roots of the estimated second moments.}
}
\description{Divides each trait index of a cumulant by its scale and applies the same transformation to its sampling covariance and available jackknife replicates.}
\value{A new \code{gcmsem_moments} object with rescaled moments and matrices. \code{scaling$scale} contains the divisors for this call; \code{scaling$cumulative_scale} records their product across successive calls. The input object is preserved and no files are written.}
\details{The scales are held fixed when computing standard errors. When scales are estimated from the input moments, this operation conditions on those values. Model restrictions such as equal factor loadings apply on the selected scale. For a fitted path from trait x to trait y, multiply its estimate and standard error by \code{scale[y] / scale[x]} to express that path on the input scale. All conditional annotation objects use the same trait scales as the pooled object.}
\seealso{\code{\link{estimate_moments}}, \code{\link{gcmsem}}, \code{\link{read_moments}}}
