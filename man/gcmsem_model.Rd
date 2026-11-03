\name{gcmsem_model}
\alias{gcmsem_model}
\alias{print.gcmsem_model}
\title{Define a genomic comoment SEM}
\usage{
gcmsem_model(syntax, observed, latent = character(),
  gaussian_observed = TRUE, gaussian_latent = TRUE)
\method{print}{gcmsem_model}(x, ...)
}
\arguments{
\item{syntax}{Single string with newline- or semicolon-separated model statements.}
\item{observed}{Unique observed trait names, in moment-key ordering.}
\item{latent}{Additional latent names; factor names on the left of \code{=~} are inferred.}
\item{gaussian_observed, gaussian_latent}{Logical scalar or names in the respective node class to keep Gaussian by default. FALSE frees marginal third and fourth source cumulants; explicit syntax overrides defaults.}
\item{x}{Model object.}
\item{...}{Unused additional print arguments.}
}
\description{Parses a small lavaan-inspired syntax into a sparse linear SEM with
source cumulants and shared parameter labels. No lavaan or MCMSEM dependency.}
\details{Operators are \code{=~} for loadings, \code{~} for regressions and
\code{~~} for source covariances. \code{x ~~~ x*y} is the source cumulant xxy;
\code{x ~~~~ x*y*y} is xxyy. Optional aliases \code{x ~-- y} and
\code{x ~~- y} denote xxy and xyy. Numeric modifiers fix values; labels tie
parameters. Latent variances default to one, observed source variances are free,
and other source cross-cumulants default to zero. Gaussian defaults affect source
disturbances, not the resulting observed marginal distribution. Mean/intercept,
threshold, group and general nonlinear-constraint syntax are not supported.}
\value{A \code{gcmsem_model} containing the parsed sparse parameter table,
observed/latent node names, bounds and default starts.}
\examples{
m <- gcmsem_model('f =~ a*x + a*y; y ~ b*x', c('x','y'),
                  gaussian_observed=FALSE, gaussian_latent=FALSE)
print(m)
}
