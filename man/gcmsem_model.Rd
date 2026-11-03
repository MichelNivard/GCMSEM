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
\code{~~} for source covariances. In these statements, numeric modifiers before
\code{*} fix values and labels name free parameters.

Higher source cumulants use explicit variable tuples:
\code{k3(x, x, y) = free("s")} declares the third source cumulant xxy with
parameter label s; \code{k3(x, y, y) = free("t")} instead declares xyy.
\code{k4(x, x, y, y) = free("q")} declares the fourth source cumulant xxyy.
Repeated variable names specify repeated indices. Tuples are symmetric and do
not encode a causal direction. Fourth cumulants are excess cumulants, rather
than fourth central moments.

A standalone \code{k3(x, x, y)} or \code{k3(x, x, y) = free()} declares an
anonymous free parameter. A finite numeric right side fixes the value, for
example \code{k3(f, f, f) = 0} or \code{k4(f, f, f, f) = -1e-6}.
Labels in \code{free("label")} may use either quote style. Variable names and
labels must begin with a letter and contain letters, digits, underscores or
periods. Labels may have the same spelling as variables without ambiguity.
Reusing a label ties parameters on their original scale, including across
different statement types. These declarations occur inside the model string;
they are parsed without evaluating R code.

Latent variances default to one, observed source variances are free, and other
source cross-cumulants default to zero. Gaussian defaults affect source
disturbances, not the resulting observed marginal distribution. Explicit
cumulant declarations override the corresponding defaults. Selecting observed
orders for fitting is separate; use the \code{orders} argument to \code{gcmsem}.

The earlier aliases \code{x ~-- y} and \code{x ~~- y} are rejected with
migration instructions: use \code{k3(x, x, y)} and \code{k3(x, y, y)},
respectively. The older \code{~~~} and \code{~~~~} operators remain temporarily
readable for saved models and emit deprecation warnings with tuple replacements.
Mean/intercept, threshold, group and general nonlinear-constraint syntax are not
supported.}
\value{A \code{gcmsem_model} containing the parsed sparse parameter table,
observed/latent node names, bounds and default starts.}
\examples{
m <- gcmsem_model('
  f =~ a*x + a*y
  y ~ b*x
  k3(f, f, f) = free("skew_f")
  k4(f, f, f, f) = free("fourth_f")
', c('x','y'))
print(m)
}
