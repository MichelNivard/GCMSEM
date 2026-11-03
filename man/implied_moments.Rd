\name{implied_moments}
\alias{implied_moments}
\alias{moment_keys}
\alias{source_diagnostics}
\title{Compute model cumulants, exact derivatives and necessary feasibility checks}
\usage{
implied_moments(model, par = NULL, orders = 2:4, jacobian = FALSE,
  keys = NULL, backend = c("auto", "R", "native"))
moment_keys(traits, orders = 2:4)
source_diagnostics(model, par = NULL, joint_max_nodes = 12L,
  tolerance = 1e-8)
}
\arguments{
\item{model}{A \code{gcmsem_model}.}
\item{par}{All free parameter values, optionally named. NULL uses model defaults.}
\item{orders}{Cumulant orders. The model computes orders 2 to 4; moment_keys also permits order 1.}
\item{jacobian}{Return the exact first derivative matrix as well as predictions.}
\item{keys}{Optional selected canonical keys, e.g. \code{c("x|x", "x|x|y")}.}
\item{backend}{Use the compiled tensor contraction by default, the independent R implementation, or require the native kernel.}
\item{traits}{Ordered unique trait names.}
\item{joint_max_nodes}{Maximum number of non-deterministic source nodes for full degree-two moment-matrix PSD checks.}
\item{tolerance}{Relative numerical tolerance for necessary PSD and Pearson checks.}
}
\description{Propagates sparse symmetric source cumulants through the linear
structural system. Computes exact structured first derivatives, including through
matrix inversion; numerical differencing is not used by the estimator.}
\details{Only unique observed moment tuples are computed, avoiding full fourth-order
arrays. Source moments are checked using nonnegative variances, covariance PSD,
marginal Pearson bounds and, for sufficiently small models, the raw moment matrix
for monomials of degree at most two. These conditions are necessary, not sufficient,
for a joint distribution. Consistent deterministic zero sources are excluded from
rank checks. No global-identification guarantee is implied.}
\value{For implied_moments, a named numeric vector, or a list with
\code{estimate} and \code{jacobian}. For moment_keys, character keys. For
source_diagnostics, a list of feasibility flags, eigenvalues, margins and limitations.}
\examples{
m <- gcmsem_model('y ~ b*x', c('x','y'), gaussian_observed=FALSE)
implied_moments(m, jacobian=TRUE)
source_diagnostics(m)
}
