\name{fit_measures}
\alias{fit_measures}
\title{Model discrepancy and approximate uncertainty}
\usage{fit_measures(fit)}
\arguments{\item{fit}{A fitted \code{gcmsem_fit} object.}}
\description{Returns the minimized DWLS discrepancy Q, an approximate standard
error, and an exploratory null reference accounting for fitted parameters and
the full supplied moment sampling covariance.}
\value{A one-row data frame with \code{Q}, \code{Q_se}, \code{Q_se_delta},
\code{Q_null_mean}, \code{Q_null_sd}, \code{Q_null_95},
\code{Q_reference_scale}, \code{Q_reference_df}, \code{Q_p_approx},
\code{residual_df} and \code{Q_status}. It is also stored in
\code{fit$fit_measures}.}
\details{Let e be residuals divided by their moment standard errors, J the
correspondingly scaled model Jacobian, C the moment sampling correlation and
P the orthogonal projection off the columns of J. Put Cr = P C P.
The local Gaussian quadratic approximation gives a null mean tr(Cr) and
variance 2 tr(Cr squared). The reported squared SE adds 4 e' Cr e to that null
variance. This substitutes observed residuals for unknown population misfit;
it is a rough plug-in sampling SD, and can overestimate uncertainty near exact
fit. Q_se_delta retains only the first-order term, which vanishes at zero
residuals. The null reference matches the first two moments to a scaled
chi-square distribution. Its approximate upper-tail probability is
Q_p_approx; its 95th percentile is Q_null_95.

V, fitting weights and the local model tangent are held fixed. The calculation
does not propagate uncertainty in V, reference LD, model selection or
nonlinear curvature under substantial misspecification. It is not a validated
WLSMV correction or a causal-direction test. Do not construct a normal
confidence interval for Q from Q_se or subtract separate Q variances for
comparisons of dependent models. Penalized/constrained fits, nonconvergence
and deficient local parameter rank return NA uncertainty. Unconstrained
Heywood solutions retain formal values with a flagged status. A saturated
tangent space has zero local variance and no nondegenerate reference p-value.}
\seealso{\code{\link{gcmsem}}, \code{\link{parameter_estimates}}}
