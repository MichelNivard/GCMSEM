\name{gcmsem}
\alias{gcmsem}
\title{Fit a genomic comoment SEM by DWLS}
\usage{
gcmsem(model, moments, V = NULL, orders = 2:4, keys = NULL,
  start = NULL, n_starts = 5L, seed = 1L, stable = FALSE,
  control = list(iter.max = 1000L, eval.max = 2000L, rel.tol = 1e-10),
  log_file = NULL, verbose = TRUE, rank_tolerance = 1e-7,
  boundary_tolerance = 1e-6, gradient_tolerance = 1e-3,
  feasibility = c("check", "penalty"), penalty_weight = 1e6)
}
\arguments{
\item{model}{A \code{gcmsem_model}.}
\item{moments}{Named numeric cumulants, a data frame with key/estimate columns, or a bundle from read_moments.}
\item{V}{Full sampling covariance of the estimated moments. May be singular PSD, with positive fitted-moment diagonal. Omitted for a bundle.}
\item{orders}{Observed cumulant orders to use in fitting. Does not change source parameters.}
\item{keys}{Specific observed moment keys to fit; overrides orders.}
\item{start}{All free parameters in original model units. NULL derives starts from model defaults and observed SNP-effect variances.}
\item{n_starts}{Number of reproducible nlminb starts.}
\item{seed}{Restart RNG seed; caller RNG state is preserved.}
\item{stable}{Require a returned solution with spectral radius of B less than one.}
\item{control}{nlminb optimizer control list.}
\item{log_file}{Optional timestamped text log.}
\item{verbose}{Write progress to the console.}
\item{rank_tolerance}{Relative singular-value rank threshold.}
\item{boundary_tolerance}{Boundary threshold in internally scaled units.}
\item{gradient_tolerance}{Maximum interior standardized gradient for inference.}
\item{feasibility}{Default unpenalized diagnostic checks, or an explicitly requested exact-gradient source-feasibility penalty sensitivity estimator.}
\item{penalty_weight}{Initial positive penalty strength, increased up to 1e12 until independent normalized source checks pass.}
}
\description{Fits DWLS using exact structured Jacobians, automatic parameter
scaling, multistart optimization and full-sampling-covariance sandwich uncertainty.}
\details{The objective uses only the diagonal of V. Sandwich uncertainty uses all
of V without an extra sample-size division. Model/parameter scale is preserved.
Selected keys subset both dimensions of V. Default unconstrained DWLS chooses the minimum converged raw Q regardless of source
admissibility and allows negative free residual variances. Classical and higher-order
Heywood cases emit prominent warnings. Point estimates and computable formal sandwich
SEs are retained; source_admissible, formal_se_available and inference_valid are
separate flags. Conditional Wald p-values are unavailable for flagged inference.
Nonconvergence, large gradients, rank deficiency or invalid sandwich variance can
make formal SEs unavailable. The optional penalty estimator uses exact-gradient
squared negative-eigenvalue penalties for necessary source moment conditions and
is explicitly labelled as a different estimator. Raw Q and penalized objective are
distinct. Numerical penalty enforcement does not guarantee exact feasibility or
global optimality. Q is uncalibrated: no WLSMV correction or chi-square
fit p-value is supplied. Inference is local and conditional on the model and solution.}
\value{A \code{gcmsem_fit} with coefficients, parameter table, sandwich covariance,
fitted/observed cumulants, residuals, discrepancy, per-start results and diagnostics.}
\examples{
m <- gcmsem_model('y ~ b*x', c('x','y'), gaussian_observed=FALSE)
p <- m$start; p['b'] <- 0.4
v <- implied_moments(m,p)
V <- diag(0.01,length(v)); dimnames(V) <- list(names(v),names(v))
fit <- gcmsem(m,v,V,n_starts=1,verbose=FALSE)
summary(fit)
}
