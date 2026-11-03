\name{parameter_estimates}
\alias{parameter_estimates}
\alias{coef.gcmsem_fit}
\alias{vcov.gcmsem_fit}
\alias{fitted.gcmsem_fit}
\alias{print.gcmsem_fit}
\alias{summary.gcmsem_fit}
\title{Inspect a fitted genomic comoment SEM}
\usage{
parameter_estimates(fit)
\method{coef}{gcmsem_fit}(object, ...)
\method{vcov}{gcmsem_fit}(object, ...)
\method{fitted}{gcmsem_fit}(object, ...)
\method{print}{gcmsem_fit}(x, ...)
\method{summary}{gcmsem_fit}(object, ...)
}
\arguments{
\item{fit, object, x}{A fitted \code{gcmsem_fit}.}
\item{...}{Unused additional arguments.}
}
\description{Extract original-unit parameters and their conditional robust
uncertainty, model-implied cumulants, or display the numerical diagnostics.}
\value{parameter_estimates returns a data frame with parameter, estimate, se,
se_status, z and conditional Wald p. coef returns a named vector; vcov a matrix;
fitted a named vector. Print and summary return the fit invisibly. Inadmissible source moments retain computable formal SEs labelled formal_only_flagged.
Numerically unavailable SEs are NA. Flagged local inference has NA conditional p-values
and is explained in diagnostics.}
