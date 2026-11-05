\name{gcmsem_methods}
\alias{print.gcmsem_fit}
\alias{summary.gcmsem_fit}
\alias{print.summary_gcmsem_fit}
\alias{print.gcmsem_moments}
\alias{summary.gcmsem_moments}
\alias{print.summary_gcmsem_moments}
\alias{as.data.frame.gcmsem_fit}
\alias{as.data.frame.gcmsem_moments}
\alias{coef.gcmsem_moments}
\alias{vcov.gcmsem_moments}
\title{Inspect GCMSEM fits and moments in R}
\usage{
\method{print}{gcmsem_fit}(x, digits = 4L, ...)
\method{summary}{gcmsem_fit}(object, digits = 4L, ...)
\method{print}{summary_gcmsem_fit}(x, digits = x$digits, ...)
\method{print}{gcmsem_moments}(x, digits = 4L, max.traits = 6L, ...)
\method{summary}{gcmsem_moments}(object, digits = 4L, ...)
\method{print}{summary_gcmsem_moments}(x, digits = x$digits, max.traits = 6L, ...)
\method{as.data.frame}{gcmsem_fit}(x, row.names = NULL, optional = FALSE, ...)
\method{as.data.frame}{gcmsem_moments}(x, row.names = NULL, optional = FALSE, ...)
\method{coef}{gcmsem_moments}(object, ...)
\method{vcov}{gcmsem_moments}(object, ...)
}
\arguments{
\item{x,object}{A GCMSEM fit, moments object, or summary of the matching class.}
\item{digits}{Number of displayed digits, from 1 to 15. Stored estimates retain full precision.}
\item{max.traits}{Maximum number of traits shown in console matrix previews. Use \code{Inf} to print every trait. All matrices in the R object remain complete.}
\item{row.names}{Optional row names for a returned data frame.}
\item{optional}{Included for compatibility with \code{as.data.frame}.}
\item{...}{Additional generic arguments.}
}
\description{Provides compact console displays, structured summary objects and ordinary R extraction methods.}
\details{Assigning \code{s <- summary(fit)} returns a summary object without printing.
Typing \code{s} or \code{summary(fit)} at the R prompt displays its fit diagnostics and
parameter table. Moment summaries contain covariance and correlation matrices with
SEs. \code{as.data.frame(fit)} returns the parameter table;
\code{as.data.frame(moments)} returns the cumulant table.
\code{coef(moments)} returns named cumulant estimates and \code{vcov(moments)}
returns their full sampling covariance.
Source-moment Heywood cases remain prominent in fit displays, with their point estimates
and available formal sandwich SEs retained. The numeric matrices and data frames in
these objects can be passed directly to other R functions.}
\value{\code{summary()} returns a \code{summary_gcmsem_fit} or
\code{summary_gcmsem_moments} list. Print methods invisibly return their input.
Data-frame methods return ordinary data frames. Coefficient and covariance methods
return a named numeric vector and a named matrix, respectively.}
