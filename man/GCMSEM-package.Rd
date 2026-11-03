\name{GCMSEM-package}
\alias{GCMSEM-package}
\alias{GCMSEM}
\docType{package}
\title{Genomic Comoment Structural Equation Models}
\description{Standalone linear structural equation models for genomic second,
third and fourth cumulants, with exact derivatives and robust DWLS inference.}
\details{The source equation is \eqn{x = B x + u}; source cumulants are
propagated by \eqn{(I-B)^{-1}}. Fourth-order parameters are excess cumulants,
not fourth central moments. The Python companion prepares LD-corrected SNP-effect
cumulants and their full sampling covariance. See \code{gcmsem_model},
\code{gcmsem} and the repository documentation for assumptions and limits.}
\author{Michel Nivard}
\seealso{\code{\link{gcmsem_model}}, \code{\link{gcmsem}}}
