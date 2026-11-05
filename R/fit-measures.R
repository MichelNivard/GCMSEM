# Local quadratic-form diagnostics for the minimized DWLS discrepancy.
# J and residuals are in units of moment SEs; C is their full correlation.
.q_uncertainty <- function(residuals, J, C, available=TRUE, reason="") {
  d <- length(residuals); p <- ncol(J)
  ans <- data.frame(Q=sum(residuals^2), Q_se=NA_real_, Q_se_delta=NA_real_,
    Q_null_mean=NA_real_, Q_null_sd=NA_real_, Q_null_95=NA_real_,
    Q_reference_scale=NA_real_, Q_reference_df=NA_real_, Q_p_approx=NA_real_,
    residual_df=d-p, Q_status=if(available) "local_gaussian_approximation" else reason,
    stringsAsFactors=FALSE)
  if(!available)return(ans)
  P <- diag(d)
  if(p) {
    U <- svd(J, nu=p, nv=0L)$u
    P <- P-tcrossprod(U)
  }
  Cr <- P %*% C %*% P
  Cr <- (Cr+t(Cr))/2
  # Remove only eigenvalue roundoff accepted by input covariance validation.
  ee <- eigen(Cr, symmetric=TRUE)
  lambda <- pmax(ee$values, 0)
  lambda[lambda < max(1,max(lambda))*1e-12] <- 0
  Cr <- tcrossprod(sweep(ee$vectors,2L,sqrt(lambda),"*"))
  mean0 <- sum(lambda); var0 <- 2*sum(lambda^2)
  delta_var <- max(0,4*as.numeric(crossprod(residuals,Cr %*% residuals)))
  ans$Q_se_delta <- sqrt(delta_var)
  ans$Q_se <- sqrt(var0+delta_var)
  ans$Q_null_mean <- mean0; ans$Q_null_sd <- sqrt(var0)
  if(mean0>0 && var0>0) {
    scale <- var0/(2*mean0); df <- 2*mean0^2/var0
    ans$Q_reference_scale <- scale; ans$Q_reference_df <- df
    ans$Q_null_95 <- scale*stats::qchisq(.95,df=df)
    ans$Q_p_approx <- stats::pchisq(ans$Q/scale,df=df,lower.tail=FALSE)
  }
  ans
}

#' Extract model discrepancy and approximate uncertainty
#' @export
fit_measures <- function(fit) {
  if(!inherits(fit,"gcmsem_fit"))stop("fit must be a gcmsem_fit")
  if(is.null(fit$fit_measures))stop("Refit this older saved object to calculate Q uncertainty")
  fit$fit_measures
}
