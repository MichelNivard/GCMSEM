# Independent Gaussian sampling and refitting checks for Q diagnostics.
local({
  library(GCMSEM)
  checks <- 0L
  check <- function(value,label) {
    checks <<- checks+1L
    if(!isTRUE(value))stop("Q uncertainty check failed: ",label)
  }
  near <- function(a,b,tol=1e-8)isTRUE(all.equal(unname(a),unname(b),tolerance=tol))
  calc <- GCMSEM:::.q_uncertainty

  # In an independent Gaussian linear problem, residual Q is chi-square(d-p).
  J <- cbind(c(1,1,1,1,1),c(-2,-1,0,1,2))
  zero <- calc(rep(0,5),J,diag(5))
  check(near(zero$Q_null_mean,3) && near(zero$Q_null_sd,sqrt(6)),"chi-square reduction")
  check(near(zero$Q_null_95,qchisq(.95,3)) && near(zero$Q_reference_df,3),"reference reduction")
  check(zero$Q_se_delta==0 && near(zero$Q_se,sqrt(6)),"zero misfit retains quadratic noise")

  # Monte Carlo independently fits each sample by ordinary least squares.
  C <- outer(1:5,1:5,function(i,j).6^abs(i-j))
  e <- residuals(lm(c(1,-1,2,-2,1) ~ J-1))
  expected <- calc(e,J,C)
  set.seed(9017)
  draws <- matrix(rnorm(5*80000),5)
  samples <- t(chol(C)) %*% draws
  sampled_q <- function(Y)colSums(qr.resid(qr(J),Y)^2)
  null_q <- sampled_q(samples)
  misfit_q <- sampled_q(samples+e)
  check(abs(mean(null_q)/expected$Q_null_mean-1)<.02,"correlated null mean by refitting")
  check(abs(sd(null_q)/expected$Q_null_sd-1)<.02,"correlated null SD by refitting")
  check(abs(sd(misfit_q)/expected$Q_se-1)<.02,"noncentral SD by refitting with known misfit")
  check(near(calc(e,J %*% matrix(c(3,1,0,2),2),C)$Q_se,expected$Q_se),"parameter-basis invariance")
  singular <- matrix(1,5,5)
  check(calc(e,J,singular)$Q_null_sd<1e-7,"singular V entirely in fitted tangent")
  saturated <- calc(rep(0,5),diag(5),C)
  check(saturated$Q_se==0 && is.na(saturated$Q_p_approx),"saturated space has no residual test")
  unavailable <- calc(e,J,C,FALSE,"unavailable_numerical")
  check(is.na(unavailable$Q_se) && unavailable$Q>0,"unavailable diagnostics retain Q")

  model <- gcmsem_model("x ~~ 1*x","x")
  m <- c("x|x"=1.18); V <- matrix(.01,1,1,dimnames=list(names(m),names(m)))
  fit <- gcmsem(model,m,V,orders=2,n_starts=1)
  q <- fit_measures(fit)
  check(near(q$Q,3.24) && near(q$Q_se,sqrt(2+4*3.24)),"fixed model integration")
  check(near(q$Q_p_approx,pchisq(3.24,1,lower.tail=FALSE)),"one-dimensional reference exact")
  check(identical(summary(fit)$fit_measures,q) && !fit$diagnostics$fit_statistic_calibrated,"summary retains approximation status")
  penalized <- gcmsem(model,m,V,orders=2,n_starts=1,feasibility="penalty")
  check(is.na(fit_measures(penalized)$Q_se),"penalized fits do not use unpenalized tangent")
  cat(checks,"Q uncertainty checks passed\n")
})
