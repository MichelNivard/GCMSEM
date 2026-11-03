# Reproducible native/R oracle and derivative check; run after R CMD INSTALL .
library(GCMSEM)
symbol <- getDLLRegisteredRoutines(getLoadedDLLs()[["GCMSEM"]])$.Call$C_gcmsem_contract
.combinations <- getFromNamespace(".combinations", "GCMSEM")
.permutations <- getFromNamespace(".permutations", "GCMSEM")
native <- function(A,dA,K,tuples,p,jacobian) .Call(symbol,A,dA,K,tuples,as.integer(p),jacobian)
reference <- function(A,dA,K,tuples,p,jacobian) {
 v <- numeric(length(tuples)); J <- matrix(0,length(tuples),p)
 for(m in seq_along(tuples)) {
  I<-tuples[[m]]; q<-length(I)
  for(r in K[[q]]) {
   co<-0; dc<-numeric(p)
   for(perm in .permutations(r$tuple)) {
    aa<-A[cbind(I,perm)]; co<-co+prod(aa)
    if(jacobian && p>0)for(k in seq_len(q))dc<-dc+prod(aa[-k])*dA[I[k],perm[k],]
   }
   v[m]<-v[m]+r$value*co
   if(jacobian) {J[m,]<-J[m,]+r$value*dc;if(r$parameter>0)J[m,r$parameter]<-J[m,r$parameter]+co}
  }
 }
 list(estimate=v,jacobian=if(jacobian)J else NULL)
}
set.seed(42496)
max_value <- max_jacobian <- max_difference <- 0
for(nobs in c(1,2,5,7))for(p in c(0L,1L,11L))for(zeros in c(FALSE,TRUE)) {
 ns<-nobs+2L; A<-matrix(rnorm(nobs*ns),nobs,ns); dA<-array(rnorm(nobs*ns*p),c(nobs,ns,p))
 if(zeros)A[A>0]<-0
 K<-vector('list',4)
 for(q in 2:4) {
  rows<-.combinations(ns,q); selected<-unique(round(seq(1,length(rows),length.out=min(25,length(rows)))))
  K[[q]]<-lapply(rows[selected],function(t)list(tuple=as.integer(t),value=sample(c(0,rnorm(2)),1),parameter=if(p)sample(0:p,1) else 0L))
 }
 tt<-unlist(lapply(2:4,function(q).combinations(nobs,q)),recursive=FALSE)
 for(jac in c(FALSE,TRUE)) {
  x<-native(A,if(jac)dA else NULL,K,tt,p,jac); y<-reference(A,dA,K,tt,p,jac)
  ve<-max(abs(x$estimate-y$estimate)); je<-if(jac && p)max(abs(x$jacobian-y$jacobian)) else 0
  max_value<-max(max_value,ve);max_jacobian<-max(max_jacobian,je)
  stopifnot(ve<1e-10,je<1e-9)
 }
 # Independent directional finite difference for 5-observed-variable cases.
 if(nobs==5 && p==11) {
  direction<-rnorm(p); eps<-1e-6; deltaA<-apply(dA,c(1,2),function(v)sum(v*direction))
  shifted<-function(sign) {
   kk<-lapply(K,function(rows)lapply(rows,function(r){if(r$parameter)r$value<-r$value+sign*eps*direction[r$parameter];r}))
   native(A+sign*eps*deltaA,NULL,kk,tt,p,FALSE)$estimate
  }
  delta<-(shifted(1)-shifted(-1))/(2*eps)
  analytic<-as.vector(native(A,dA,K,tt,p,TRUE)$jacobian%*%direction)
  relative<-max(abs(delta-analytic))/max(1,max(abs(analytic)))
  max_difference<-max(max_difference,relative); stopifnot(relative<1e-7)
 }
}
# Repeated index multiplicity: K1112 produces 4, K1122 produces 6.
K<-vector('list',4);K[[2]]<-K[[3]]<-list();K[[4]]<-list(list(tuple=c(1L,1L,1L,2L),value=1,parameter=0L))
stopifnot(native(matrix(1,1,2),NULL,K,list(rep(1L,4)),0L,FALSE)$estimate==4)
K[[4]][[1]]$tuple<-c(1L,1L,2L,2L)
stopifnot(native(matrix(1,1,2),NULL,K,list(rep(1L,4)),0L,FALSE)$estimate==6)
cat(sprintf('PASS: 48 random value/Jacobian comparisons, N=1,2,5,7; p=0,1,11; exact zero loading stress. Max value %.3g; max Jacobian %.3g; finite-difference relative error %.3g\n',max_value,max_jacobian,max_difference))
# Moderately sized benchmark of the same computational work, excluding model state.
nobs<-8L;ns<-10L;p<-20L;A<-matrix(rnorm(nobs*ns),nobs,ns);dA<-array(rnorm(nobs*ns*p),c(nobs,ns,p));K<-vector('list',4)
for(q in 2:4)K[[q]]<-lapply(1:ns,function(i)list(tuple=rep(as.integer(i),q),value=rnorm(1),parameter=sample(0:p,1)))
tt<-unlist(lapply(2:4,function(q).combinations(nobs,q)),recursive=FALSE)
a<-system.time(for(k in 1:3)native(A,dA,K,tt,p,TRUE))['elapsed']/3
b<-system.time(reference(A,dA,K,tt,p,TRUE))['elapsed']
cat(sprintf('Benchmark N=8 sources=10 parameters=20 moments=%d: native %.4fs R %.4fs speedup %.1fx (single environment)\n',length(tt),a,b,b/max(a,1e-9)))
