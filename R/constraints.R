# Squared distance to the PSD cone, evaluated on necessary source moment matrices.
# The negative spectral part is differentiable even at repeated eigenvalues;
# its derivative is the Frobenius product with 2 U diag(min(lambda,0)) U'.
.constraint_setup <- function(model,initial,joint_max_nodes=12L) {
  tab<-.source_tables(model,initial)
  scales<-sqrt(pmax(abs(diag(tab$C)),1e-300))
  scales[diag(tab$C)==0]<-1
  list(scales=scales,joint=length(model$nodes)<=joint_max_nodes)
}
.source_penalty <- function(model,par,setup,derivative=TRUE) {
  n<-length(model$nodes);p<-length(par);sc<-setup$scales
  K<-vector('list',4L);for(q in 2:4)K[[q]]<-list()
  for(r in model$rows)if(r$kind=='cumulant') {
    den<-prod(sc[r$tuple]);val<-if(r$parameter>0)par[r$parameter] else r$value
    d<-numeric(p);if(derivative && r$parameter>0)d[r$parameter]<-1/den
    K[[r$order]][[paste(r$tuple,collapse=',')]]<-list(value=unname(val/den),d=d)
  }
  get<-function(t) {
    z<-K[[length(t)]][[paste(sort(t),collapse=',')]]
    if(is.null(z))list(value=0,d=numeric(p)) else z
  }
  plus<-function(a,b)list(value=a$value+b$value,d=a$d+b$d)
  mul<-function(a,b)list(value=a$value*b$value,d=a$d*b$value+a$value*b$d)
  raw<-function(t) {
    q<-length(t)
    if(q==0)return(list(value=1,d=numeric(p)))
    if(q==1)return(list(value=0,d=numeric(p)))
    z<-get(t)
    if(q==4)for(pair in list(c(1,2,3,4),c(1,3,2,4),c(1,4,2,3)))
      z<-plus(z,mul(get(t[pair[1:2]]),get(t[pair[3:4]])))
    z
  }
  matrix_parts<-function(monomials) {
    d<-length(monomials)
    if(derivative && d*d*p>5e7)stop("Source-constraint Jacobian workspace exceeds 50 million entries; reduce source model or use feasibility=check with explicit diagnostics")
    M<-matrix(0,d,d);D<-if(derivative)array(0,c(d,d,p)) else NULL
    for(i in seq_len(d))for(j in seq_len(i)) {
      z<-raw(c(monomials[[i]],monomials[[j]]))
      M[i,j]<-M[j,i]<-z$value
      if(derivative && p) {D[i,j,]<-z$d;D[j,i,]<-z$d}
    }
    list(M=M,D=D)
  }
  # Fixed normalization avoids any unit-dependent parameter floors and retains
  # an exact gradient. Row/column constants do not change PSD feasibility.
  if(setup$joint) {
    mono<-c(list(integer()),lapply(seq_len(n),function(i)i),.combinations(n,2))
    parts<-list(matrix_parts(mono))
  } else {
    parts<-list(matrix_parts(lapply(seq_len(n),function(i)i)))
    for(i in seq_len(n))parts[[length(parts)+1L]]<-matrix_parts(list(integer(),i,c(i,i)))
  }
  loss<-0;gradient<-numeric(p);minimum<-Inf
  for(part in parts) {
    e<-eigen(part$M,symmetric=TRUE);neg<-which(e$values<0)
    minimum<-min(minimum,min(e$values))
    if(!length(neg))next
    loss<-loss+sum(e$values[neg]^2)
    if(derivative && p) {
      U<-e$vectors[,neg,drop=FALSE]
      G<-sweep(U,2,2*e$values[neg],'*')%*%t(U)
      for(j in seq_len(p))gradient[j]<-gradient[j]+sum(G*part$D[,,j])
    }
  }
  list(value=loss,gradient=gradient,min_eigen=minimum,joint=setup$joint)
}
