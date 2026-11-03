# Symmetric tuple indexing is sparse: no p^4 array is materialized.
.combinations <- function(n,q,start=1L) {
  if(q==0L)return(list(integer()))
  out<-list()
  for(i in seq.int(start,n)) for(t in .combinations(n,q-1L,i))out[[length(out)+1L]]<-c(i,t)
  out
}
.permutations <- function(x) {
  if(length(x)==1L)return(list(x))
  out<-list()
  for(v in unique(x))for(z in .permutations(x[-match(v,x)]))out[[length(out)+1L]]<-c(v,z)
  out
}

#' Canonical cumulant keys
#' @export
moment_keys <- function(traits,orders=2:4) {
  if(!length(traits) || anyDuplicated(traits) || any(!.name_ok(traits)))stop("Invalid traits")
  if(!length(orders) || any(!orders %in% 1:4) || anyDuplicated(orders))stop("orders must be unique values from 1:4")
  total<-sum(vapply(orders,function(q)choose(length(traits)+q-1,q),numeric(1)))
  if(!is.finite(total)||total>1e6)stop("More than one million symmetric moment keys requested; use selected keys or fewer orders/traits")
  unlist(lapply(sort(orders),function(q)vapply(.combinations(length(traits),q),function(t)paste(traits[t],collapse="|"),"")),use.names=FALSE)
}

.model_state <- function(model,par,derivatives=FALSE) {
  par<-.resolve_par(model,par); n<-length(model$nodes); p<-length(par)
  B<-matrix(0,n,n); K<-vector("list",4L); dB<-list()
  for(q in 2:4)K[[q]]<-list()
  for(r in model$rows) {
    value<-if(r$parameter>0L)par[r$parameter] else r$value
    if(r$kind=="path")B[r$tuple[1],r$tuple[2]]<-value else {
      r$value<-unname(value); K[[r$order]][[paste(r$tuple,collapse=",")]]<-r
    }
  }
  T<-tryCatch(solve(diag(n)-B),error=function(e)NULL)
  if(is.null(T) || any(!is.finite(T)))stop("Structural system I-B is singular")
  A<-T[seq_along(model$observed),,drop=FALSE]
  dA<-NULL
  if(derivatives) {
    if(length(A)*p>5e7)stop("Jacobian workspace exceeds 50 million entries; reduce model or select moments")
    dA<-array(0,c(dim(A),p))
    for(r in model$rows)if(r$kind=="path" && r$parameter>0L)
      dA[,,r$parameter]<-dA[,,r$parameter]+outer(A[,r$tuple[1]],T[r$tuple[2],])
  }
  list(B=B,T=T,A=A,dA=dA,K=K,par=par)
}

.parse_keys <- function(keys,observed) {
  if(!is.character(keys) || !length(keys) || anyDuplicated(keys))stop("keys must be unique moment keys")
  tuples<-strsplit(keys,"|",fixed=TRUE)
  lapply(tuples,function(t) {
    i<-match(t,observed)
    if(anyNA(i) || !length(t)%in%2:4 || is.unsorted(i))stop("Keys must contain 2-4 observed names in declared trait order: ",paste(t,collapse="|"))
    i
  })
}

#' Model-implied observed cumulants and exact Jacobian
#' @export
implied_moments <- function(model,par=NULL,orders=2:4,jacobian=FALSE,keys=NULL,backend=c("auto","R","native")) {
  if(!inherits(model,"gcmsem_model"))stop("model must be a gcmsem_model")
  if(is.null(keys)) {
    if(any(!orders%in%2:4))stop("GCMSEM fits centered cumulants of orders 2:4")
    keys<-moment_keys(model$observed,orders)
  }
  tuples<-.parse_keys(keys,model$observed)
  backend<-match.arg(backend)
  native_available<-exists("C_gcmsem_contract",inherits=TRUE)
  if(backend=="native" && !native_available)stop("Native kernel unavailable; install the package with a C++ compiler or use backend=R")
  if(jacobian && length(keys)*length(model$start)>5e7)stop("Requested moment Jacobian exceeds 50 million entries; select moment keys or simplify the model")
  s<-.model_state(model,par,jacobian); A<-s$A; p<-length(s$par)
  if(backend!="R" && native_available) {
    ans<-.Call(C_gcmsem_contract,A,s$dA,s$K,tuples,as.integer(p),as.logical(jacobian))
    names(ans$estimate)<-keys
    if(jacobian) {dimnames(ans$jacobian)<-list(keys,names(s$par));return(ans)}
    return(ans$estimate)
  }
  value<-setNames(numeric(length(keys)),keys)
  J<-matrix(0,length(keys),p,dimnames=list(keys,names(s$par)))
  for(m in seq_along(keys)) {
    I<-tuples[[m]]; q<-length(I)
    for(r in s$K[[q]]) {
      if(r$value==0 && (!jacobian || r$parameter==0L))next
      co<-0; dc<-numeric(p)
      for(perm in .permutations(r$tuple)) {
        aa<-A[cbind(I,perm)]; co<-co+prod(aa)
        if(jacobian && p>0L)for(k in seq_len(q)) {
          # Product rule remains correct when any loading is exactly zero.
          dc<-dc+prod(aa[-k])*s$dA[I[k],perm[k],]
        }
      }
      value[m]<-value[m]+r$value*co
      if(jacobian) {
        J[m,]<-J[m,]+r$value*dc
        if(r$parameter>0L)J[m,r$parameter]<-J[m,r$parameter]+co
      }
    }
  }
  if(jacobian)list(estimate=value,jacobian=J) else value
}

.source_tables <- function(model,par) {
  n<-length(model$nodes); K<-vector("list",4)
  for(q in 2:4)K[[q]]<-list()
  for(r in model$rows)if(r$kind=="cumulant")K[[r$order]][[paste(r$tuple,collapse=",")]]<-if(r$parameter>0)unname(par[r$parameter]) else r$value
  at<-function(t) {v<-K[[length(t)]][[paste(sort(t),collapse=",")]];if(is.null(v))0 else v}
  C<-outer(seq_len(n),seq_len(n),Vectorize(function(i,j)at(c(i,j))))
  list(K=K,at=at,C=C)
}

#' Necessary source moment feasibility diagnostics
#' @export
source_diagnostics <- function(model,par=NULL,joint_max_nodes=12L,tolerance=1e-8) {
  par<-.resolve_par(model,par); st<-.source_tables(model,par); C<-st$C; at<-st$at; n<-nrow(C)
  variances<-diag(C); scale<-sqrt(pmax(abs(variances),1e-300))
  active<-which(variances!=0)
  zero<-which(variances==0)
  zero_consistent<-TRUE
  if(length(zero))for(q in 2:4)for(key in names(st$K[[q]])) {
    ids<-as.integer(strsplit(key,",",fixed=TRUE)[[1]])
    if(any(ids %in% zero) && st$K[[q]][[key]]!=0)zero_consistent<-FALSE
  }
  cn<-C[active,active,drop=FALSE]/outer(scale[active],scale[active])
  ce<-if(length(active))eigen(cn,symmetric=TRUE,only.values=TRUE)$values else 1
  pearson<-vapply(seq_len(n),function(i) {
    if(variances[i]<=0)return(if(variances[i]==0 && at(rep(i,3))==0 && at(rep(i,4))==0)Inf else -Inf)
    at(rep(i,4))/variances[i]^2+2-at(rep(i,3))^2/variances[i]^3
  },numeric(1))
  names(pearson)<-model$nodes
  joint_min<-NA_real_; tested<-length(active)<=joint_max_nodes
  if(tested) {
    monomials<-c(list(integer()),lapply(active,function(i)i),
                  if(length(active))lapply(.combinations(length(active),2),function(t)active[t]) else list())
    get_raw<-function(t) {
      q<-length(t)
      val<-switch(as.character(q),"0"=1,"1"=0,"2"=at(t),"3"=at(t),
        "4"=at(t)+at(t[c(1,2)])*at(t[c(3,4)])+at(t[c(1,3)])*at(t[c(2,4)])+at(t[c(1,4)])*at(t[c(2,3)]))
      if(q)val/prod(scale[t]) else val
    }
    M<-outer(seq_along(monomials),seq_along(monomials),Vectorize(function(i,j)get_raw(c(monomials[[i]],monomials[[j]]))))
    # Normalize degree-two monomials too: eigen tolerance must be scale invariant.
    d<-sqrt(pmax(abs(diag(M)),1e-300)); M<-M/outer(d,d)
    joint_min<-min(eigen(M,symmetric=TRUE,only.values=TRUE)$values)
  }
  feasible<-all(variances>=0) && zero_consistent && min(ce)>=-tolerance && all(pearson>=-tolerance) && (!tested || joint_min>=-tolerance)
  list(necessary_conditions_pass=feasible,covariance_min_eigen=min(ce),
       source_variances=setNames(variances,model$nodes),negative_variances=model$nodes[variances<0],
       deterministic_nodes=model$nodes[zero],deterministic_zero_consistent=zero_consistent,
       pearson_margin=pearson,joint_moment_min_eigen=joint_min,joint_moment_checked=tested,
       note="PSD degree-two moment matrix and marginal Pearson bounds are necessary, not sufficient, for a joint distribution.")
}
