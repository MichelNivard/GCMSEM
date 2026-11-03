.validate_input <- function(moments,V,orders,keys) {
  if(inherits(moments,"gcmsem_moments")) {
    if(is.null(V))V<-moments$V
    moments<-moments$moments
  }
  if(is.data.frame(moments)) {
    if(!all(c("key","estimate")%in%names(moments)))stop("Moment table requires key and estimate")
    m<-setNames(moments$estimate,moments$key)
  } else m<-moments
  if(!is.numeric(m) || is.null(names(m)) || anyDuplicated(names(m)) || any(!is.finite(m)))stop("moments must be finite numeric cumulants with unique keys")
  order<-lengths(strsplit(names(m),"|",fixed=TRUE))
  if(is.null(keys))keys<-names(m)[order%in%orders]
  if(!length(keys) || any(!keys%in%names(m)) || anyDuplicated(keys))stop("No valid selected moment keys")
  if(is.null(V))stop("Supply the sampling covariance V (not covariance of SNP effects)")
  V<-as.matrix(V)
  if(!is.numeric(V) || nrow(V)!=ncol(V) || any(!is.finite(V)))stop("V must be a finite square numeric matrix")
  if(!is.null(rownames(V)) || !is.null(colnames(V))) {
    if(is.null(rownames(V)) || is.null(colnames(V)) || anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)) ||
       !all(keys%in%rownames(V)) || !all(keys%in%colnames(V)))stop("V row and column keys must include all fitted moments")
    V<-V[keys,keys,drop=FALSE]
  } else {
    if(nrow(V)!=length(m))stop("Unnamed V must follow all input moments in their exact order")
    i<-match(keys,names(m));V<-V[i,i,drop=FALSE]
    dimnames(V)<-list(keys,keys)
  }
  if(any(diag(V)<=0))stop("Every fitted moment requires positive finite sampling variance")
  sn<-sqrt(diag(V)); C<-V/outer(sn,sn)
  if(max(abs(C-t(C)))>1e-8)stop("V must be symmetric")
  C<-(C+t(C))/2; ev<-eigen(C,symmetric=TRUE,only.values=TRUE)$values
  if(min(ev) < -1e-7*max(ev))stop("V is not positive semidefinite")
  list(m=m[keys],V=(V+t(V))/2,se=sn,keys=keys,corr=C,V_rank=sum(ev>max(ev)*1e-9),V_eigen=ev)
}

.parameter_scales <- function(model,m) {
  scales<-setNames(rep(1,length(model$nodes)),model$nodes)
  for(v in model$observed) {
    key<-paste(v,v,sep="|")
    if(key%in%names(m) && m[key]!=0)scales[v]<-sqrt(abs(m[key]))
  }
  p<-length(model$start); candidates<-vector("list",p)
  for(r in model$rows)if(r$parameter>0L) {
    ss<-scales[model$nodes[r$tuple]]
    val<-if(r$kind=="path")ss[1]/ss[2] else prod(ss)
    candidates[[r$parameter]]<-c(candidates[[r$parameter]],val)
  }
  setNames(vapply(candidates,function(x)exp(mean(log(pmax(x,1e-250)))),numeric(1)),names(model$start))
}

.log_writer <- function(path,verbose) {
  if(!is.null(path)) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    cat("",file=path)
  }
  function(...) {
    line<-paste(format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),paste(...,collapse=" "))
    if(verbose)cat(line,"\n")
    if(!is.null(path))cat(line,"\n",file=path,append=TRUE)
  }
}

#' Fit a genomic comoment structural equation model
#' @export
gcmsem <- function(model,moments,V=NULL,orders=2:4,keys=NULL,start=NULL,
                   n_starts=5L,seed=1L,stable=FALSE,
                   control=list(iter.max=1000L,eval.max=2000L,rel.tol=1e-10),
                   log_file=NULL,verbose=TRUE,rank_tolerance=1e-7,
                   boundary_tolerance=1e-6,gradient_tolerance=1e-3,
                   feasibility=c("check","penalty"),penalty_weight=1e6) {
  if(!inherits(model,"gcmsem_model"))stop("Build the model using gcmsem_model()")
  feasibility<-match.arg(feasibility)
  if(!is.numeric(penalty_weight)||length(penalty_weight)!=1||!is.finite(penalty_weight)||penalty_weight<=0)stop("penalty_weight must be positive and finite")
  if(!length(orders) || any(!orders%in%2:4) || anyDuplicated(orders))stop("orders must be selected from 2:4")
  dat<-.validate_input(moments,V,orders,keys); .parse_keys(dat$keys,model$observed)
  p<-length(model$start); d<-length(dat$m)
  if(!is.numeric(n_starts)||length(n_starts)!=1||n_starts<1||n_starts!=as.integer(n_starts))stop("n_starts must be a positive integer")
  log<-.log_writer(log_file,verbose)
  log("GCMSEM 0.1.0 | DWLS with full-V sandwich | R",getRversion())
  log("MODEL",gsub("[\r\n]+"," ; ",model$syntax))
  log("SETTINGS seed",seed,"starts",n_starts,"stable",stable,"feasibility",feasibility)
  log("Observed",length(model$observed),"latent",length(model$latent),"moments",d,"parameters",p,"V rank",dat$V_rank)
  log("Moment orders",paste(sort(unique(lengths(strsplit(dat$keys,"|",fixed=TRUE)))),collapse=","),"; exact structured first derivatives")
  scale<-.parameter_scales(model,dat$m)
  initial<-if(is.null(start))model$start*scale else .resolve_par(model,start)
  effective_lower<-if(feasibility=="check")rep(-Inf,p) else model$lower
  if(any(initial<effective_lower))stop("start violates parameter bounds")
  lo<-effective_lower/scale; hi<-model$upper/scale
  constraint_setup<-.constraint_setup(model,initial)
  weight<-if(feasibility=="penalty")penalty_weight else 0
  log("Source feasibility",feasibility,"; full joint moment penalty",constraint_setup$joint)
  cache_z<-NULL;cache<-NULL
  evaluate<-function(z) {
    if(identical(z,cache_z))return(cache)
    val<-tryCatch({
      par<-setNames(z*scale,names(scale))
      im<-implied_moments(model,par,jacobian=TRUE,keys=dat$keys)
      e<-(im$estimate-dat$m)/dat$se
      J<-sweep(sweep(im$jacobian,1,dat$se,"/"),2,scale,"*")
      pen<-if(weight>0).source_penalty(model,par,constraint_setup) else list(value=0,gradient=rep(0,p))
      ans<-list(objective=sum(e*e)+weight*pen$value,
        dwls=sum(e*e),penalty=pen$value,
        gradient=2*as.vector(crossprod(J,e))+weight*pen$gradient*scale,
        dwls_gradient=2*as.vector(crossprod(J,e)),implied=im$estimate,J=J)
      if(!is.finite(ans$objective)||any(!is.finite(ans$gradient)))stop("Non-finite objective")
      ans
    },error=function(e)list(objective=1e100,gradient=rep(0,p),error=conditionMessage(e)))
    cache_z<<-z;cache<<-val;val
  }
  fun<-function(z)evaluate(z)$objective
  grad<-function(z)evaluate(z)$gradient
  # Preserve caller RNG state, so optimizer restarts never change simulations.
  old_seed<-if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))get(".Random.seed",.GlobalEnv) else NULL
  on.exit(if(is.null(old_seed)) {if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)} else assign(".Random.seed",old_seed,.GlobalEnv),add=TRUE)
  set.seed(seed)
  candidates<-list(); trace<-list()
  for(k in seq_len(n_starts)) {
    z<-initial/scale
    if(k>1L && p) {
      z<-z+stats::rnorm(p,sd=0.25)
      z<-pmax(lo,pmin(hi,z))
    }
    weight<-if(feasibility=="penalty")penalty_weight else 0
    cache_z<-NULL;cache<-NULL
    opt<-if(p)tryCatch(stats::nlminb(z,objective=fun,gradient=grad,lower=lo,upper=hi,control=control),
                           error=function(e)list(par=z,objective=Inf,convergence=99L,message=conditionMessage(e),iterations=0L)) else
      list(par=numeric(),objective=fun(numeric()),convergence=0L,message="No free parameters",iterations=0L)
    par<-setNames(opt$par*scale,names(scale))
    src<-tryCatch(source_diagnostics(model,par),error=function(e)list(necessary_conditions_pass=FALSE,error=conditionMessage(e)))
    # Quadratic penalties are refined until the independent normalized PSD
    # diagnostics pass or a documented maximum penalty is reached.
    while(feasibility=="penalty" && p && !isTRUE(src$necessary_conditions_pass) && weight<1e12 && is.finite(opt$objective)) {
      weight<-min(weight*1e3,1e12);cache_z<-NULL;cache<-NULL
      refined_control<-utils::modifyList(control,list(x.tol=1e-14,rel.tol=1e-13))
      opt<-tryCatch(stats::nlminb(opt$par,objective=fun,gradient=grad,lower=lo,upper=hi,control=refined_control),
        error=function(e)list(par=opt$par,objective=Inf,convergence=99L,message=conditionMessage(e),iterations=0L))
      par<-setNames(opt$par*scale,names(scale))
      src<-tryCatch(source_diagnostics(model,par),error=function(e)list(necessary_conditions_pass=FALSE,error=conditionMessage(e)))
    }
    result<-evaluate(opt$par)
    dwls<-if(is.null(result$dwls))Inf else result$dwls
    pg<-result$gradient
    if(p) {
      pg[is.finite(lo)&opt$par<=lo+boundary_tolerance&pg>0]<-0
      pg[is.finite(hi)&opt$par>=hi-boundary_tolerance&pg<0]<-0
    }
    stationarity<-if(p)max(abs(pg)) else 0
    st<-tryCatch(.model_state(model,par),error=function(e)NULL)
    radius<-if(is.null(st))Inf else max(Mod(eigen(st$B,only.values=TRUE)$values))
    accepted<-is.finite(opt$objective)&&opt$objective<1e99&&(feasibility=="check"||src$necessary_conditions_pass)&&(!stable||radius<1-boundary_tolerance)
    trace[[k]]<-data.frame(start=k,objective=dwls,penalized_objective=opt$objective,penalty_weight=weight,convergence=opt$convergence,projected_gradient=stationarity,
                          moment_conditions=src$necessary_conditions_pass,spectral_radius=radius,
                          accepted=accepted,message=opt$message,stringsAsFactors=FALSE)
    candidates[[k]]<-list(opt=opt,par=par,source=src,radius=radius,accepted=accepted,dwls=dwls,weight=weight,stationarity=stationarity)
    log("Start",k,"Q",format(dwls,digits=8),"convergence",opt$convergence,"moment checks",src$necessary_conditions_pass,"spectral radius",format(radius,digits=4))
  }
  valid<-which(vapply(candidates,function(z)z$accepted && z$opt$convergence==0L,logical(1)))
  if(!length(valid))valid<-which(vapply(candidates,function(z)z$accepted,logical(1)))
  if(!length(valid))valid<-seq_along(candidates)
  best<-candidates[[valid[which.min(vapply(candidates[valid],function(z)z$dwls,numeric(1)))]]]
  weight<-best$weight;cache_z<-NULL;cache<-NULL
  par<-best$par; final<-evaluate(best$opt$par)
  if(is.null(final$J))stop("All starts failed: ",if(!is.null(final$error))final$error else paste(unique(vapply(candidates,function(z)z$opt$message,"")),collapse="; "))
  J<-final$J
  singular<-if(p)svd(J,nu=0,nv=0)$d else numeric()
  rank<-if(length(singular))sum(singular>max(singular)*rank_tolerance) else 0L
  rank<-as.integer(rank)
  gradmax<-if(p)max(abs(final$gradient)) else 0
  raw_gradmax<-if(p)max(abs(final$dwls_gradient)) else 0
  active<-if(p)which(is.finite(lo)&abs(best$opt$par-lo)<boundary_tolerance | is.finite(hi)&abs(best$opt$par-hi)<boundary_tolerance) else integer()
  # Moment-boundary checks are scale invariant; fixed zero nuisance components are not free boundaries.
  src_boundary<-isTRUE(abs(best$source$covariance_min_eigen)<boundary_tolerance) ||
    (!is.null(best$source$pearson_margin)&&any(abs(best$source$pearson_margin)<boundary_tolerance)) ||
    (isTRUE(best$source$joint_moment_checked)&&isTRUE(abs(best$source$joint_moment_min_eigen)<boundary_tolerance))
  reasons<-character()
  if(best$opt$convergence!=0L)reasons<-c(reasons,"optimizer did not converge")
  if(gradmax>gradient_tolerance && !length(active))reasons<-c(reasons,"gradient remains large")
  if(rank<p)reasons<-c(reasons,"model Jacobian is rank deficient")
  numerical_reasons<-reasons
  if(length(active)||src_boundary)reasons<-c(reasons,"parameter or source moment boundary")
  if(!isTRUE(best$source$necessary_conditions_pass))reasons<-c(reasons,"source moment feasibility checks failed")
  if(stable && best$radius>=1-boundary_tolerance)reasons<-c(reasons,"requested stable branch not found")
  Cpar<-matrix(NA_real_,p,p,dimnames=list(names(par),names(par)))
  sandwich_rank<-NA_integer_
  if(rank==p && p) {
    # SVD inversion of J'J avoids squaring its condition number numerically.
    s<-svd(J); bread<-s$v%*%diag(1/(s$d^2),nrow=p)%*%t(s$v)
    meat<-crossprod(J,dat$corr%*%J)
    cz<-bread%*%meat%*%bread;cz<-(cz+t(cz))/2
    eig<-eigen(cz,symmetric=TRUE,only.values=TRUE)$values
    sandwich_rank<-sum(eig>max(abs(eig))*rank_tolerance)
    if(any(diag(cz) < -1e-10*max(abs(diag(cz))))) {
      reasons<-c(reasons,"negative sandwich variance")
      numerical_reasons<-c(numerical_reasons,"negative sandwich variance")
    }
    if(!length(numerical_reasons))Cpar<-cz*outer(scale,scale)
  }
  se<-if(p)sqrt(pmax(diag(Cpar),0)) else numeric()
  formal_se_available<-!length(numerical_reasons)
  parameters<-data.frame(parameter=names(par),estimate=unname(par),se=unname(se),
    se_status=rep(if(!formal_se_available)"unavailable_numerical" else if(length(reasons))"formal_only_flagged" else "regular_local",p),
    z=unname(par/se),p=if(length(reasons))rep(NA_real_,p) else 2*stats::pnorm(-abs(unname(par/se))),stringsAsFactors=FALSE)
  diagnostics<-list(converged=best$opt$convergence==0L,jacobian_rank=rank,n_parameters=p,
    n_moments=d,residual_df=d-p,V_rank=dat$V_rank,sandwich_rank=sandwich_rank,
    gradient_max=gradmax,dwls_gradient_max=raw_gradmax,active_parameters=names(par)[active],spectral_radius=best$radius,
    stable=best$radius<1,stability_required=stable,source=best$source,
    inference_valid=!length(reasons),inference_reasons=unique(reasons),
    formal_se_available=formal_se_available,formal_se_reasons=unique(numerical_reasons),
    source_admissible=isTRUE(best$source$necessary_conditions_pass),
    local_identification_only=TRUE,fit_statistic_calibrated=FALSE,
    feasibility_method=feasibility,penalty_weight=best$weight,
    condition_number=if(p&&rank==p)max(singular)/min(singular) else Inf)
  fit<-list(model=model,coefficients=par,parameters=parameters,vcov=Cpar,
    fitted=final$implied,observed=dat$m,V=dat$V,residuals=dat$m-final$implied,
    objective=final$dwls,penalized_objective=best$opt$objective,
    constraint_penalty=final$penalty,estimator=if(feasibility=="check")"DWLS" else "DWLS with source-feasibility penalty",diagnostics=diagnostics,
    starts=do.call(rbind,trace),parameter_scale=scale,call=match.call())
  class(fit)<-"gcmsem_fit"
  log("Selected Q",format(fit$objective,digits=9),"Jacobian rank",paste0(rank,"/",p),"inference",diagnostics$inference_valid)
  for(i in seq_len(nrow(parameters)))log("PARAMETER",parameters$parameter[i],"estimate",format(parameters$estimate[i],digits=10),"formal_se",format(parameters$se[i],digits=10),"status",parameters$se_status[i])
  if(length(numerical_reasons))log("Formal SEs unavailable:",paste(unique(numerical_reasons),collapse="; "))
  if(length(reasons))log("Inference flagged; point estimates and computable formal SEs retained:",paste(unique(reasons),collapse="; "))
  if(!isTRUE(best$source$necessary_conditions_pass)) {
    flags<-character()
    if(length(best$source$negative_variances))flags<-c(flags,paste0("CLASSICAL HEYWOOD CASE: negative source variance for ",paste(best$source$negative_variances,collapse=", ")))
    if(isTRUE(best$source$covariance_min_eigen < -1e-8))flags<-c(flags,paste0("CLASSICAL HEYWOOD CASE: source covariance is not PSD (min eigenvalue ",format(best$source$covariance_min_eigen,digits=5),")"))
    bad_pearson<-best$source$pearson_margin[is.finite(best$source$pearson_margin)&best$source$pearson_margin < -1e-8]
    if(length(bad_pearson))flags<-c(flags,paste0("HIGHER-ORDER HEYWOOD CASE: Pearson bound fails for ",paste(names(bad_pearson),collapse=", ")))
    if(isTRUE(best$source$joint_moment_min_eigen < -1e-8))flags<-c(flags,paste0("HIGHER-ORDER HEYWOOD CASE: joint source moment matrix is not PSD (min eigenvalue ",format(best$source$joint_moment_min_eigen,digits=5),")"))
    if(!length(flags))flags<-"SOURCE MOMENT HEYWOOD CASE: necessary source moment conditions fail"
    msg<-paste(c(flags,"All point estimates and computable formal sandwich SEs are retained. These diagnostics can reflect sampling variation or model problems; Gaussian/Wald inference is flagged."),collapse="; ")
    log(msg);warning(msg,call.=FALSE)
  }
  log("Q is an uncalibrated discrepancy; no chi-square/WLSMV p-value. SEs are conditional on the specified model and selected solution.")
  fit
}
