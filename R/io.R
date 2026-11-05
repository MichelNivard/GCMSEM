#' Read a Python GCMSEM moment bundle
#' @export
read_moments <- function(prefix) {
  if(!is.character(prefix) || length(prefix)!=1L || is.na(prefix) || !nzchar(prefix))
    stop("prefix must be one nonempty path")
  path<-paste0(prefix,".moments.tsv");vp<-paste0(prefix,".sampling_cov.tsv")
  m<-utils::read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
  if(!all(c("key","order","estimate","se")%in%names(m)))stop("Invalid moments.tsv schema")
  if(anyDuplicated(m$key) || any(!is.finite(m$estimate)) || any(!is.finite(m$se)) || any(m$se<0))stop("Invalid moment values")
  if(any(m$order!=lengths(strsplit(m$key,"|",fixed=TRUE))))stop("Moment order/key mismatch")
  v<-utils::read.delim(vp,check.names=FALSE,stringsAsFactors=FALSE)
  if(names(v)[1]!="key" || anyDuplicated(v$key))stop("Invalid sampling covariance keys")
  V<-as.matrix(v[,-1,drop=FALSE]);rownames(V)<-v$key
  if(!is.numeric(V)||any(!is.finite(V))||!setequal(rownames(V),m$key)||!setequal(colnames(V),m$key))stop("Sampling covariance and moment keys disagree")
  V<-V[m$key,m$key,drop=FALSE]
  if(any(abs(diag(V)-m$se^2)>1e-7*pmax(abs(diag(V)),m$se^2,1e-300)))stop("Moment SE and sampling covariance diagonal disagree")
  if(max(abs(V-t(V)))>1e-8*max(abs(V),1e-300))stop("Sampling covariance must be symmetric")
  meta_path<-paste0(prefix,".metadata.json")
  metadata<-NULL
  if(file.exists(meta_path)) {
    if(requireNamespace("jsonlite",quietly=TRUE))metadata<-jsonlite::fromJSON(meta_path) else {
      txt<-paste(readLines(meta_path,warn=FALSE),collapse="\n")
      metadata<-list(raw_json=txt)
      if(grepl('"r_scaling"\\s*:',txt,perl=TRUE))
        stop("Install jsonlite to read a bundle with moment scaling metadata")
      if(!grepl('"moment_type"\\s*:\\s*"cumulant"',txt,perl=TRUE))stop("Metadata must declare moment_type cumulant")
      version<-regmatches(txt,regexpr('"schema_version"\\s*:\\s*"[^"]+"',txt,perl=TRUE))
      if(!length(version)||!nzchar(version)||!grepl('"1.0"',version,fixed=TRUE))stop("Unsupported bundle schema version")
    }
    if(is.null(metadata$raw_json) && (is.null(metadata$moment_type)||metadata$moment_type!="cumulant"))stop("Only cumulant bundles are supported")
    if(is.null(metadata$raw_json) && (is.null(metadata$schema_version)||metadata$schema_version!="1.0"))stop("Unsupported bundle schema version")
  } else stop("Missing metadata.json: cumulant convention must be declared")
  traits<-unique(unlist(strsplit(m$key[m$order==1],"|",fixed=TRUE)))
  if(!length(traits))traits<-unique(unlist(strsplit(m$key,"|",fixed=TRUE)))
  if(!is.null(metadata$traits))traits<-metadata$traits
  log_path<-paste0(prefix,".log")
  jackknife<-NULL
  if(file.exists(paste0(prefix,".jackknife.tsv"))) {
    jk<-utils::read.delim(paste0(prefix,".jackknife.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
    if(!ncol(jk)||names(jk)[1L]!="block"||anyDuplicated(jk$block)||!setequal(names(jk)[-1L],m$key))
      stop("Invalid jackknife moment keys or block identifiers")
    jackknife<-as.matrix(jk[,m$key,drop=FALSE]);rownames(jackknife)<-as.character(jk$block)
    if(!is.numeric(jackknife)||any(!is.finite(jackknife)))stop("Invalid jackknife moment values")
  }
  scaling<-.restore_moment_scaling(metadata$r_scaling,traits)
  .enrich_moments(structure(list(moments=m,V=V,metadata=metadata,traits=traits,scaling=scaling,
    prefix=prefix,jackknife=jackknife,log=if(file.exists(log_path))readLines(log_path,warn=FALSE) else character()),
    class="gcmsem_moments"))
}

# JSON objects preserve trait names; earlier array encodings follow traits.
.restore_moment_scaling <- function(scaling,traits) {
  if(is.null(scaling))return(NULL)
  if(!is.list(scaling)||!all(c("scale","cumulative_scale","estimated_from_input","held_fixed")%in%names(scaling)))
    stop("Invalid moment scaling metadata")
  restore<-function(value) {
    if(is.list(value))value<-unlist(value,recursive=FALSE,use.names=TRUE)
    if(!is.numeric(value)||length(value)!=length(traits)||any(!is.finite(value))||any(value<=0))
      stop("Moment scaling metadata requires positive finite trait scales")
    if(is.null(names(value)))names(value)<-traits
    if(anyNA(names(value))||anyDuplicated(names(value))||!setequal(names(value),traits))
      stop("Moment scaling metadata does not match the trait names")
    value[traits]
  }
  for(nm in c("estimated_from_input","held_fixed"))
    if(!is.logical(scaling[[nm]])||length(scaling[[nm]])!=1L||is.na(scaling[[nm]]))
      stop("Invalid moment scaling flags")
  if(!isTRUE(scaling$held_fixed))stop("Only fixed-scale moment transformations are supported")
  list(scale=restore(scaling$scale),cumulative_scale=restore(scaling$cumulative_scale),
    estimated_from_input=scaling$estimated_from_input,held_fixed=TRUE)
}

# Add convenient second-order matrices without changing any moment estimate.
.enrich_moments <- function(x) {
  traits<-x$traits
  if(!is.character(traits)||!length(traits)||anyNA(traits)||anyDuplicated(traits)||any(!nzchar(traits)))
    stop("Moment traits must be unique nonempty names")
  m<-x$moments; V<-x$V; n<-length(traits)
  S<-S_SE<-R<-R_SE<-matrix(NA_real_,n,n,dimnames=list(traits,traits))
  key_matrix<-matrix(NA_character_,n,n,dimnames=list(traits,traits))
  tuples<-strsplit(m$key,"|",fixed=TRUE)
  if(any(!unlist(tuples)%in%traits))stop("Moment keys contain an undeclared trait")
  for(k in which(m$order==2L)) {
    ij<-match(tuples[[k]],traits);i<-ij[1L];j<-ij[2L]
    if(!is.na(key_matrix[i,j]))stop("Duplicate symmetric second-order moment")
    S[i,j]<-S[j,i]<-m$estimate[k]
    S_SE[i,j]<-S_SE[j,i]<-m$se[k]
    key_matrix[i,j]<-key_matrix[j,i]<-m$key[k]
  }
  for(i in seq_len(n))for(j in i:n) {
    if(!is.finite(S[i,i])||!is.finite(S[j,j])||S[i,i]<=0||S[j,j]<=0||!is.finite(S[i,j]))next
    if(i==j) {R[i,j]<-1;R_SE[i,j]<-0;next}
    r<-S[i,j]/sqrt(S[i,i]*S[j,j]);R[i,j]<-R[j,i]<-r
    keys<-c(key_matrix[i,j],key_matrix[i,i],key_matrix[j,j])
    gradient<-c(1/sqrt(S[i,i]*S[j,j]),-r/(2*S[i,i]),-r/(2*S[j,j]))
    vv<-as.numeric(crossprod(gradient,V[keys,keys,drop=FALSE]%*%gradient))
    # Round-off can produce a tiny negative quadratic form in a PSD V.
    scale<-sum(abs(outer(gradient,gradient)*V[keys,keys,drop=FALSE]))
    if(is.finite(vv) && vv>= -1e-10*max(scale,1e-300))R_SE[i,j]<-R_SE[j,i]<-sqrt(max(vv,0))
  }
  blocks<-x$metadata$jackknife_blocks
  n_blocks<-if(is.data.frame(blocks))nrow(blocks) else if(is.list(blocks))length(blocks) else NA_integer_
  x$S<-S;x$S_SE<-S_SE;x$R<-R;x$R_SE<-R_SE
  x$diagnostics<-list(n_traits=n,n_moments=nrow(m),orders=sort(unique(m$order)),
    n_snps=if(is.null(x$metadata$n_tags))NA_integer_ else x$metadata$n_tags,
    n_blocks=n_blocks,sampling_covariance_rank=x$metadata$sampling_covariance_rank,
    nonpositive_variances=traits[is.finite(diag(S))&diag(S)<=0],
    undefined_correlations=sum(is.na(R[upper.tri(R)])),
    correlation_se_method="Delta method using the full moment sampling covariance V")
  x
}

#' Extract parameter estimates
#' @export
parameter_estimates <- function(fit) {
  if(!inherits(fit,"gcmsem_fit"))stop("fit must be a gcmsem_fit")
  fit$parameters
}
#' @export
coef.gcmsem_fit <- function(object,...)object$coefficients
#' @export
vcov.gcmsem_fit <- function(object,...)object$vcov
#' @export
fitted.gcmsem_fit <- function(object,...)object$fitted
#' @export
coef.gcmsem_moments <- function(object,...)setNames(object$moments$estimate,object$moments$key)
#' @export
vcov.gcmsem_moments <- function(object,...)object$V
#' @export
as.data.frame.gcmsem_fit <- function(x,row.names=NULL,optional=FALSE,...) {
  ans<-x$parameters
  if(!is.null(row.names))rownames(ans)<-row.names
  ans
}
#' @export
as.data.frame.gcmsem_moments <- function(x,row.names=NULL,optional=FALSE,...) {
  ans<-x$moments
  if(!is.null(row.names))rownames(ans)<-row.names
  ans
}

.print_digits <- function(digits) {
  if(!is.numeric(digits)||length(digits)!=1L||!is.finite(digits)||digits<1L||digits>15L||digits!=as.integer(digits))
    stop("digits must be an integer between 1 and 15")
  as.integer(digits)
}
.yes_no <- function(x)if(isTRUE(x))"yes" else if(identical(x,FALSE))"no" else "unavailable"
.display_number <- function(x,digits=4L)trimws(formatC(x,format="g",digits=digits))
.matrix_preview <- function(x,max.traits) {
  if(!is.numeric(max.traits)||length(max.traits)!=1L||is.na(max.traits)||max.traits<1||
     (is.finite(max.traits)&&max.traits!=as.integer(max.traits)))
    stop("max.traits must be a positive integer or Inf")
  n<-min(nrow(x),max.traits)
  x[seq_len(n),seq_len(n),drop=FALSE]
}

.fit_header <- function(x,digits) {
  d<-x$diagnostics
  cat("GCMSEM ",x$estimator," fit\n",sep="")
  cat("Observed traits: ",paste(x$model$observed,collapse=", "),"\n",sep="")
  cat("Moments: ",d$n_moments,"   Free parameters: ",d$n_parameters,
      "   Jacobian rank: ",d$jacobian_rank,"/",d$n_parameters,"\n",sep="")
  cat("Converged: ",.yes_no(d$converged),"   Uncalibrated Q: ",.display_number(x$objective,digits),"\n",sep="")
  q<-x$fit_measures
  if(!is.null(q)) {
    cat("Q approximate SE: ",.display_number(q$Q_se,digits),
      "   Approximate null mean (SD): ",.display_number(q$Q_null_mean,digits),
      " (",.display_number(q$Q_null_sd,digits),")\n",sep="")
    cat("Q moment-matched p (approximate): ",.display_number(q$Q_p_approx,digits),
      "   Status: ",q$Q_status,"\n",sep="")
  }
  cat("Formal sandwich SEs: ",if(isTRUE(d$formal_se_available))"available" else "unavailable",
      "   Inference checks: ",if(isTRUE(d$inference_valid))"passed" else "flagged","\n",sep="")
  if(isTRUE(d$stability_required))cat("Stable branch: ",.yes_no(d$stable),"\n",sep="")
  if(isTRUE(d$source_admissible))cat("Source moment checks: passed\n") else {
    cat("*** HEYWOOD CASE: source moment checks failed ***\n")
    if(length(d$source$negative_variances))cat("Negative source variance: ",paste(d$source$negative_variances,collapse=", "),"\n",sep="")
    cat("Estimates and available formal SEs retained.\n")
  }
  if(length(d$inference_reasons))cat("Flags: ",paste(d$inference_reasons,collapse="; "),"\n",sep="")
  invisible(x)
}

.parameter_terms <- function(model) {
  ans<-setNames(rep("",length(model$start)),names(model$start))
  for(r in model$rows)if(r$parameter>0L) {
    nm<-names(model$start)[r$parameter];v<-model$nodes[r$tuple]
    term<-if(r$kind=="path")paste(v[1L],"~",v[2L]) else if(r$order==2L)
      paste(v[1L],"~~",v[2L]) else paste0("k",r$order,"(",paste(v,collapse=", "),")")
    ans[nm]<-if(nzchar(ans[nm]))paste(ans[nm],term,sep="; ") else term
  }
  ans
}

#' @export
summary.gcmsem_fit <- function(object,digits=4L,...) {
  digits<-.print_digits(digits)
  table<-object$parameters
  terms<-.parameter_terms(object$model)
  table<-data.frame(term=unname(terms[table$parameter]),table,check.names=FALSE,stringsAsFactors=FALSE)
  structure(list(call=object$call,model=object$model,estimator=object$estimator,
    objective=object$objective,fit_measures=object$fit_measures,parameters=table,diagnostics=object$diagnostics,
    starts=object$starts,digits=digits),class="summary_gcmsem_fit")
}
#' @export
print.gcmsem_fit <- function(x,digits=4L,...) {
  .fit_header(x,.print_digits(digits))
  invisible(x)
}
#' @export
print.summary_gcmsem_fit <- function(x,digits=x$digits,...) {
  digits<-.print_digits(digits);.fit_header(x,digits)
  cat("\nParameter estimates\n")
  tab<-x$parameters
  for(nm in intersect(c("estimate","se","z","p"),names(tab)))tab[[nm]]<-.display_number(tab[[nm]],digits)
  names(tab)[names(tab)=="parameter"]<-"label"
  names(tab)[names(tab)=="se"]<-"formal SE"
  print(tab,row.names=FALSE,right=TRUE)
  invisible(x)
}

.moments_header <- function(x) {
  d<-x$diagnostics
  cat("GCMSEM moments\n")
  cat("Traits: ",paste(x$traits,collapse=", "),"\n",sep="")
  if(length(d$n_snps)&&is.finite(d$n_snps))cat("SNPs: ",format(d$n_snps,big.mark=",",scientific=FALSE,trim=TRUE),"\n",sep="")
  counts<-table(x$moments$order)
  cat("Cumulants: ",nrow(x$moments)," (",paste(paste0("order ",names(counts),": ",as.integer(counts)),collapse=", "),")\n",sep="")
  cat("Sampling covariance V: ",nrow(x$V)," x ",ncol(x$V),sep="")
  if(length(d$n_blocks)&&is.finite(d$n_blocks))cat("   Jackknife groups: ",d$n_blocks,sep="")
  cat("\n")
  if(length(d$nonpositive_variances))cat("Correlation undefined for nonpositive variances: ",paste(d$nonpositive_variances,collapse=", "),"\n",sep="")
  invisible(x)
}
#' @export
print.gcmsem_moments <- function(x,digits=4L,max.traits=6L,...) {
  digits<-.print_digits(digits)
  if(is.null(x$R))x<-.enrich_moments(x)
  .moments_header(x)
  cat("\nGenetic correlations ($R)\n")
  shown<-.matrix_preview(x$R,max.traits)
  print(round(shown,digits),na.print="NA")
  if(nrow(shown)<nrow(x$R))cat("Showing ",nrow(shown)," of ",nrow(x$R)," traits; the full matrix is in $R.\n",sep="")
  invisible(x)
}
#' @export
summary.gcmsem_moments <- function(object,digits=4L,...) {
  digits<-.print_digits(digits)
  if(is.null(object$R))object<-.enrich_moments(object)
  structure(list(traits=object$traits,moments=object$moments,V=object$V,
    S=object$S,S_SE=object$S_SE,R=object$R,R_SE=object$R_SE,
    diagnostics=object$diagnostics,digits=digits),class="summary_gcmsem_moments")
}
#' @export
print.summary_gcmsem_moments <- function(x,digits=x$digits,max.traits=6L,...) {
  digits<-.print_digits(digits);.moments_header(x)
  for(nm in c("S","S_SE","R","R_SE")) {
    heading<-switch(nm,S="SNP-effect covariance",S_SE="Covariance SEs",R="Genetic correlations",R_SE="Correlation SEs (delta method using full V)")
    cat("\n",heading," ($",nm,")\n",sep="")
    shown<-.matrix_preview(x[[nm]],max.traits)
    print(signif(shown,digits),na.print="NA")
  }
  if(nrow(shown)<nrow(x$R))cat("Showing ",nrow(shown)," of ",nrow(x$R)," traits; all matrices are complete in the summary object.\n",sep="")
  invisible(x)
}
