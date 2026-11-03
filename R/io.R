#' Read a Python GCMSEM moment bundle
#' @export
read_moments <- function(prefix) {
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
  meta_path<-paste0(prefix,".metadata.json")
  metadata<-NULL
  if(file.exists(meta_path)) {
    if(requireNamespace("jsonlite",quietly=TRUE))metadata<-jsonlite::fromJSON(meta_path) else {
      txt<-paste(readLines(meta_path,warn=FALSE),collapse="\n")
      # Keep complete JSON verbatim without pretending to parse arbitrary JSON.
      metadata<-list(raw_json=txt)
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
  structure(list(moments=m,V=V,metadata=metadata,traits=traits,prefix=prefix),class="gcmsem_moments")
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
summary.gcmsem_fit <- function(object,...) {
  print(object)
  cat("\nParameter estimates (source cumulants; order 4 is excess):\n")
  print(object$parameters,row.names=FALSE)
  cat("\nSource moment conditions:\n")
  print(object$diagnostics$source)
  invisible(object)
}
#' @export
print.gcmsem_fit <- function(x,...) {
  cat("GCMSEM",x$estimator,"fit\n")
  cat("Discrepancy Q:",format(x$objective,digits=7),"(uncalibrated)\n")
  d<-x$diagnostics
  cat("Moments:",d$n_moments," Parameters:",d$n_parameters," Jacobian rank:",d$jacobian_rank,"\n")
  cat("Converged:",d$converged," Stable:",d$stable," Conditional inference:",d$inference_valid,"\n")
  cat("Source admissible:",d$source_admissible," Formal sandwich SEs available:",d$formal_se_available,"\n")
  if(!d$inference_valid)cat("Inference flagged (point estimates and computable formal SEs retained):",paste(d$inference_reasons,collapse="; "),"\n")
  cat("Local identification and moment-feasibility checks do not establish a unique causal interpretation.\n")
  invisible(x)
}
