# Explicit file export is separate from the ordinary in-memory R workflow.
.validate_export_options <- function(write.to.tsv,out,overwrite) {
  if(!is.logical(write.to.tsv)||length(write.to.tsv)!=1L||is.na(write.to.tsv))stop("write.to.tsv must be TRUE or FALSE")
  if(!is.logical(overwrite)||length(overwrite)!=1L||is.na(overwrite))stop("overwrite must be TRUE or FALSE")
  if(write.to.tsv && (!is.character(out)||length(out)!=1L||is.na(out)||!nzchar(out)))
    stop("Supply an output prefix in out when write.to.tsv = TRUE")
  invisible(NULL)
}

.export_table <- function(x,path) {
  x<-as.data.frame(x,check.names=FALSE,stringsAsFactors=FALSE)
  # Seventeen significant digits preserve binary-double values on re-reading.
  for(nm in names(x))if(is.numeric(x[[nm]])&&!is.integer(x[[nm]]))
    x[[nm]]<-trimws(formatC(x[[nm]],format="g",digits=17L))
  utils::write.table(x,path,sep="\t",quote=TRUE,row.names=FALSE,na="NA")
}
.export_matrix <- function(x,label) {
  rn<-rownames(x)
  if(is.null(rn))rn<-as.character(seq_len(nrow(x)))
  ans<-data.frame(setNames(list(rn),label),x,check.names=FALSE,stringsAsFactors=FALSE)
  ans
}

.export_plan <- function(x,prefix) {
  if(inherits(x,"gcmsem_moments")) {
    files<-paste0(prefix,c(".moments.tsv",".sampling_cov.tsv",".covariance.tsv",
      ".covariance_se.tsv",".correlation.tsv",".correlation_se.tsv",".metadata.json"))
    if(!is.null(x$jackknife))files<-c(files,paste0(prefix,".jackknife.tsv"))
    if(length(x$log))files<-c(files,paste0(prefix,".log"))
    children<-x$conditional
    if(length(children)) {
      if(is.null(names(children))||anyNA(names(children))||anyDuplicated(names(children))||
         any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",names(children))))
        stop("Conditional bundles require unique simple names")
      for(nm in names(children))files<-c(files,.export_plan(children[[nm]],paste0(prefix,".",nm)))
    }
  } else if(inherits(x,"gcmsem_fit")) {
    files<-paste0(prefix,c(".parameters.tsv",".fitted.tsv",".sampling_cov.tsv",
      ".vcov.tsv",".starts.tsv",".diagnostics.tsv",".model.txt"))
    if(!is.null(x$fit_measures))files<-c(files,paste0(prefix,".fit_measures.tsv"))
    if(length(x$log))files<-c(files,paste0(prefix,".log"))
  } else stop("x must be a gcmsem_fit or gcmsem_moments object")
  files
}

# These optional members can otherwise survive an overwrite and be attached
# to a later analysis by read_moments(). Never remove arbitrary prefix matches.
.export_optional_paths <- function(x,prefix) {
  paths<-paste0(prefix,".log")
  if(inherits(x,"gcmsem_moments")) {
    paths<-c(paths,paste0(prefix,".jackknife.tsv"))
    if(length(x$conditional))for(nm in names(x$conditional))
      paths<-c(paths,.export_optional_paths(x$conditional[[nm]],paste0(prefix,".",nm)))
  }
  paths
}

.structured_export_metadata <- function(x) {
  is.null(x$metadata$raw_json)||length(setdiff(names(x$metadata),"raw_json"))>0L||
    !is.null(x$scaling)||length(x$conditional)>0L
}

.write_export <- function(x,prefix) {
  if(inherits(x,"gcmsem_moments")) {
    x<-.enrich_moments(x)
    .export_table(x$moments,paste0(prefix,".moments.tsv"))
    .export_table(.export_matrix(x$V,"key"),paste0(prefix,".sampling_cov.tsv"))
    for(nm in c("S","S_SE","R","R_SE")) {
      suffix<-switch(nm,S=".covariance.tsv",S_SE=".covariance_se.tsv",R=".correlation.tsv",R_SE=".correlation_se.tsv")
      .export_table(.export_matrix(x[[nm]],"trait"),paste0(prefix,suffix))
    }
    if(!is.null(x$jackknife)).export_table(.export_matrix(x$jackknife,"block"),paste0(prefix,".jackknife.tsv"))
    meta<-x$metadata
    if(!.structured_export_metadata(x))writeLines(meta$raw_json,paste0(prefix,".metadata.json"),useBytes=TRUE) else {
      if(!is.null(meta$raw_json)) {
        parsed<-jsonlite::fromJSON(meta$raw_json)
        for(nm in setdiff(names(meta),"raw_json"))parsed[[nm]]<-meta[[nm]]
        meta<-parsed
      }
      scaling<-if(!is.null(x$scaling))x$scaling else meta$r_scaling
      if(!is.null(scaling)) {
        scaling<-.restore_moment_scaling(scaling,x$traits)
        # Named numeric vectors are arrays by default in jsonlite. Lists
        # encode objects and retain trait labels, including a one-trait model.
        scaling$scale<-as.list(scaling$scale)
        scaling$cumulative_scale<-as.list(scaling$cumulative_scale)
        meta$r_scaling<-scaling
      }
      meta$schema_version<-"1.0";meta$moment_type<-"cumulant";meta$traits<-x$traits
      meta$output_prefix<-prefix
      meta$output_files<-basename(.export_plan(x,prefix))
      meta$conditional_prefixes<-if(length(x$conditional))paste0(basename(prefix),".",names(x$conditional)) else list()
      jsonlite::write_json(meta,paste0(prefix,".metadata.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA,null="null",na="null")
    }
    if(length(x$conditional))for(nm in names(x$conditional)).write_export(x$conditional[[nm]],paste0(prefix,".",nm))
  } else {
    .export_table(x$parameters,paste0(prefix,".parameters.tsv"))
    .export_table(data.frame(key=names(x$observed),observed=as.numeric(x$observed),
      fitted=as.numeric(x$fitted),residual=as.numeric(x$residuals)),paste0(prefix,".fitted.tsv"))
    .export_table(.export_matrix(x$V,"key"),paste0(prefix,".sampling_cov.tsv"))
    .export_table(.export_matrix(x$vcov,"parameter"),paste0(prefix,".vcov.tsv"))
    .export_table(x$starts,paste0(prefix,".starts.tsv"))
    if(!is.null(x$fit_measures)).export_table(x$fit_measures,paste0(prefix,".fit_measures.tsv"))
    d<-x$diagnostics
    fields<-c("converged","jacobian_rank","n_parameters","n_moments","residual_df","V_rank",
      "sandwich_rank","gradient_max","dwls_gradient_max","spectral_radius","stable",
      "stability_required","source_admissible","formal_se_available","inference_valid",
      "fit_statistic_calibrated","feasibility_method","condition_number")
    values<-vapply(fields,function(nm)if(is.null(d[[nm]]))NA_character_ else paste(d[[nm]],collapse="; "),character(1))
    diag_table<-data.frame(diagnostic=c(fields,"inference_reasons","formal_se_reasons","negative_source_variances"),
      value=c(values,paste(d$inference_reasons,collapse="; "),paste(d$formal_se_reasons,collapse="; "),
        paste(d$source$negative_variances,collapse="; ")),stringsAsFactors=FALSE)
    .export_table(diag_table,paste0(prefix,".diagnostics.tsv"))
    writeLines(x$model$syntax,paste0(prefix,".model.txt"),useBytes=TRUE)
  }
  if(length(x$log))writeLines(x$log,paste0(prefix,".log"),useBytes=TRUE)
  invisible(NULL)
}

#' Export an existing GCMSEM object
#' @export
export_gcmsem <- function(x,prefix,overwrite=FALSE) {
  .validate_export_options(TRUE,prefix,overwrite)
  files<-.export_plan(x,prefix)
  if(anyDuplicated(files))stop("Export paths are not unique")
  managed<-unique(c(files,.export_optional_paths(x,prefix)))
  if(any(dir.exists(managed)))stop("An export path is an existing directory")
  existing<-managed[file.exists(managed)]
  if(length(existing)&&!overwrite)stop("Refusing to overwrite existing files: ",paste(existing,collapse=", "))
  needs_json<-function(z)inherits(z,"gcmsem_moments") &&
    (.structured_export_metadata(z)||any(vapply(z$conditional,needs_json,logical(1))))
  if(needs_json(x)&&!requireNamespace("jsonlite",quietly=TRUE))
    stop("Install jsonlite to export structured moment metadata")
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(dirname(prefix)))stop("Could not create export directory")
  .write_export(x,prefix)
  stale<-setdiff(existing,files)
  if(length(stale)&&!all(file.remove(stale)))stop("Could not remove stale optional export files")
  invisible(files)
}
