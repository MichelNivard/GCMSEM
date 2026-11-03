#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
work <- if(length(args))args[1] else "validation/work"
out <- if(length(args)>1)args[2] else "validation/results"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Source the checked-out implementation so this script always tests this commit.
# Installed users may replace the next line by library(GCMSEM).
source_dir <- Sys.getenv("GCMSEM_SOURCE_DIR", "R")
for(f in list.files(source_dir,pattern="\\.R$",full.names=TRUE))source(f)
manifest <- read.delim(file.path(work,"manifest.tsv"),stringsAsFactors=FALSE)
fits <- list(); pars <- list(); exact <- list(); alternatives <- list()
make_model <- function(scenario,gaussian=FALSE) {
  syntax <- switch(scenario,
    chain="B ~ ab*A\nC ~ bc*B",
    null="A ~~ A\nB ~~ B\nC ~~ C",
    factor="F =~ a*A + a*B + a*C\nF ~~ 1*F")
  gcmsem_model(syntax,observed=c("A","B","C"),
               gaussian_observed=gaussian,gaussian_latent=gaussian)
}
for(i in seq_len(nrow(manifest))) {
  rec <- manifest[i,]; prefix <- file.path(work,"fits",rec$name)
  if(!file.exists(paste0(prefix,".moments.tsv")))next
  cat("Fitting",rec$name,"\n")
  result <- tryCatch({
    mm <- read_moments(prefix); model <- make_model(rec$scenario)
    fit <- gcmsem(model,mm,n_starts=if(rec$replicate==1)8 else 3,seed=rec$seed,
                  control=list(iter.max=1500L,eval.max=3000L,rel.tol=1e-11),
                  log_file=file.path(work,"fits",paste0(rec$name,".sem.log")),verbose=FALSE)
    dg <- fit$diagnostics
    fits[[length(fits)+1L]] <- data.frame(run=rec$name,scenario=rec$scenario,
      replicate=rec$replicate,status="ok",Q=fit$objective,converged=dg$converged,
      inference_valid=dg$inference_valid,formal_se_available=dg$formal_se_available,
      source_admissible=dg$source_admissible,jacobian_rank=dg$jacobian_rank,
      n_parameters=dg$n_parameters,gradient_max=dg$gradient_max,
      inference_reasons=paste(dg$inference_reasons,collapse="; "))
    pp <- parameter_estimates(fit);pp$run<-rec$name;pp$scenario<-rec$scenario
    pp$truth <- NA_real_
    if(rec$scenario=="chain")pp$truth[match(c("ab","bc"),pp$parameter)]<-c(.4,.25)
    if(rec$scenario=="factor")pp$truth[pp$parameter=="a"]<-.4
    pp$covered<-abs(pp$estimate-pp$truth)<=1.96*pp$se
    pars[[length(pars)+1L]]<-pp
    if(rec$replicate==1) {
      target <- scan(paste0(prefix,".truth.tsv"),quiet=TRUE)
      names(target)<-mm$moments$key
      oracle <- gcmsem(model,target,V=diag(length(target)),n_starts=8,seed=rec$seed,verbose=FALSE,
                        control=list(iter.max=3000L,eval.max=6000L,rel.tol=1e-12))
      # Identity weights here are a numerical oracle check, not sampling data.
      # Do not publish the resulting arbitrary SEs/p-values as inference.
      op<-parameter_estimates(oracle)[,c("parameter","estimate")]
      op$scenario<-rec$scenario;op$Q<-oracle$objective
      stopifnot(oracle$objective < 1e-9)
      if(rec$scenario=="chain")stopifnot(max(abs(coef(oracle)[c("ab","bc")]-c(.4,.25)))<1e-5)
      if(rec$scenario=="factor")stopifnot(abs(abs(coef(oracle)["a"])-.4)<1e-5)
      exact[[length(exact)+1L]]<-op
      gauss <- gcmsem(make_model(rec$scenario,TRUE),mm,n_starts=3,seed=rec$seed,verbose=FALSE)
      # Covariance-only fitting is an independent observed-moment toggle.
      # When unused higher-source parameters remain free rank diagnostics should
      # withhold ordinary uncertainty, rather than silently change the model.
      covonly <- gcmsem(model,mm,orders=2L,n_starts=3,seed=rec$seed,verbose=FALSE)
      alternatives[[length(alternatives)+1L]] <- data.frame(scenario=rec$scenario,
        nonGaussian_Q=fit$objective,Gaussian_Q=gauss$objective,
        covariance_only_Q=covonly$objective,covariance_only_rank=covonly$diagnostics$jacobian_rank,
        covariance_only_parameters=covonly$diagnostics$n_parameters,
        covariance_only_inference_valid=covonly$diagnostics$inference_valid)
    }
    TRUE
  },error=function(e){
    fits[[length(fits)+1L]] <<- data.frame(run=rec$name,scenario=rec$scenario,
      replicate=rec$replicate,status=paste("failed:",conditionMessage(e)),Q=NA,converged=FALSE,
      inference_valid=FALSE,formal_se_available=FALSE,source_admissible=FALSE,
      jacobian_rank=NA,n_parameters=NA,gradient_max=NA,inference_reasons=conditionMessage(e))
    message(rec$name,": ",conditionMessage(e));FALSE
  })
  write.table(do.call(rbind,fits),file.path(out,"sem_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  if(length(pars))write.table(do.call(rbind,pars),file.path(out,"sem_parameters.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  if(length(exact))write.table(do.call(rbind,exact),file.path(out,"exact_model_recovery.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  if(length(alternatives))write.table(do.call(rbind,alternatives),file.path(out,"moment_toggles.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
capture.output(sessionInfo(),file=file.path(out,"fitting_session.txt"))
cat("All GCMSEM fits complete.\n")
