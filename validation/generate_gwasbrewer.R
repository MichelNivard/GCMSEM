#!/usr/bin/env Rscript
# Actual GWASBrewer simulations. Large generated inputs are deliberately untracked.
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[1] else "validation/work"
nrep <- if (length(args)>1) as.integer(args[2]) else 20L
stopifnot(nrep > 0L)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(GWASBrewer))
data("ld_mat_list", package="GWASBrewer")
data("AF", package="GWASBrewer")
# Preserve the real HapMap3 chromosome-19 LD pattern and its varying block sizes.
# Recode BOTH LD and AF to the minor allele before drawing asymmetric effects.
sizes <- vapply(ld_mat_list, nrow, integer(1))
ends <- cumsum(sizes); starts <- c(1L, head(ends,-1L)+1L)
ld <- lapply(seq_along(sizes), function(i) {
  ix <- starts[i]:ends[i]; signs <- ifelse(AF[ix]>.5,-1,1)
  D <- Diagonal(x=signs)
  as(forceSymmetric(D %*% ld_mat_list[[i]] %*% D), "dsCMatrix")
})
af <- pmin(AF,1-AF)
lddir <- file.path(out,"ld"); dir.create(lddir,showWarnings=FALSE)
for(i in seq_along(ld)) {
  filename <- file.path(lddir,sprintf("block_%02d.mtx",i))
  valid <- FALSE
  if(file.exists(filename)) {
    contents <- readLines(filename,warn=FALSE)
    declared <- if(length(contents)>=2) as.numeric(strsplit(contents[2]," ")[[1]][3]) else NA_real_
    valid <- is.finite(declared) && length(contents)==declared+2
  }
  if(!valid) {
    temporary <- paste0(filename,".tmp")
    writeMM(ld[[i]], temporary)
    contents <- readLines(temporary,warn=FALSE)
    declared <- as.numeric(strsplit(contents[2]," ")[[1]][3])
    stopifnot(length(contents)==declared+2)
    stopifnot(file.rename(temporary,filename))
  }
}
write.table(data.frame(block=seq_along(sizes),size=sizes,start=starts,end=ends),
            file.path(lddir,"blocks.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(af=af),file.path(lddir,"af.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("Generating LD eigendecompositions once.\n")
eigfile <- file.path(out,"ld_eigen.rds")
if(file.exists(eigfile)) eig <- readRDS(eigfile) else {
  eig <- lapply(ld,function(x) eigen(as.matrix(x),symmetric=TRUE))
  saveRDS(eig,eigfile)
}
repeats <- 2L
J <- sum(sizes)*repeats
N <- if(length(args)>2) as.numeric(args[3]) else 50000
stopifnot(is.finite(N), N>0)
gamma_effect <- function(n,sd,snp_info) (rgamma(n,shape=2)-2)/sqrt(2)*sd/sqrt(n)
normal_effect <- function(n,sd,snp_info) rnorm(n)*sd/sqrt(n)
scenarios <- c(rep("chain",nrep), rep("null",max(3L,ceiling(nrep/4))),
               rep("factor",max(3L,ceiling(nrep/4))))
manifest <- list()
only_run <- Sys.getenv("GCMSEM_RUN_NAME", "")
for(k in seq_along(scenarios)) {
  scenario <- scenarios[k]; rep_id <- sum(scenarios[seq_len(k)]==scenario)
  name <- sprintf("%s_%03d",scenario,rep_id)
  if(nzchar(only_run) && name!=only_run)next
  seed <- 20261103L+k; set.seed(seed)
  p <- if(scenario=="factor") 4L else 3L
  G <- matrix(0,p,p); h2 <- rep(.2,p)
  if(scenario=="chain") { G[1,2]<-.4;G[2,3]<-.25 }
  if(scenario=="factor") G[1,2:4] <- .4
  # Fixed overlap, known observational correlation, and independent sampling
  # errors across the two replicated sets of genomic blocks.
  Nm <- matrix(N*.5,p,p);diag(Nm)<-N
  R_obs <- matrix(.2,p,p);diag(R_obs)<-1
  d <- sim_mv(N=Nm,J=J,h2=h2,pi=rep(.15,p),G=G,R_obs=R_obs,
              R_LD=ld,R_LD_eigen=eig,af=af,snp_effect_function=gamma_effect,
              h2_exact=FALSE,pi_exact=FALSE)
  keep <- if(scenario=="factor") 2:4 else 1:3
  trait_names <- c("A","B","C")
  sx <- sqrt(2*d$snp_info$AF*(1-d$snp_info$AF))
  # All exported coefficients are per allele. Python CLI also exercises conversion.
  dd <- data.frame(SNP=paste0("sim",seq_len(J)),A1="A",A2="C",AF=d$snp_info$AF)
  for(t in seq_along(keep)) {
    dd[[paste0("beta_",trait_names[t])]] <- d$beta_hat[,keep[t]]
    dd[[paste0("se_",trait_names[t])]] <- d$se_beta_hat[,keep[t]]
    dd[[paste0("joint_",trait_names[t])]] <- d$beta_joint[,keep[t]]*sx
  }
  filename <- file.path(out,paste0(name,".tsv.gz"))
  temporary <- paste0(filename,".tmp")
  con <- gzfile(temporary,"wt")
  write.table(dd,con,sep="\t",row.names=FALSE,quote=FALSE);close(con)
  con <- gzfile(temporary,"rt");n_written<-length(readLines(con));close(con)
  stopifnot(n_written==J+1L)
  stopifnot(file.rename(temporary,filename))
  write.table(d$R[keep,keep],file.path(out,paste0(name,".noise.tsv")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  # Expected direct effect variances come from the DAG specification, not the
  # realized effects. Their distribution includes the zero-effect SNPs.
  dir_h2 <- GWASBrewer:::check_G(G,h2)$dir_h2
  T <- solve(diag(p)-t(G))
  truth <- list(scenario=scenario,seed=seed,J=J,N=N,pi=.15,h2=h2,
                dir_h2=dir_h2,T=T[keep,,drop=FALSE],G=G,keep=keep,
                noise=d$R[keep,keep],source_k3=2/sqrt(2),source_excess=3)
  saveRDS(truth,file.path(out,paste0(name,".truth.rds")))
  write.table(T[keep,,drop=FALSE],file.path(out,paste0(name,".T.tsv")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(data.frame(dir_h2=dir_h2),file.path(out,paste0(name,".variance.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  manifest[[k]] <- data.frame(name=name,scenario=scenario,replicate=rep_id,seed=seed,
                              J=J,ld_patterns=length(ld),blocks=length(ld)*repeats,
                              min_block=min(sizes),max_block=max(sizes),N=N)
  cat(sprintf("Completed %s: %d SNPs and %d observed traits.\n",name,J,length(keep)))
}
if(!nzchar(only_run))write.table(do.call(rbind,manifest),file.path(out,"manifest.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"generation_session.txt"))
cat("All GWASBrewer generation complete.\n")
