#!/usr/bin/env Rscript
# Fits each real two-trait bundle with standalone GCMSEM. No MCMSEM dependency.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 4L)
  stop("Usage: Rscript tutorials/fit_causal.R PREFIX OUTDIR [N_STARTS=32] [SEED=20261105]")
library(GCMSEM)
prefix <- args[1L]; outdir <- args[2L]
n_starts <- if(length(args)>=3L) as.integer(args[3L]) else 32L
seed <- if(length(args)>=4L) as.integer(args[4L]) else 20261105L
if(is.na(n_starts) || n_starts<1L || is.na(seed)) stop("Invalid starts or seed")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
input <- read_moments(prefix)
if(length(input$traits)!=2L) stop("This tutorial expects exactly two traits")
x <- input$traits[1L]; y <- input$traits[2L]

# Numerical conditioning: divide every cumulant index by the same fixed SNP
# effect SD. Apply the identical linear transform to both axes of V. These
# constants are held fixed throughout fitting and the sandwich calculation.
m <- setNames(input$moments$estimate, input$moments$key)
sd_effect <- setNames(sqrt(m[paste(input$traits, input$traits, sep="|")]), input$traits)
if(any(!is.finite(sd_effect)) || any(sd_effect<=0))
  stop("Positive marginal effect variances are required for this numerical rescaling")
multiplier <- vapply(strsplit(names(m), "|", fixed=TRUE),
                     function(v) 1/prod(sd_effect[v]), numeric(1))
names(multiplier) <- names(m)
m <- m * multiplier
V <- input$V * outer(multiplier, multiplier)
write.table(data.frame(trait=input$traits, fixed_effect_sd=sd_effect),
            file.path(outdir,"scaling.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

source_syntax <- paste(
  sprintf("%s ~~ variance_x*%s",x,x), sprintf("%s ~~ variance_y*%s",y,y),
  sprintf('k3(%s, %s, %s) = free("skew_x")',x,x,x),
  sprintf('k3(%s, %s, %s) = free("skew_y")',y,y,y),
  sprintf('k4(%s, %s, %s, %s) = free("excess_x")',x,x,x,x),
  sprintf('k4(%s, %s, %s, %s) = free("excess_y")',y,y,y,y),sep="\n")
common_source <- paste(sprintf("F =~ 1*%s + 1*%s",x,y),
                       "F ~~ shared_variance*F", "k3(F, F, F) = 0", "k4(F, F, F, F) = 0",sep="\n")
specs <- list(
  forward_gaussian_confounder=paste(sprintf("%s ~ forward*%s",y,x),common_source,source_syntax,sep="\n"),
  reverse_gaussian_confounder=paste(sprintf("%s ~ reverse*%s",x,y),common_source,source_syntax,sep="\n"),
  reciprocal_no_confounder=paste(sprintf("%s ~ forward*%s",y,x),
                                sprintf("%s ~ reverse*%s",x,y),source_syntax,sep="\n"))
writeLines(c("Equal common-source direct loadings are imposed on the fixed SNP-effect-SD scale.",
  "This is an explicit modeling restriction, not an invariant across trait rescalings.",
  "The equal factor loadings are fixed at one and its variance is freely estimated.",
  "A negative common-source variance is retained and prominently flagged as a Heywood case.",
  "This variance parameterization avoids imposing a nonnegative squared-loading restriction.",
  "The reciprocal no-confounder model is a prespecified sensitivity fit; it is always run.",
  "Source moments and residual variances are unconstrained; infeasible estimates are retained.",
  "Only the reciprocal model requests the dynamically stable branch (spectral radius < 1).",
  "These models estimate associations between SNP-effect distributions. Causal interpretation",
  "requires assumptions about independent source effects, shared causes and homogeneous effects.",
  "Q is descriptive and uncalibrated; differences in Q do not yield a causal-direction p-value.",
  "Formal sandwich SEs do not resolve source-assignment equivalence, model uncertainty or",
  "the known higher-order interval-coverage limitations of this research release."),
  file.path(outdir,"interpretation.txt"))
status <- list(); paths <- list()
write_tsv <- function(z, name) write.table(z,file.path(outdir,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
for(i in seq_along(specs)) {
  name <- names(specs)[i]
  writeLines(specs[[i]],file.path(outdir,paste0(name,".model.txt")))
  model <- gcmsem_model(specs[[i]],observed=c(x,y),gaussian_observed=FALSE,gaussian_latent=TRUE)
  fit <- gcmsem(model,m,V,orders=2:4,n_starts=n_starts,seed=seed+i-1L,
    stable=(name=="reciprocal_no_confounder"),feasibility="check",
    control=list(iter.max=5000L,eval.max=10000L,rel.tol=1e-11),
    log_file=file.path(outdir,paste0(name,".log")))
  saveRDS(fit,file.path(outdir,paste0(name,".fit.rds")))
  write_tsv(parameter_estimates(fit),paste0(name,".parameters.tsv"))
  write_tsv(fit$starts,paste0(name,".starts.tsv"))
  write_tsv(data.frame(key=names(fit$observed), observed=fit$observed, implied=fit$fitted,
                      residual=fit$residuals),paste0(name,".fitted.tsv"))
  dput(fit$diagnostics,file.path(outdir,paste0(name,".diagnostics.R")))
  d <- fit$diagnostics
  status[[name]] <- data.frame(model=name,Q=fit$objective,converged=d$converged,
    jacobian_rank=d$jacobian_rank,n_parameters=d$n_parameters,n_moments=d$n_moments,
    source_admissible=d$source_admissible,formal_se_available=d$formal_se_available,
    inference_valid=d$inference_valid,stable=d$stable,stability_required=d$stability_required,
    spectral_radius=d$spectral_radius,gradient_max=d$gradient_max,
    reasons=paste(d$inference_reasons,collapse="; "),stringsAsFactors=FALSE)
  for(path in intersect(c("forward","reverse"),names(coef(fit)))) {
    from <- if(path=="forward") x else y
    to <- if(path=="forward") y else x
    factor <- unname(sd_effect[to]/sd_effect[from])
    p <- fit$parameters[fit$parameters$parameter==path,,drop=FALSE]
    paths[[paste(name,path)]] <- data.frame(model=name,from=from,to=to,
      effect_sd_estimate=p$estimate,effect_sd_formal_se=p$se,
      original_scale_estimate=p$estimate*factor,original_scale_formal_se=p$se*factor,
      se_status=p$se_status,source_admissible=d$source_admissible,
      inference_valid=d$inference_valid,stringsAsFactors=FALSE)
  }
  # The reciprocal model has an exactly equivalent source-swapped solution.
  # Report it explicitly; stability chooses a branch, not a better-fitting model.
  if(name=="reciprocal_no_confounder" && all(abs(coef(fit)[c("forward","reverse")])>1e-10)) {
    p <- coef(fit); alt <- p; b <- p[["forward"]]; c <- p[["reverse"]]
    alt["forward"] <- 1/c; alt["reverse"] <- 1/b
    for(order in 2:4) {
      left <- c("variance_x","skew_x","excess_x")[order-1L]
      right <- c("variance_y","skew_y","excess_y")[order-1L]
      alt[left] <- (-1/b)^order*p[right]
      alt[right] <- (-1/c)^order*p[left]
    }
    alt_m <- implied_moments(model,alt,keys=names(fit$observed))
    write_tsv(data.frame(parameter=names(alt),estimate=as.numeric(alt)),"source_swapped.parameters.tsv")
    write_tsv(data.frame(max_absolute_implied_difference=max(abs(alt_m-fit$fitted)),
      relative_implied_difference=max(abs(alt_m-fit$fitted))/max(1,max(abs(fit$fitted))),
      selected_spectral_radius=d$spectral_radius,
      swapped_spectral_radius=sqrt(abs(alt[["forward"]]*alt[["reverse"]]))),"source_swapped.equivalence.tsv")
  }
  summary(fit)
}
write_tsv(do.call(rbind,status),"model_status.tsv")
write_tsv(do.call(rbind,paths),"paths.tsv")
capture.output(sessionInfo(),file=file.path(outdir,"session-info.txt"))
