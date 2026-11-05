# Interactive helpers for the archived matched-input GCMSEM / S-LD4M diagnostic.
# Use original_polygenicity.R for the author-reference pipeline comparison.
# Run from the GCMSEM repository root. Results stay in your R session.

read_polygenicity_comparison <- function(out = "tutorial-results/polygenicity") {
  if (!requireNamespace("GCMSEM", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE))
    stop("Install GCMSEM and jsonlite")
  comparison <- read.delim(file.path(out, "comparison.tsv"), check.names = FALSE)
  strata <- read.delim(file.path(out, "strata.tsv"), check.names = FALSE)
  moments <- setNames(lapply(comparison$trait, function(trait)
    GCMSEM::read_moments(file.path(out, trait))), comparison$trait)
  structure(list(comparison = comparison, strata = strata, moments = moments,
    analysis = jsonlite::fromJSON(file.path(out, "analysis.json")),
    agreement = jsonlite::fromJSON(file.path(out, "comparison.json"))),
    class = "gcmsem_polygenicity_comparison")
}

run_polygenicity_example <- function(data_dir = "tutorial-data",
    out = "tutorial-results/polygenicity", sld4m_dir = "external/SLD4M",
    python = Sys.which("python"), octave = Sys.which("octave")) {
  if (!nzchar(python) || !nzchar(octave)) stop("Provide the Python and Octave executable paths")
  expected <- "dd947104cb3cb86bb8a651bf3b5abf5e5cbf46f3"
  actual <- system2("git", c("-C", shQuote(sld4m_dir), "rev-parse", "HEAD"), stdout = TRUE)
  if (!identical(actual, expected)) stop("Check out the pinned SLD4M commit shown in the tutorial")
  if (system2("git", c("-C", shQuote(sld4m_dir), "diff", "--quiet", "HEAD", "--", "SLD4M.m")) != 0)
    stop("Use an unmodified upstream SLD4M.m")
  run <- function(exe, args) {
    status <- system2(exe, args)
    if (status != 0L) stop("Analysis command failed with exit status ", status)
  }
  run(python, c("tutorials/prepare_polygenicity.py", "--data-dir", shQuote(data_dir), "--out", shQuote(out)))
  matlab_string <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
  code <- paste0("addpath('tutorials'); run_sld4m(", matlab_string(file.path(out,"sld4m_input.mat")),
    ",", matlab_string(sld4m_dir), ",", matlab_string(file.path(out,"sld4m.mat")), ");")
  run(octave, c("--no-gui", "--quiet", "--eval", shQuote(code)))
  zcode <- paste0("addpath('tutorials'); run_sld4m(", matlab_string(file.path(out,"sld4m_input.mat")),
    ",", matlab_string(sld4m_dir), ",", matlab_string(file.path(out,"sld4m_zscore.mat")), ",'zscore');")
  run(octave, c("--no-gui", "--quiet", "--eval", shQuote(zcode)))
  run(python, c("tutorials/summarize_polygenicity.py", "--out", shQuote(out), "--sld4m", shQuote(sld4m_dir)))
  read_polygenicity_comparison(out)
}

print.gcmsem_polygenicity_comparison <- function(x, ...) {
  cat("GCMSEM and matched-input S-LD4M: effective polygenicity\n")
  empirical <- is.null(x$analysis$ld_reference)
  cat(nrow(x$comparison), "traits;", format(x$analysis$n_tags, big.mark=","),
      if(empirical) "common HM3 tags\n" else "simulated source/tag variants\n")
  cat(if(empirical) "Source-weighted MAF-bin aggregation; HM3 proxy for source LD averages\n" else "Exact generating LD; every source is also a tag\n")
  cat("S-LD4M input:", x$agreement$sld4m_input, "\n")
  print(x$comparison[, c("trait","gcmsem","gcmsem_se","sld4m","sld4m_se")],
        digits=4, row.names=FALSE)
  cat("\nPearson r:", format(x$agreement$pearson, digits=4), "\n")
  cat("RMSE:", format(x$agreement$rmse, digits=5), "effective variants\n")
  flagged <- x$comparison$flags != "none"
  cat("Traits with one or more diagnostic flags:", sum(flagged), "of", length(flagged), "\n")
  cat("Inspect $comparison$flags and $strata; signed estimates are retained.\n")
  invisible(x)
}

plot.gcmsem_polygenicity_comparison <- function(x, uncertainty = TRUE, input = c("primary", "beta"), ...) {
  input <- match.arg(input)
  d <- x$comparison
  # Linear axes preserve negative estimates. Error bars are +/- one reported
  # jackknife SE, not confidence intervals for these potentially unstable ratios.
  xx <- d[[if(input=="beta") "sld4m_beta" else "sld4m"]]/1000; yy <- d$gcmsem/1000
  xe <- if (uncertainty) d[[if(input=="beta") "sld4m_beta_se" else "sld4m_se"]]/1000 else rep(0,length(xx))
  ye <- if (uncertainty) d$gcmsem_se/1000 else rep(0,length(yy))
  lim <- range(c(0,xx-xe,xx+xe,yy-ye,yy+ye), finite=TRUE)
  flagged <- d[[if(input=="beta") "flags_beta" else "flags"]] != "none"
  col <- ifelse(flagged, "#B45309", "#16697A")
  plot(xx, yy, xlim=lim, ylim=lim, asp=1, pch=ifelse(flagged,17,19), col=col,
       xlab=if (identical(x$agreement$sld4m_input,"author reference"))
         "S-LD4M with author reference (thousands)" else
         paste0("S-LD4M with shared inputs, ",if(input=="beta") "BETA squared" else x$agreement$sld4m_input," (thousands)"),
       ylab="GCMSEM (thousands of effective variants)",
       main=paste(nrow(d), if(uncertainty) "traits: estimates and one-SE bars" else "traits: effective polygenicity"), ...)
  abline(0,1,lty=2,col="grey45")
  if (uncertainty) {
    segments(xx-xe,yy,xx+xe,yy,col=adjustcolor(col,.5))
    segments(xx,yy-ye,xx,yy+ye,col=adjustcolor(col,.5))
  }
  points(xx,yy,pch=ifelse(flagged,17,19),col=col)
  # Small deterministic label repulsion keeps clustered traits readable.
  usr <- par("usr"); pad <- diff(usr[1:2])*.012
  widths <- strwidth(d$trait,cex=.70); height <- strheight("M",cex=.70)*1.5
  lx <- pmin(pmax(xx+pad,usr[1]+pad),usr[2]-widths-pad); ly <- yy
  for(iter in seq_len(200)) for(i in seq_along(lx)) for(j in seq_len(i-1L)) {
    overlap <- lx[i] < lx[j]+widths[j]+pad && lx[j] < lx[i]+widths[i]+pad
    if(overlap && abs(ly[i]-ly[j]) < height) {
      direction <- if(ly[i] >= ly[j]) 1 else -1
      move <- (height-abs(ly[i]-ly[j]))/2+.001*height
      ly[i] <- ly[i]+direction*move; ly[j] <- ly[j]-direction*move
    }
  }
  ly <- pmin(pmax(ly,usr[3]+height),usr[4]-height)
  segments(xx,yy,lx,ly,col=adjustcolor(col,.5))
  text(lx,ly,labels=d$trait,adj=c(0,.5),cex=.70,col=col)
  legend("topleft",c("No diagnostic flag","One or more flags","Equality"),
         pch=c(19,17,NA),lty=c(NA,NA,2),col=c("#16697A","#B45309","grey45"),bty="n",cex=.8)
  invisible(x)
}

# An algebraic transformation of a pooled cumulant object. This has NO LD term.
# It is a zero-mean Gaussian-mixture equivalent count, not a causal-SNP count.
independent_effective_count <- function(moments, trait) {
  k2 <- paste(rep(trait,2),collapse="|"); k4 <- paste(rep(trait,4),collapse="|")
  keys <- c(k2,k4); ii <- match(keys,moments$moments$key)
  if(anyNA(ii)) stop("The object must contain second and fourth cumulants for this trait")
  M <- moments$metadata$n_reference_variants
  if(length(M)!=1L || !is.finite(M) || M<=0) stop("Missing source-variant count")
  v <- moments$moments$estimate[ii[1]]; k <- moments$moments$estimate[ii[2]]
  den <- k+3*v^2; estimate <- 3*M*v^2/den
  gradient <- c(6*M*v*k/den^2,-3*M*v^2/den^2)
  variance <- as.numeric(crossprod(gradient,moments$V[keys,keys,drop=FALSE]%*%gradient))
  data.frame(trait=trait,source_variants=M,excess_kurtosis=k/v^2,
             independent_effective_count=estimate,delta_se=sqrt(max(variance,0)),
             interpretation="No LD adjustment; assumes a zero-mean Gaussian variance mixture")
}
