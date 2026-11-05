# Run from the repository root. Computed results are returned as R objects.
source("tutorials/polygenicity.R")

read_original_polygenicity <- function(out = "tutorial-results/original-polygenicity",
                                      gcmsem_out = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite")
  d <- read.delim(file.path(out,"comparison.tsv"),check.names=FALSE)
  moments <- if (is.null(gcmsem_out)) NULL else setNames(lapply(d$trait, function(trait)
    GCMSEM::read_moments(file.path(gcmsem_out,trait))),d$trait)
  structure(list(comparison=d,
    annotations=read.delim(file.path(out,"annotations.tsv"),check.names=FALSE),
    moments=moments, agreement=jsonlite::fromJSON(file.path(out,"comparison.json")),
    provenance=jsonlite::fromJSON(file.path(out,"provenance.json"))),
    class="gcmsem_original_polygenicity")
}

run_original_polygenicity <- function(data_dir="tutorial-data",
    gcmsem_out="tutorial-results/gcmsem-polygenicity",
    native_dir="tutorial-results/original-sld4m",
    out="tutorial-results/original-polygenicity", sld4m_dir="external/SLD4M",
    python=Sys.which("python"), octave=Sys.which("octave"), jobs=1L,
    compute_gcmsem=TRUE) {
  if (!nzchar(python) || !nzchar(octave)) stop("Provide Python and Octave paths")
  run <- function(exe,args) {
    status <- system2(exe,args)
    if (status != 0L) stop("Analysis failed with status ",status)
  }
  if (compute_gcmsem) run(python,c("tutorials/prepare_polygenicity.py",
      "--data-dir",shQuote(data_dir),"--out",shQuote(gcmsem_out)))
  run(python,c("tutorials/run_sld4m_panel.py","--reference",shQuote(native_dir),
      "--upstream",shQuote(sld4m_dir),"--octave",shQuote(octave),"--jobs",as.integer(jobs)))
  run(python,c("tutorials/summarize_original_sld4m.py","--gcmsem",shQuote(gcmsem_out),
      "--native",shQuote(native_dir),"--upstream",shQuote(sld4m_dir),"--out",shQuote(out)))
  read_original_polygenicity(out,gcmsem_out)
}

print.gcmsem_original_polygenicity <- function(x,...) {
  cat("GCMSEM versus S-LD4M with the author's reference inputs\n")
  cat(nrow(x$comparison),"traits; S-LD4M:",x$agreement$n_annotations,"annotations,",
      x$agreement$n_blocks,"jackknife blocks\n")
  cat("Native baseline-category estimates; separate SNP selection and weights\n")
  print(x$comparison[,c("trait","gcmsem","gcmsem_se","sld4m","sld4m_se")],
      digits=5,row.names=FALSE)
  cat("\nPearson r:",format(x$agreement$pearson,digits=5),"\n")
  cat("RMSE:",format(x$agreement$rmse,digits=6),"effective variants\n")
  cat("Traits with diagnostic flags:",sum(x$comparison$flags!="none"),"of",nrow(x$comparison),"\n")
  cat("Inspect $comparison$flags and $annotations; signed estimates are retained.\n")
  invisible(x)
}

plot.gcmsem_original_polygenicity <- function(x,uncertainty=TRUE,...) {
  x$agreement$sld4m_input <- "author reference"
  plot.gcmsem_polygenicity_comparison(x,uncertainty=uncertainty,input="primary",...)
}
