# The interactive R API returns complete objects and writes only on request.
local({
  library(GCMSEM)
  checks <- 0L
  check <- function(value, label) {
    checks <<- checks + 1L
    if (!isTRUE(value)) stop("Presentation check failed: ", label)
  }
  near <- function(x, y, tolerance=1e-12)
    isTRUE(all.equal(unname(x), unname(y), tolerance=tolerance, check.attributes=FALSE))
  fails <- function(expr) inherits(tryCatch(force(expr), error=identity), "error")
  work <- tempfile("gcmsem-presentation-")
  dir.create(work)
  old <- setwd(work)
  on.exit({setwd(old); unlink(work, recursive=TRUE)}, add=TRUE)

  model <- gcmsem_model("y ~ b*x", c("x", "y"), gaussian_observed=FALSE)
  truth <- model$start
  truth[] <- c(.8, .3, 1.2, .5, .1, .7, .4)
  values <- implied_moments(model, truth)
  V <- diag(.001, length(values)); dimnames(V) <- list(names(values), names(values))
  output <- capture.output(fit <- gcmsem(model, values, V, n_starts=1))
  check(inherits(fit, "gcmsem_fit") && !length(output), "quiet assigned fit returns an R object")
  check(!length(list.files(all.files=TRUE, no..=TRUE)), "default fit creates no files")
  check(is.character(fit$log) && any(grepl("Selected Q", fit$log, fixed=TRUE)), "in-memory fit retains log")
  output <- capture.output(result <- summary(fit))
  check(inherits(result, "summary_gcmsem_fit") && !length(output), "summary assignment is quiet and structured")
  check(near(result$parameters$estimate, coef(fit)), "summary preserves numeric estimates")
  check(identical(as.data.frame(fit), parameter_estimates(fit)), "data frame extractor retains complete parameter table")
  rendered <- capture.output(print(result))
  check(any(grepl("Parameter estimates", rendered, fixed=TRUE)) && any(grepl("y ~ x", rendered, fixed=TRUE)), "summary prints readable SEM terms")
  check(any(grepl("Formal sandwich SEs: available", rendered, fixed=TRUE)), "summary labels formal SE status")
  check(fails(gcmsem(model, values, V, n_starts=1, write.to.tsv=TRUE)), "explicit export requires out")
  check(!length(list.files(all.files=TRUE, no..=TRUE)), "invalid export options do not create files")

  paths <- export_gcmsem(fit, "fit/result")
  check(length(paths) == 9L && all(file.exists(paths)), "explicit fit export creates declared files")
  exported <- utils::read.delim("fit/result.parameters.tsv", check.names=FALSE)
  check(near(exported$estimate, coef(fit)), "fit parameter export retains precision")
  hashes <- tools::md5sum(paths)
  check(fails(export_gcmsem(fit, "fit/result")), "fit export refuses accidental overwrite")
  check(identical(hashes, tools::md5sum(paths)), "refused export leaves files unchanged")
  check(length(export_gcmsem(fit, "fit/result", overwrite=TRUE)) == 9L, "explicit overwrite works")
  output <- capture.output(written <- gcmsem(model, values, V, n_starts=1,
    write.to.tsv=TRUE, out="written/result"))
  check(inherits(written, "gcmsem_fit") && all(file.exists(written$output_files)) &&
    near(coef(written), coef(fit)) && !length(output), "write.to.tsv still returns the same fitted object")
  output <- capture.output(logged <- gcmsem(model, values, V, n_starts=1, log_file="explicit.log"))
  check(!length(output) && length(readLines("explicit.log")) == length(logged$log), "log_file is explicit and preserves in-memory log")
  fixed <- gcmsem_model("x ~~ 1*x", "x")
  fixed_values <- implied_moments(fixed)
  fixed_V <- diag(.01, length(fixed_values)); dimnames(fixed_V) <- list(names(fixed_values), names(fixed_values))
  fixed_fit <- gcmsem(fixed, fixed_values, fixed_V, n_starts=1)
  check(all(file.exists(export_gcmsem(fixed_fit, "fixed/result"))), "fixed model export supports zero free parameters")
  check(inherits(summary(fixed_fit), "summary_gcmsem_fit"), "fixed model summary remains structured")

  keys <- c("A", "B", "C", "A|A", "A|B", "A|C", "B|B", "B|C", "C|C")
  estimates <- c(0, 0, 0, 4, 3, 1, 9, 2, 16)
  sampling <- diag(.01, length(keys)); dimnames(sampling) <- list(keys, keys)
  sampling["A|B", "A|A"] <- sampling["A|A", "A|B"] <- .002
  sampling["A|B", "B|B"] <- sampling["B|B", "A|B"] <- .003
  sampling["A|A", "B|B"] <- sampling["B|B", "A|A"] <- -.001
  moments <- GCMSEM:::.enrich_moments(structure(list(
    moments=data.frame(key=keys, order=lengths(strsplit(keys, "|", fixed=TRUE)),
      estimate=estimates, se=sqrt(diag(sampling))),
    V=sampling, traits=c("A", "B", "C"),
    metadata=list(schema_version="1.0", moment_type="cumulant", traits=c("A", "B", "C"),
      n_tags=1234, jackknife_blocks=data.frame(block=c(0L, 1L))),
    jackknife=matrix(rep(estimates, each=2L), 2L, dimnames=list(c("0", "1"), keys)),
    log=c("example log line", "complete")), class="gcmsem_moments"))
  check(identical(dim(moments$S), c(3L, 3L)) && identical(rownames(moments$R), moments$traits), "N-trait matrices retain dimensions and trait names")
  check(moments$R["A", "B"] == .5, "correlation derives from covariance")
  g <- c(1/6, -.5/8, -.5/18)
  k <- c("A|B", "A|A", "B|B")
  expected <- sqrt(sum(g^2 * diag(sampling[k, k])) +
    2*g[1]*g[2]*.002 + 2*g[1]*g[3]*.003 - 2*g[2]*g[3]*.001)
  check(near(moments$R_SE["A", "B"], expected), "correlation SE uses the full delta-method sampling covariance")
  check(identical(vcov(moments), sampling) && near(coef(moments), estimates), "moment coef and vcov extractors preserve values")
  check(identical(as.data.frame(moments), moments$moments), "moment data frame extractor returns complete table")
  output <- capture.output(ms <- summary(moments))
  check(inherits(ms, "summary_gcmsem_moments") && !length(output) && identical(ms$V, sampling), "moment summary is quiet and retains full matrices")
  rendered <- capture.output(print(moments, max.traits=2L))
  check(any(grepl("Showing 2 of 3 traits", rendered, fixed=TRUE)), "compact matrix preview identifies complete object matrix")
  check(any(grepl("delta method using full V", capture.output(print(ms)), fixed=TRUE)), "correlation uncertainty method is printed explicitly")
  unrestricted <- moments
  unrestricted$moments$estimate[unrestricted$moments$key == "A|B"] <- 12
  unrestricted <- GCMSEM:::.enrich_moments(unrestricted)
  check(unrestricted$R["A", "B"] == 2, "estimated correlations are not clipped")
  unrestricted$moments$estimate[unrestricted$moments$key == "A|A"] <- -4
  unrestricted <- GCMSEM:::.enrich_moments(unrestricted)
  check(unrestricted$S["A", "A"] == -4 && all(is.na(unrestricted$R["A", ])), "nonpositive variance retained with undefined correlations")

  # Structured metadata export requires the optional jsonlite dependency.
  if (requireNamespace("jsonlite", quietly=TRUE)) {
    paths <- export_gcmsem(moments, "moments/result")
    reread <- read_moments("moments/result")
    check(inherits(reread, "gcmsem_moments") && near(reread$moments$estimate, estimates) &&
      near(reread$V, sampling) && near(reread$R_SE, moments$R_SE), "moment bundle numerical round trip")
    check(near(reread$jackknife, moments$jackknife) && identical(dimnames(reread$jackknife), dimnames(moments$jackknife)) &&
      identical(reread$log, moments$log), "moment export preserves jackknife draws and logs")
    check(reread$metadata$n_tags == 1234 && reread$diagnostics$n_blocks == 2L, "moment export preserves sample counts")
    moments$conditional <- list(`stratum-1`=moments)
    paths <- export_gcmsem(moments, "conditional/result")
    check(all(file.exists(paths)) && file.exists("conditional/result.stratum-1.moments.tsv"), "conditional moment bundles are exported together")
    dir.create("collision")
    writeLines("existing child", "collision/result.stratum-1.moments.tsv")
    check(fails(export_gcmsem(moments, "collision/result")) &&
      !file.exists("collision/result.moments.tsv"), "child collision is detected before any parent write")
  }

  model_bad <- gcmsem_model("x ~~ v*x", "x", gaussian_observed=FALSE)
  target <- c("x|x"=-1, "x|x|x"=0, "x|x|x|x"=-4)
  bad_V <- diag(.1, 3L); dimnames(bad_V) <- list(names(target), names(target))
  warnings <- character()
  bad_fit <- withCallingHandlers(gcmsem(model_bad, target, bad_V, n_starts=1),
    warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
  check(any(grepl("CLASSICAL HEYWOOD CASE", warnings, fixed=TRUE)), "quiet default preserves blaring Heywood warning")
  check(coef(bad_fit)["v"] < 0 && all(is.finite(parameter_estimates(bad_fit)$se)), "unconstrained estimates and computable formal SEs remain available")
  rendered <- capture.output(print(summary(bad_fit)))
  check(any(grepl("*** HEYWOOD CASE", rendered, fixed=TRUE)) &&
    any(grepl("formal_only_flagged", rendered, fixed=TRUE)), "printed summary retains source warning and uncertainty flags")
  cat("PASS presentation:", checks, "checks\n")
})
