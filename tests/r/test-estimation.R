# Integration tests use an independently generated identity-LD reference.
# R-only installations can still fit SEMs; skip this bridge test if the
# optional Python engine is not installed in the test environment.
library(GCMSEM)
.run_estimation_tests <- function() {
  python <- Sys.getenv("GCMSEM_TEST_PYTHON", Sys.which("python"))
  available <- nzchar(python) && requireNamespace("jsonlite", quietly = TRUE) &&
    identical(as.integer(suppressWarnings(system2(python,
      c("-c", shQuote("import numpy; import gcmsem.cli")), stdout = FALSE, stderr = FALSE))), 0L)
  if (!available) {
    cat("SKIP estimation bridge: install jsonlite and the Python GCMSEM engine to run integration tests\n")
    return(invisible(NULL))
  }
  checks <- 0L
  check <- function(x, message) {
    checks <<- checks + 1L
    if (!isTRUE(x)) stop("FAILED estimation bridge: ", message)
    cat("ok estimation", checks, message, "\n")
  }
  failure <- function(expr) tryCatch(force(expr), error = identity)
  fixture <- tempfile("gcmsem bridge space's & ")
  dir.create(fixture)
  old_wd <- getwd()
  on.exit({setwd(old_wd); unlink(fixture, recursive = TRUE)}, add = TRUE)
  script <- file.path(fixture, "make fixture.py")
  writeLines(c(
    "from pathlib import Path",
    "import csv, sys",
    "import numpy as np",
    "from gcmsem.scores import build_scores, save_scores",
    "p = Path(sys.argv[1])",
    "blocks = []",
    "for b in range(20):",
    "    n = 8",
    "    blocks.append(dict(id=f'b{b}', ld=np.eye(n), snp=np.array([f'rs{b*n+j}' for j in range(n)]), a1=np.repeat('A', n), a2=np.repeat('C', n), af=np.tile([.2, .3], n//2), chrom=np.repeat('1', n), bp=np.arange(b*1000+1, b*1000+n+1)))",
    "s = build_scores(blocks, maf_edges=(0, .25, .5))",
    "save_scores(s, p/'signed scores.npz')",
    "save_scores(build_scores(blocks, maf_edges=(0, .5)), p/'one-bin.npz')",
    "rng = np.random.default_rng(417)",
    "y = rng.normal(size=(160, 2)); y[:, 1] = .4*y[:, 0] + .8*y[:, 1]",
    "def table(name, header, rows):",
    "    with (p/name).open('w', newline='') as f:",
    "        w = csv.writer(f, delimiter='\\t'); w.writerow(header); w.writerows(rows)",
    "for k, trait in enumerate(['A', 'B']):",
    "    rows = []",
    "    for j, snp in enumerate(s['snp']):",
    "        af = s['af'][j]; sd = np.sqrt(2*af*(1-af)); flip = j % 3 == 0",
    "        rows.append([snp, 'C' if flip else 'A', 'A' if flip else 'C', (-1 if flip else 1)*y[j,k]/sd, .1/sd, 1-af if flip else af])",
    "    table(trait+' gwas.tsv', ['SNP','A1','A2','BETA','SE','EAF'], rows)",
    "table('manifest.tsv', ['trait','path','effect_scale'], [['A','A gwas.tsv','allele'],['B','B gwas.tsv','allele']])",
    "table('true.tsv', ['A','B'], y)"
  ), script)
  status <- system2(python, shQuote(c(script, fixture)),
                    stdout = file.path(fixture, "setup.log"), stderr = file.path(fixture, "setup.log"))
  if (status != 0L) stop(paste(readLines(file.path(fixture, "setup.log")), collapse = "\n"))
  scores <- file.path(fixture, "signed scores.npz")
  manifest <- file.path(fixture, "manifest.tsv")
  work <- file.path(fixture, "interactive session")
  dir.create(work); setwd(work)
  before <- list.files(all.files = TRUE, recursive = TRUE)
  temp_before <- list.files(tempdir(), pattern = "^gcmsem-estimate-", full.names = TRUE)
  zero <- matrix(0, 2, 2, dimnames = list(c("B", "A"), c("B", "A")))
  x <- estimate_moments(manifest, scores, n_blocks = 20, noise_intercept = zero,
                        python = python)
  truth <- utils::read.delim(file.path(fixture, "true.tsv"))
  expected <- stats::cov(truth) * (nrow(truth) - 1) / nrow(truth)
  check(inherits(x, "gcmsem_moments"), "returns a moment object")
  check(identical(before, list.files(all.files = TRUE, recursive = TRUE)), "default creates no files in the R working directory")
  check(identical(temp_before, list.files(tempdir(), pattern = "^gcmsem-estimate-", full.names = TRUE)), "temporary engine files are removed")
  check(is.null(x$prefix) && is.null(x$metadata$output_prefix), "temporary output paths are cleared")
  values <- stats::setNames(x$moments$estimate, x$moments$key)
  check(max(abs(values[c("A", "B")] - colMeans(truth))) < 1e-12, "identity-LD means agree with an independent R calculation")
  check(max(abs(values[c("A|A", "A|B", "B|B")] - c(expected[1,1], expected[1,2], expected[2,2]))) < 1e-12,
        "identity-LD covariance agrees with an independent R calculation")
  check(identical(dim(x$jackknife), c(20L, 14L)) && identical(colnames(x$jackknife), x$moments$key), "all deletion estimates are retained in R")
  check(identical(names(x$conditional), c("stratum-1", "stratum-2")) &&
        all(vapply(x$conditional, inherits, logical(1), "gcmsem_moments")), "annotation-specific moment objects are retained")
  check(all(vapply(x$conditional, function(z) identical(dim(z$jackknife), c(20L, 14L)) && is.null(z$prefix), logical(1))),
        "conditional deletions survive temporary cleanup")
  check(any(grepl("SUCCESS", x$log, fixed = TRUE)) && x$engine$status == 0L, "engine log and invocation remain available")
  files <- c(A = file.path(fixture, "A gwas.tsv"), B = file.path(fixture, "B gwas.tsv"))
  y <- estimate_moments(files, scores, effect_scale = c(B = "allele", A = "allele"),
    max_order = 2, n_blocks = 20, noise_intercept = zero, conditional = FALSE, python = python)
  check(max(abs(y$moments$estimate - values[y$moments$key])) < 1e-12 && !length(y$conditional),
        "named files align scales and second-order selection")
  # A virtual environment commonly has bin/python -> /usr/bin/python3.
  # Resolving that symlink would execute the wrong environment.
  venv <- file.path(fixture, "virtual environment")
  venv_status <- system2(python, c("-m", "venv", "--without-pip", "--system-site-packages", shQuote(venv)),
                        stdout = FALSE, stderr = FALSE)
  if (venv_status == 0L) {
    venv_python <- file.path(venv, if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
    original_path <- Sys.getenv("PYTHONPATH", unset = NA_character_)
    on.exit(if (is.na(original_path)) Sys.unsetenv("PYTHONPATH") else Sys.setenv(PYTHONPATH = original_path), add = TRUE)
    # Expose the already tested dependencies if the parent is itself a venv.
    search_path <- system2(python, c("-c", shQuote("import os, sys; print(os.pathsep.join(p for p in sys.path if p))")), stdout = TRUE)
    Sys.setenv(PYTHONPATH = search_path)
    from_venv <- estimate_moments(files, scores, effect_scale = "allele", max_order = 2,
      n_blocks = 20, noise_intercept = zero, conditional = FALSE, python = venv_python)
    if (is.na(original_path)) Sys.unsetenv("PYTHONPATH") else Sys.setenv(PYTHONPATH = original_path)
    check(identical(from_venv$moments, y$moments) &&
          identical(from_venv$engine$python, file.path(normalizePath(dirname(venv_python), winslash = "/"), basename(venv_python))),
          "virtual-environment executable symlinks retain their environment")
  } else cat("SKIP virtual-environment check: Python venv module is unavailable\n")
  manifest_frame <- data.frame(trait = names(files), path = unname(files), effect_scale = "allele")
  from_manifest <- estimate_moments(manifest_frame, scores, max_order = 2,
    n_blocks = 20, noise_intercept = zero, conditional = FALSE, python = python)
  check(identical(from_manifest$moments, y$moments) &&
        from_manifest$metadata$r_input$type == "manifest_data_frame",
        "an R manifest data frame returns the same numerical result")
  poly <- suppressWarnings(estimate_polygenicity(files, scores, effect_scale = "allele",
    n_blocks = 20, noise_intercept = 0, python = python))
  check(inherits(poly, "gcmsem_polygenicity") && identical(poly$traits, c("A", "B")), "standalone polygenicity returns separate trait results")
  check(identical(before, list.files(all.files = TRUE, recursive = TRUE)), "polygenicity creates no working-directory files by default")
  check(all(poly$estimates$n_tags == 160) && all(poly$estimates$n_sources == 160), "polygenicity distinguishes source and tag counts")
  for (trait in poly$traits) {
    d <- poly$details[[trait]]
    h <- sum(d$strata$source_snps*d$strata$variance)
    denominator <- sum(d$strata$source_snps*(d$strata$fourth_cumulant+3*d$strata$variance^2))
    check(abs(coef(poly)[trait] - 3*h^2/denominator) < 1e-9, "identity-LD polygenicity agrees with independent R algebra")
    check(nrow(d$jackknife) == 20 && identical(dim(d$conditional_covariance), c(4L,4L)), "deletions and joint stratum covariance retained")
  }
  check(nrow(as.data.frame(poly)) == 2 && any(grepl("intermediate estimates", capture.output(summary(poly,trait="A")))), "readable polygenicity tables and detailed summary")
  one_bin <- suppressWarnings(estimate_polygenicity(files["A"],file.path(fixture,"one-bin.npz"),
    effect_scale="allele",n_blocks=20,noise_intercept=0,python=python))
  check(identical(dim(one_bin$details$A$ld_mean),c(1L,1L)) && nrow(one_bin$details$A$strata)==1,
    "single source stratum preserves matrix and table dimensions")
  selected <- suppressWarnings(estimate_polygenicity(files, scores, traits="A", effect_scale="allele",
    n_blocks=20, noise_intercept=c(A=0), python=python))
  check(identical(coef(selected),coef(poly)["A"]), "trait selection preserves standalone analysis")
  saveRDS(poly,file.path(fixture,"poly.rds"))
  check(identical(readRDS(file.path(fixture,"poly.rds")),poly), "polygenicity RDS round trip")
  exports <- export_gcmsem(poly,file.path(fixture,"poly-export"))
  check(all(file.exists(exports)) && any(grepl("conditional_covariance",exports)), "polygenicity export includes covariance and audit tables")
  check(inherits(failure(export_gcmsem(poly,file.path(fixture,"poly-export"))),"error"), "polygenicity export protects existing files")
  grDevices::pdf(file.path(fixture,"poly.pdf")); plot(poly,trait="A"); grDevices::dev.off()
  check(file.info(file.path(fixture,"poly.pdf"))$size > 0, "polygenicity influence plot renders")
  check(inherits(failure(estimate_polygenicity(files,scores,traits="missing",effect_scale="allele")),"error"), "unknown polygenicity trait is rejected")
  check(inherits(failure(estimate_polygenicity(files,scores,noise_intercept=-1,effect_scale="allele")),"error"), "invalid noise intercept is rejected")
  tables <- lapply(files, utils::read.delim)
  z <- estimate_moments(tables, scores, effect_scale = "allele", n_blocks = 20,
    noise_intercept = zero, python = python)
  check(max(abs(z$moments$estimate - x$moments$estimate)) < 1e-12 && max(abs(z$V - x$V)) < 1e-12,
        "in-memory GWAS data frames preserve double precision and numerical results")
  check(z$metadata$r_input$type == "named_data_frames" &&
        all(grepl("R data frame: ", z$metadata$r_input$inputs$path, fixed = TRUE)), "data-frame provenance is readable after cleanup")
  third <- tables$A
  third$BETA <- tables$A$BETA + .5 * tables$B$BETA
  three <- c(list(C = third), tables)
  zero_three <- matrix(0, 3, 3, dimnames = list(c("B", "C", "A"), c("B", "C", "A")))
  three_result <- estimate_moments(three, scores, effect_scale = "allele", max_order = 2,
    n_blocks = 20, noise_intercept = zero_three, conditional = FALSE, python = python)
  three_truth <- cbind(C = truth$A + .5 * truth$B, truth)
  check(identical(dim(three_result$S), c(3L, 3L)) &&
        identical(rownames(three_result$S), c("C", "A", "B")) &&
        max(abs(three_result$R - stats::cor(three_truth))) < 1e-12,
        "three-trait matrices retain manifest order and match independent correlations")
  selected <- estimate_moments(files, scores, effect_scale = "allele", max_order = 3,
    moment_keys = "A|A|B", n_blocks = 20, noise_intercept = zero, python = python)
  check(nrow(selected$moments) == 6L && "A|A|B" %in% selected$moments$key,
        "selected higher-order tuples retain required lower orders")
  output <- file.path(work, "explicit results", "pair")
  e <- estimate_moments(files, scores, effect_scale = "allele", n_blocks = 20,
    noise_intercept = zero, python = python, write.to.tsv = TRUE, out = output)
  restored <- read_moments(output)
  check(inherits(e, "gcmsem_moments") && max(abs(restored$V - x$V)) < 1e-12 &&
        max(abs(restored$moments$estimate - x$moments$estimate)) < 1e-12,
        "explicit export also returns an R object and round-trips")
  check(file.exists(paste0(output, ".stratum-1.metadata.json")) &&
        file.exists(paste0(output, ".stratum-2.metadata.json")), "explicit export retains conditional bundles")
  old_file <- readLines(paste0(output, ".moments.tsv"))
  collision <- failure(estimate_moments(files, scores, effect_scale = "allele", n_blocks = 20,
    noise_intercept = zero, python = python, write.to.tsv = TRUE, out = output))
  check(inherits(collision, "error") && identical(old_file, readLines(paste0(output, ".moments.tsv"))),
        "explicit export refuses to overwrite an existing result")
  bad <- tables
  bad$A <- rbind(bad$A, bad$A[1L, , drop = FALSE])
  failure_out <- file.path(work, "failed output", "pair")
  err <- failure(estimate_moments(bad, scores, effect_scale = "allele", n_blocks = 20,
    noise_intercept = zero, python = python, write.to.tsv = TRUE, out = failure_out))
  check(inherits(err, "gcmsem_engine_error") && err$status == 2L &&
        any(grepl("duplicate SNP", err$log, fixed = TRUE)), "engine failure exposes status and diagnostic log in an R condition")
  check(!dir.exists(dirname(failure_out)) &&
        identical(temp_before, list.files(tempdir(), pattern = "^gcmsem-estimate-", full.names = TRUE)),
        "failed estimation leaves no partial export or temporary bundle")
  check(inherits(failure(estimate_moments(files, scores, python = python)), "error"),
        "named files require an explicit effect scale")
  check(inherits(failure(estimate_moments(files, scores, effect_scale = "allele",
    python = python, write.to.tsv = TRUE)), "error"), "explicit export requires an output prefix")
  check(inherits(failure(estimate_moments(manifest, scores, python = python,
    out = "unused")), "error"), "an output prefix requires explicit file export")
  cat("PASS estimation bridge:", checks, "checks\n")
}
.run_estimation_tests()
rm(.run_estimation_tests)
