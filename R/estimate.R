# R entry point to the Python summary-statistics engine. Engine interchange
# files are an implementation detail unless the caller requests an export.

.gcmsem_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  x
}

.gcmsem_number <- function(x, name, minimum = 0, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      (integer && x != floor(x)))
    stop(name, " must be a finite ", if (integer) "integer " else "number ",
         "at least ", minimum, call. = FALSE)
  x
}

.gcmsem_existing_file <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      !file.exists(x) || dir.exists(x))
    stop(name, " must name an existing file", call. = FALSE)
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

.gcmsem_traits <- function(x) {
  if (!is.character(x) || !length(x) || anyNA(x) || anyDuplicated(x) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", x)))
    stop("Trait names must be unique, start with a letter, and contain letters, digits, _ or .", call. = FALSE)
  x
}

.gcmsem_effect_scales <- function(x, traits) {
  if (is.null(x))
    stop("Supply effect_scale = 'allele' or 'genotype_standardized' for named files or R data frames", call. = FALSE)
  if (!is.character(x) || anyNA(x)) stop("effect_scale must be character", call. = FALSE)
  if (length(x) == 1L && is.null(names(x))) x <- rep(x, length(traits)) else {
    if (is.null(names(x)) || anyDuplicated(names(x)) || !setequal(names(x), traits))
      stop("Use one effect_scale for all traits, or a named vector with one value per trait", call. = FALSE)
    x <- x[traits]
  }
  if (any(!x %in% c("allele", "genotype_standardized")))
    stop("effect_scale values must be 'allele' or 'genotype_standardized'", call. = FALSE)
  unname(x)
}

.gcmsem_write_input_table <- function(x, path) {
  # Seventeen significant digits preserve double precision across the bridge.
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_along(x)) {
    if (is.factor(x[[i]])) x[[i]] <- as.character(x[[i]])
    if (is.numeric(x[[i]]) && !is.integer(x[[i]]))
      x[[i]] <- format(x[[i]], digits = 17L, scientific = TRUE, trim = TRUE)
  }
  utils::write.table(x, path, sep = "\t", quote = TRUE, qmethod = "double",
                     row.names = FALSE, na = "NA")
}

.gcmsem_prepare_sumstats <- function(sumstats, effect_scale, workdir) {
  source <- NULL
  if (is.character(sumstats) && length(sumstats) == 1L && is.null(names(sumstats))) {
    if (!is.null(effect_scale)) stop("The manifest supplies effect_scale; omit the separate effect_scale argument", call. = FALSE)
    source <- .gcmsem_existing_file(sumstats, "sumstats")
    con <- if (grepl("\\.gz$", source)) gzfile(source, "rt") else file(source, "rt")
    first <- readLines(con, n = 1L, warn = FALSE); close(con)
    if (!length(first)) stop("The summary-statistics manifest is empty", call. = FALSE)
    sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(",", first, fixed = TRUE)) "," else ""
    manifest <- utils::read.table(source, header = TRUE, sep = sep, quote = "\"",
      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
      colClasses = "character")
    kind <- "manifest_file"
  } else if (is.data.frame(sumstats)) {
    if (!is.null(effect_scale)) stop("The manifest supplies effect_scale; omit the separate effect_scale argument", call. = FALSE)
    manifest <- sumstats
    kind <- "manifest_data_frame"
  } else if (is.character(sumstats) && !is.null(names(sumstats))) {
    traits <- .gcmsem_traits(names(sumstats))
    manifest <- data.frame(trait = traits, path = unname(sumstats),
      effect_scale = .gcmsem_effect_scales(effect_scale, traits), stringsAsFactors = FALSE)
    kind <- "named_files"
  } else if (is.list(sumstats) && !is.null(names(sumstats)) &&
             all(vapply(sumstats, is.data.frame, logical(1)))) {
    traits <- .gcmsem_traits(names(sumstats))
    scale <- .gcmsem_effect_scales(effect_scale, traits)
    paths <- file.path(workdir, paste0("sumstats-", seq_along(sumstats), ".tsv"))
    for (i in seq_along(sumstats)) {
      if (!nrow(sumstats[[i]]) || !ncol(sumstats[[i]]) ||
          anyDuplicated(tolower(names(sumstats[[i]]))))
        stop("Each GWAS data frame must have rows and uniquely named columns", call. = FALSE)
      .gcmsem_write_input_table(sumstats[[i]], paths[i])
    }
    manifest <- data.frame(trait = traits, path = paths, effect_scale = scale,
                           stringsAsFactors = FALSE)
    kind <- "named_data_frames"
  } else stop("sumstats must be a manifest path/data.frame, a named vector of paths, or a named list of GWAS data.frames", call. = FALSE)
  names(manifest) <- tolower(trimws(names(manifest)))
  if (anyDuplicated(names(manifest)) ||
      !all(c("trait", "path", "effect_scale") %in% names(manifest)) || !nrow(manifest))
    stop("The manifest requires trait, path and effect_scale columns and at least one row", call. = FALSE)
  for (nm in c("trait", "path", "effect_scale")) manifest[[nm]] <- trimws(as.character(manifest[[nm]]))
  .gcmsem_traits(manifest$trait)
  if (anyNA(manifest$effect_scale) || any(!manifest$effect_scale %in% c("allele", "genotype_standardized")))
    stop("Manifest effect_scale values must be 'allele' or 'genotype_standardized'", call. = FALSE)
  if (anyNA(manifest$path) || any(!nzchar(manifest$path))) stop("Manifest paths must be nonempty", call. = FALSE)
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", manifest$path)
  if (!is.null(source)) manifest$path[!absolute] <- file.path(dirname(source), manifest$path[!absolute])
  manifest$path <- vapply(manifest$path, .gcmsem_existing_file, "", name = "GWAS path")
  path <- source
  if (is.null(path)) {
    path <- file.path(workdir, "sumstats.tsv")
    .gcmsem_write_input_table(manifest, path)
  }
  descriptor <- manifest
  if (kind == "named_data_frames") descriptor$path <- paste0("R data frame: ", manifest$trait)
  list(path = path, manifest = manifest,
       descriptor = list(type = kind, manifest = source, inputs = descriptor))
}

.gcmsem_read_engine_bundle <- function(prefix, temporary = TRUE) {
  x <- read_moments(prefix)
  jk <- utils::read.delim(paste0(prefix, ".jackknife.tsv"), check.names = FALSE,
                         stringsAsFactors = FALSE)
  if (!identical(names(jk), c("block", x$moments$key)) || anyDuplicated(jk$block))
    stop("Invalid jackknife output from the Python engine", call. = FALSE)
  x$jackknife <- as.matrix(jk[, -1L, drop = FALSE])
  rownames(x$jackknife) <- as.character(jk$block)
  if (!is.numeric(x$jackknife) || any(!is.finite(x$jackknife)))
    stop("Non-finite jackknife output from the Python engine", call. = FALSE)
  if (temporary) {
    x$prefix <- NULL
    x$metadata$output_prefix <- NULL
    x$metadata$output_files <- NULL
  }
  x
}

#' Estimate genomic moments from GWAS summary statistics
#' @export
estimate_moments <- function(sumstats, scores, effect_scale = NULL,
    max_order = 4L, n_blocks = 200L, moment_keys = NULL, noise_intercept = NULL,
    min_info = NULL, min_n = NULL, af_tolerance = 0.1,
    allow_strand_flips = FALSE, keep_palindromic = FALSE,
    max_moments = 5000L, max_memory_mb = 2048, max_condition = 1e6,
    conditional = TRUE, python = Sys.which("python"),
    write.to.tsv = FALSE, out = NULL, overwrite = FALSE, verbose = FALSE) {
  ans <- do.call(.estimate_moments_engine, as.list(environment()))
  ans$call <- match.call()
  ans
}

.estimate_moments_engine <- function(sumstats, scores, effect_scale = NULL,
    max_order = 4L, n_blocks = 200L, moment_keys = NULL, noise_intercept = NULL,
    min_info = NULL, min_n = NULL, af_tolerance = 0.1,
    allow_strand_flips = FALSE, keep_palindromic = FALSE,
    max_moments = 5000L, max_memory_mb = 2048, max_condition = 1e6,
    conditional = TRUE, python = Sys.which("python"),
    write.to.tsv = FALSE, out = NULL, overwrite = FALSE, verbose = FALSE, .polygenicity = FALSE) {
  call <- match.call()
  for (nm in c("allow_strand_flips", "keep_palindromic", "conditional",
               "write.to.tsv", "overwrite", "verbose"))
    .gcmsem_scalar_flag(get(nm), nm)
  .gcmsem_number(max_order, "max_order", 2, TRUE)
  if (!max_order %in% 2:4) stop("max_order must be 2, 3 or 4", call. = FALSE)
  .gcmsem_number(n_blocks, "n_blocks", 2, TRUE)
  .gcmsem_number(max_moments, "max_moments", 1, TRUE)
  .gcmsem_number(max_memory_mb, "max_memory_mb", .Machine$double.eps)
  .gcmsem_number(max_condition, "max_condition", 1)
  .gcmsem_number(af_tolerance, "af_tolerance", 0)
  if (af_tolerance >= .5) stop("af_tolerance must be less than 0.5", call. = FALSE)
  if (!is.null(min_info)) .gcmsem_number(min_info, "min_info", 0)
  if (!is.null(min_n)) .gcmsem_number(min_n, "min_n", 0)
  if (write.to.tsv && (is.null(out) || !is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out)))
    stop("Supply an output prefix in out when write.to.tsv = TRUE", call. = FALSE)
  if (!write.to.tsv && !is.null(out))
    stop("Set write.to.tsv = TRUE to use out, or omit out to return an R object only", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("estimate_moments() requires jsonlite; install it with install.packages('jsonlite')", call. = FALSE)
  scores <- .gcmsem_existing_file(scores, "scores")
  if (!is.character(python) || length(python) != 1L || is.na(python) || !nzchar(python))
    stop("Supply python as the path to a Python executable with the gcmsem engine installed", call. = FALSE)
  if (!file.exists(python)) python <- unname(Sys.which(python))
  if (!nzchar(python) || !file.exists(python) || dir.exists(python))
    stop("python must name an existing executable", call. = FALSE)
  # Keep the executable symlink: resolving .venv/bin/python to the base
  # interpreter would silently bypass the user's virtual environment.
  python <- path.expand(python)
  python <- file.path(normalizePath(dirname(python), winslash = "/", mustWork = TRUE), basename(python))
  if (file.access(python, 1L) != 0L) stop("python is not executable", call. = FALSE)
  workdir <- tempfile("gcmsem-estimate-")
  if (!dir.create(workdir)) stop("Could not create a temporary working directory", call. = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  input <- .gcmsem_prepare_sumstats(sumstats, effect_scale, workdir)
  traits <- input$manifest$trait
  prefix <- file.path(workdir, "moments")
  args <- c("-m", "gcmsem.cli", if (.polygenicity) "polygenicity" else "estimate", "--sumstats", input$path,
    "--scores", scores, "--out", prefix, "--max-order", as.character(max_order),
    "--blocks", as.character(n_blocks), "--af-tolerance", as.character(af_tolerance),
    "--max-moments", as.character(max_moments), "--max-memory-mb", as.character(max_memory_mb),
    "--max-condition", as.character(max_condition))
  for (pair in list(c("min_info", "--min-info"), c("min_n", "--min-n")))
    if (!is.null(get(pair[1L]))) args <- c(args, pair[2L], as.character(get(pair[1L])))
  if (allow_strand_flips) args <- c(args, "--allow-strand-flips")
  if (keep_palindromic) args <- c(args, "--keep-palindromic")
  if (!conditional) args <- c(args, "--no-conditional")
  if (!is.null(moment_keys)) {
    if (!is.character(moment_keys) || !length(moment_keys) || anyNA(moment_keys) ||
        any(!nzchar(moment_keys)) || any(grepl("[\r\n]", moment_keys)))
      stop("moment_keys must be a character vector such as 'Height|Height|Weight'", call. = FALSE)
    selection <- file.path(workdir, "selected-moments.txt")
    writeLines(moment_keys, selection)
    args <- c(args, "--moments", selection)
  }
  if (!is.null(noise_intercept)) {
    if (is.character(noise_intercept)) noise_path <- .gcmsem_existing_file(noise_intercept, "noise_intercept") else {
      z <- as.matrix(noise_intercept)
      if (!is.numeric(z) || !identical(dim(z), c(length(traits), length(traits))) || any(!is.finite(z)))
        stop("noise_intercept must be a finite traits-by-traits numeric matrix or a labeled TSV path", call. = FALSE)
      if (is.null(rownames(z)) && is.null(colnames(z))) dimnames(z) <- list(traits, traits)
      if (is.null(rownames(z)) || is.null(colnames(z)) || anyDuplicated(rownames(z)) || anyDuplicated(colnames(z)) ||
          !setequal(rownames(z), traits) || !setequal(colnames(z), traits))
        stop("noise_intercept row and column names must match the manifest traits", call. = FALSE)
      z <- z[traits, traits, drop = FALSE]
      noise_path <- file.path(workdir, "noise-intercept.tsv")
      .gcmsem_write_input_table(data.frame(trait = traits, z, check.names = FALSE), noise_path)
    }
    args <- c(args, "--noise-intercept", noise_path)
  }
  console <- file.path(workdir, "engine-console.log")
  if (verbose) cat("GCMSEM: estimating moments for ", paste(traits, collapse = ", "), "\n", sep = "")
  status <- tryCatch(suppressWarnings(system2(python, args = shQuote(args),
    stdout = console, stderr = console)), error = function(e) {
      structure(127L, startup_message = conditionMessage(e))
    })
  log <- if (file.exists(console)) readLines(console, warn = FALSE) else character()
  if (verbose && length(log)) cat(paste(log, collapse = "\n"), "\n")
  if (!identical(as.integer(status), 0L) || !file.exists(paste0(prefix, ".metadata.json"))) {
    detail <- if (length(log)) paste(utils::tail(log, 12L), collapse = "\n") else attr(status, "startup_message")
    if (is.null(detail)) detail <- "The engine did not produce a completed moment bundle."
    hint <- if (as.integer(status) == 127L || any(grepl("No module named .*gcmsem", log)))
      " Use a Python environment with the gcmsem engine installed." else ""
    stop(structure(list(message = paste0("GCMSEM Python estimation failed (exit status ", status,
      ").", hint, "\n", detail), call = call,
      status = as.integer(status), log = log), class = c("gcmsem_engine_error", "error", "condition")))
  }
  result <- .gcmsem_read_engine_bundle(prefix)
  if (.polygenicity) {
    polyfile <- paste0(prefix, ".polygenicity.json")
    if (!file.exists(polyfile)) stop("Update the Python GCMSEM engine: missing polygenicity output")
    result$polygenicity <- jsonlite::fromJSON(polyfile, simplifyMatrix = TRUE)
  }
  conditional_prefixes <- result$metadata$conditional_prefixes
  result$conditional <- list()
  if (length(conditional_prefixes)) {
    if (any(!grepl("^moments\\.stratum-[1-9][0-9]*$", conditional_prefixes)))
      stop("Unexpected conditional bundle names from the Python engine", call. = FALSE)
    for (cp in conditional_prefixes)
      result$conditional[[sub("^moments\\.", "", cp)]] <- .gcmsem_read_engine_bundle(file.path(workdir, cp))
  }
  result$metadata$conditional_prefixes <- NULL
  result$metadata$r_input <- input$descriptor
  result$metadata$r_interface <- "estimate_moments"
  result$log <- log
  result$engine <- list(python = python, command = c(python, args), status = as.integer(status))
  result$call <- call
  if (write.to.tsv) {
    export_gcmsem(result, out, overwrite = overwrite)
    result$prefix <- normalizePath(out, winslash = "/", mustWork = FALSE)
  }
  if (length(result$metadata$warnings))
    warning(paste(unique(unlist(result$metadata$warnings)), collapse = "\n"), call. = FALSE)
  result
}
