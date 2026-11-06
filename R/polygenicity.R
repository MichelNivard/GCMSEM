#' Estimate effective polygenicity from GWAS summary statistics
#' @export
estimate_polygenicity <- function(sumstats, scores, effect_scale = NULL,
    traits = NULL, n_blocks = 200L, noise_intercept = NULL,
    min_info = NULL, min_n = NULL, af_tolerance = 0.1,
    allow_strand_flips = FALSE, keep_palindromic = FALSE,
    max_memory_mb = 2048, max_condition = 1e6, python = Sys.which("python"),
    write.to.tsv = FALSE, out = NULL, overwrite = FALSE, verbose = FALSE) {
  call <- match.call()
  .validate_export_options(write.to.tsv, out, overwrite)
  if (!write.to.tsv && !is.null(out)) stop("Set write.to.tsv = TRUE to use out")
  workdir <- tempfile("gcmsem-polygenicity-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
  input <- .gcmsem_prepare_sumstats(sumstats, effect_scale, workdir)
  manifest <- input$manifest
  if (!is.null(traits)) {
    .gcmsem_traits(traits)
    if (!all(traits %in% manifest$trait)) stop("traits must name rows in sumstats")
    manifest <- manifest[match(traits, manifest$trait), , drop = FALSE]
  }
  if (!is.null(noise_intercept)) {
    if (!is.numeric(noise_intercept) || any(!is.finite(noise_intercept)) || any(noise_intercept < 0))
      stop("noise_intercept must be nonnegative numeric values, or NULL to estimate it")
    if (length(noise_intercept) == 1L && is.null(names(noise_intercept)))
      noise_intercept <- setNames(rep(noise_intercept, nrow(manifest)), manifest$trait)
    if (is.null(names(noise_intercept)) || anyDuplicated(names(noise_intercept)) ||
        !setequal(names(noise_intercept), manifest$trait))
      stop("Supply one unnamed noise_intercept or a named value for each selected trait")
  }
  moments <- details <- tables <- list()
  log <- character()
  for (i in seq_len(nrow(manifest))) {
    trait <- manifest$trait[i]
    x <- .estimate_moments_engine(manifest[i, , drop = FALSE], scores,
      max_order = 4L, n_blocks = n_blocks,
      noise_intercept = if (is.null(noise_intercept)) NULL else
        matrix(noise_intercept[[trait]], 1L, 1L, dimnames = list(trait, trait)),
      min_info = min_info, min_n = min_n, af_tolerance = af_tolerance,
      allow_strand_flips = allow_strand_flips, keep_palindromic = keep_palindromic,
      max_memory_mb = max_memory_mb, max_condition = max_condition,
      python = python, verbose = verbose, .polygenicity = TRUE)
    p <- x$polygenicity
    tables[[trait]] <- as.data.frame(lapply(p$summary, function(v) if (is.null(v)) NA_real_ else v),
                                    stringsAsFactors = FALSE)
    dimnames(p$covariance) <- list(p$covariance_keys, p$covariance_keys)
    dimnames(p$conditional_covariance) <- list(p$conditional_covariance_keys, p$conditional_covariance_keys)
    dimnames(p$ld_mean) <- list(paste0("stratum-", p$strata$stratum), paste0("stratum-", p$strata$stratum))
    p$summary <- NULL
    x$polygenicity <- NULL
    x$metadata$r_input <- input$descriptor
    x$metadata$r_interface <- "estimate_polygenicity"
    moments[[trait]] <- x
    details[[trait]] <- p
    log <- c(log, paste0("TRAIT ", trait), x$log)
    if (length(p$warnings)) warning(paste0("Polygenicity ", trait, ": ",
      paste(p$warnings, collapse = "; "), ". Signed estimates and all deletions retained."), call. = FALSE)
  }
  estimates <- do.call(rbind, tables); rownames(estimates) <- NULL
  ans <- structure(list(estimates = estimates, details = details, moments = moments,
    traits = manifest$trait, log = log, call = call,
    metadata = list(estimator = "LD-adjusted Gaussian-mixture effective count",
      snp_selection = "Independent QC and eligible tags for each trait; no shared intersection or tail removal",
      se_method = "Unequal-group block jackknife of the nonlinear count",
      delta_se_method = "Full conditional variance/fourth-cumulant covariance; mean LD held fixed",
      reference_uncertainty = "Not propagated", input = input$descriptor)),
    class = "gcmsem_polygenicity")
  if (write.to.tsv) export_gcmsem(ans, out, overwrite = overwrite)
  ans
}

#' @export
print.gcmsem_polygenicity <- function(x, digits = 5L, ...) {
  digits <- .print_digits(digits)
  cat("GCMSEM effective polygenicity\n")
  cat("LD-adjusted Gaussian-mixture count; genomic-block jackknife SE\n")
  tab <- x$estimates[, c("trait", "effective_count", "se", "delta_se", "n_tags", "n_sources", "n_blocks")]
  print(tab, row.names = FALSE, digits = digits)
  flags <- x$estimates$flags != "none"
  if (any(flags)) {
    cat("*** POLYGENICITY DIAGNOSTICS: interpret flagged counts and SEs with care ***\n")
    for (i in which(flags)) cat(x$estimates$trait[i], ": ", x$estimates$flags[i], "\n", sep = "")
  }
  cat("Signed estimates retained. Reference-panel uncertainty is not included.\n")
  invisible(x)
}

#' @export
summary.gcmsem_polygenicity <- function(object, trait = NULL, digits = 5L, ...) {
  if (is.null(trait)) trait <- object$traits
  if (!is.character(trait) || !length(trait) || anyNA(trait) || anyDuplicated(trait) ||
      !all(trait %in% object$traits)) stop("trait must name one or more analysed traits")
  object$traits <- trait
  object$estimates <- object$estimates[match(trait, object$estimates$trait), , drop = FALSE]
  structure(list(object = object, digits = .print_digits(digits)), class = "summary_gcmsem_polygenicity")
}

#' @export
print.summary_gcmsem_polygenicity <- function(x, ...) {
  print(x$object, digits = x$digits)
  for (trait in x$object$traits) {
    d <- x$object$details[[trait]]
    cat("\n", trait, ": intermediate estimates and jackknife SEs\n", sep = "")
    print(d$stages, row.names = FALSE, digits = x$digits)
    cat("\nSource strata: cumulants and SEs\n")
    print(d$strata, row.names = FALSE, digits = x$digits)
    cat("\nMost influential deletions (share of estimated count variance)\n")
    b <- d$jackknife[order(d$jackknife$variance_fraction, decreasing = TRUE, na.last = TRUE), ]
    print(utils::head(b[, c("block", "chrom", "start_bp", "end_bp", "denominator", "effective_count", "variance_fraction")], 5L),
      row.names = FALSE, digits = x$digits)
    cat("\nAssumptions:\n", paste(d$assumptions, collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.gcmsem_polygenicity <- function(x, row.names = NULL, optional = FALSE, ...) {
  ans <- x$estimates
  if (!is.null(row.names)) rownames(ans) <- row.names
  ans
}

#' @export
coef.gcmsem_polygenicity <- function(object, ...) setNames(object$estimates$effective_count, object$estimates$trait)

#' @export
plot.gcmsem_polygenicity <- function(x, trait = x$traits[1L], ...) {
  if (length(trait) != 1L || is.na(trait) || !trait %in% x$traits) stop("Choose one analysed trait")
  d <- x$details[[trait]]
  b <- d$jackknife
  ok <- is.finite(b$denominator) & is.finite(b$effective_count)
  if (!any(ok)) stop("No finite deletion estimates to plot")
  graphics::plot(b$denominator[ok], b$effective_count[ok],
    xlab = "Denominator after block deletion", ylab = "Effective polygenicity after block deletion",
    main = trait, ...)
  graphics::abline(v = 0, lty = 3, col = "grey")
  graphics::points(d$stages$estimate[d$stages$quantity == "denominator"],
    x$estimates$effective_count[match(trait, x$estimates$trait)], pch = 19, col = "red")
  graphics::legend("topright", c("Block deletions", "Full data"), pch = c(1, 19), col = c("black", "red"), bty = "n")
  invisible(x)
}
