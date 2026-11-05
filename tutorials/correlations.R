# Compare the GCMSEM and standard LDSC correlations in R.
#
# source("tutorials/correlations.R") defines the functions below. It does not
# read data or write files. The worked commands are in
# tutorials/Fifteen-trait-correlations.md.

compare_correlations <- function(gcmsem, ldsc, ldsc_se = NULL,
                                 write.to.tsv = FALSE, out = NULL) {
  if (!inherits(gcmsem, "gcmsem_moments") || !is.matrix(gcmsem$R))
    stop("gcmsem must be a moment object returned by read_moments() or estimate_moments()")
  if (!is.logical(write.to.tsv) || length(write.to.tsv) != 1L || is.na(write.to.tsv))
    stop("write.to.tsv must be TRUE or FALSE")
  if (write.to.tsv && (is.null(out) || length(out) != 1L || !nzchar(out)))
    stop("Supply an output prefix in out when write.to.tsv = TRUE")

  validate_matrix <- function(x, label) {
    if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x) || nrow(x) < 2L)
      stop(label, " must be a square numeric matrix with at least two traits")
    if (is.null(rownames(x)) || is.null(colnames(x)) ||
        anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
        !setequal(rownames(x), colnames(x)))
      stop(label, " must have unique, matching row and column trait names")
    x <- x[rownames(x), rownames(x), drop = FALSE]
    if (any(!is.finite(x)) || max(abs(x - t(x))) > 1e-8)
      stop(label, " must contain finite, symmetric estimates")
    x
  }

  gcmsem_rg <- validate_matrix(gcmsem$R, "GCMSEM correlations")
  ldsc <- validate_matrix(ldsc, "LDSC correlations")
  traits <- rownames(gcmsem_rg)
  if (!setequal(traits, rownames(ldsc)))
    stop("GCMSEM and LDSC matrices must contain the same named traits")
  ldsc <- ldsc[traits, traits, drop = FALSE]
  if (!is.null(ldsc_se)) {
    ldsc_se <- validate_matrix(ldsc_se, "LDSC standard errors")
    if (!setequal(traits, rownames(ldsc_se)) || any(ldsc_se < 0))
      stop("LDSC standard errors must cover the same traits and be nonnegative")
    ldsc_se <- ldsc_se[traits, traits, drop = FALSE]
  }

  # The lower triangle contains one entry for every distinct trait pair.
  index <- which(lower.tri(gcmsem_rg), arr.ind = TRUE)
  pairs <- data.frame(
    trait1 = traits[index[, "col"]],
    trait2 = traits[index[, "row"]],
    gcmsem = gcmsem_rg[index],
    ldsc = ldsc[index],
    stringsAsFactors = FALSE
  )
  pairs$difference <- pairs$gcmsem - pairs$ldsc
  if (!is.null(ldsc_se)) pairs$ldsc_se <- ldsc_se[index]
  statistics <- data.frame(
    pairs = nrow(pairs),
    pearson = if (nrow(pairs) > 1L && stats::sd(pairs$gcmsem) > 0 &&
                  stats::sd(pairs$ldsc) > 0) stats::cor(pairs$gcmsem, pairs$ldsc) else NA_real_,
    rmse = sqrt(mean(pairs$difference^2)),
    mean_difference = mean(pairs$difference),
    max_absolute_difference = max(abs(pairs$difference))
  )
  answer <- structure(list(
    traits = traits, gcmsem = gcmsem_rg, ldsc = ldsc,
    ldsc_se = ldsc_se, pairs = pairs, statistics = statistics
  ), class = "gcmsem_correlation_comparison")

  if (write.to.tsv) {
    paths <- paste0(out, c(".pairs.tsv", ".statistics.tsv"))
    if (any(file.exists(paths))) stop("Output files already exist; choose a new out prefix")
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(pairs, paths[1L], sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(statistics, paths[2L], sep = "\t", quote = FALSE, row.names = FALSE)
  }
  answer
}

print.gcmsem_correlation_comparison <- function(x, digits = 6L, ...) {
  cat("GCMSEM and standard LDSC correlation comparison\n")
  cat(length(x$traits), "traits;", nrow(x$pairs), "distinct trait pairs\n\n")
  labels <- c("Pearson correlation", "RMSE", "Mean GCMSEM - LDSC", "Largest absolute difference")
  values <- unlist(x$statistics[c("pearson", "rmse", "mean_difference", "max_absolute_difference")])
  for (i in seq_along(labels))
    cat(sprintf("%-30s %s\n", labels[i], format(values[i], digits = digits, nsmall = 6L)))
  invisible(x)
}

plot.gcmsem_correlation_comparison <- function(x, label_largest = TRUE, ...) {
  pairs <- x$pairs
  limits <- range(c(pairs$gcmsem, pairs$ldsc))
  padding <- max(diff(limits) * 0.06, 0.03)
  limits <- limits + c(-padding, padding)
  old <- graphics::par(mar = c(5.0, 5.0, 4.4, 1.3), las = 1)
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot(pairs$ldsc, pairs$gcmsem,
    xlim = limits, ylim = limits, asp = 1, pch = 21,
    bg = grDevices::adjustcolor("#2374AB", alpha.f = 0.72),
    col = "white", cex = 1.05,
    xlab = "Standard LDSC genetic correlation",
    ylab = "GCMSEM SNP-effect correlation",
    main = sprintf("Genetic correlations across %d traits", length(x$traits)), ...)
  graphics::abline(a = 0, b = 1, lty = 2, col = "#777777", lwd = 1.2)
  graphics::mtext(sprintf("%d pairs  |  Pearson r = %.3f  |  RMSE = %.3f",
    nrow(pairs), x$statistics$pearson, x$statistics$rmse),
    side = 3, line = 0.45, cex = 0.88, col = "#444444")
  if (label_largest) {
    i <- which.max(abs(pairs$difference))
    graphics::points(pairs$ldsc[i], pairs$gcmsem[i], pch = 21,
      bg = "#D66B27", col = "white", cex = 1.35)
    label <- gsub("([a-z])([A-Z])", "\\1 \\2",
      paste(pairs$trait1[i], pairs$trait2[i], sep = " / "))
    graphics::text(pairs$ldsc[i], pairs$gcmsem[i], labels = label,
      pos = 2, offset = 0.65, cex = 0.78, col = "#A64B13")
  }
  invisible(x)
}
