#' Rescale a moment object by fixed trait scales
#' @export
standardize_moments <- function(x, scale = NULL) {
  if (!inherits(x, "gcmsem_moments")) stop("x must be a gcmsem_moments object")
  x <- .enrich_moments(x)
  traits <- x$traits
  estimated <- is.null(scale)
  if (estimated) {
    variances <- diag(x$S)
    if (any(!is.finite(variances)) || any(variances <= 0))
      stop("Automatic scaling requires positive trait variances; supply positive fixed scale values")
    scale <- stats::setNames(sqrt(variances), traits)
  }
  if (!is.numeric(scale) || length(scale) != length(traits) ||
      any(!is.finite(scale)) || any(scale <= 0))
    stop("scale must contain one positive finite number per trait")
  if (is.null(names(scale))) names(scale) <- traits
  if (anyDuplicated(names(scale)) || !setequal(names(scale), traits))
    stop("Named scale values must match the moment object's traits")
  scale <- scale[traits]
  keys <- x$moments$key
  multiplier <- vapply(strsplit(keys, "|", fixed = TRUE), function(tuple) {
    if (any(!tuple %in% traits)) stop("Moment key contains an unknown trait")
    1 / prod(scale[tuple])
  }, numeric(1))
  if (any(!is.finite(multiplier)) || any(multiplier <= 0))
    stop("Scale values produce a non-finite moment transformation")
  names(multiplier) <- keys
  out <- x
  out$moments$estimate <- x$moments$estimate * multiplier
  out$moments$se <- x$moments$se * multiplier
  out$V <- x$V[keys, keys, drop = FALSE] * outer(multiplier, multiplier)
  if (!is.null(x$jackknife)) {
    if (!is.matrix(x$jackknife) || !all(keys %in% colnames(x$jackknife)))
      stop("Jackknife columns must identify every moment key")
    out$jackknife <- sweep(x$jackknife[, keys, drop = FALSE], 2L, multiplier, "*")
  }
  if (length(x$conditional))
    out$conditional <- lapply(x$conditional, standardize_moments, scale = scale)
  previous <- x$scaling$cumulative_scale
  cumulative <- if (is.null(previous)) scale else previous[traits] * scale
  out$scaling <- list(scale = scale, cumulative_scale = cumulative,
                      estimated_from_input = estimated, held_fixed = TRUE)
  out$metadata$r_scaling <- out$scaling
  out$prefix <- NULL
  .enrich_moments(out)
}
