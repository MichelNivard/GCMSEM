library(GCMSEM)
keys <- c("X", "Y", "X|X", "X|Y", "Y|Y", "X|X|Y", "X|Y|Y|Y")
point <- c(0.2, -0.1, 4, 3, 9, 0.8, -2)
V <- diag(seq_along(keys) / 100)
dimnames(V) <- list(keys, keys)
x <- structure(list(moments = data.frame(key = keys,
  order = lengths(strsplit(keys, "|", fixed = TRUE)), estimate = point,
  se = sqrt(diag(V))), V = V, traits = c("X", "Y"), metadata = list(),
  jackknife = rbind(b1 = stats::setNames(point + 0.01, keys),
                    b2 = stats::setNames(point - 0.01, keys))), class = "gcmsem_moments")
before <- serialize(x, NULL)
z <- standardize_moments(x)
f <- c(1/2, 1/3, 1/4, 1/6, 1/9, 1/12, 1/54)
stopifnot(identical(serialize(x, NULL), before),
  isTRUE(all.equal(z$moments$estimate, point * f, check.attributes = FALSE)),
  isTRUE(all.equal(z$V, V * outer(f, f))),
  isTRUE(all.equal(diag(z$S), c(X = 1, Y = 1))),
  isTRUE(all.equal(unname(z$S["X", "Y"]), 0.5)),
  isTRUE(all.equal(z$jackknife, sweep(x$jackknife, 2L, f, "*"))),
  isTRUE(all.equal(z$scaling$scale, c(X = 2, Y = 3))))
z2 <- standardize_moments(z, c(Y = 2, X = 4))
stopifnot(isTRUE(all.equal(z2$scaling$cumulative_scale, c(X = 8, Y = 6))))
x$conditional <- list(stratum.1 = x)
zc <- standardize_moments(x, c(X = 2, Y = 3))
stopifnot(isTRUE(all.equal(zc$conditional[[1L]]$moments$estimate, z$moments$estimate)))
bad <- try(standardize_moments(x, c(X = 0, Y = 3)), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
if (requireNamespace("jsonlite", quietly = TRUE)) local({
  work <- tempfile("gcmsem-transform-export-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prefix <- file.path(work, "scaled")
  export_gcmsem(z, prefix)
  restored <- read_moments(prefix)
  stopifnot(isTRUE(all.equal(restored$scaling, z$scaling)),
    isTRUE(all.equal(standardize_moments(restored, c(Y=2, X=4))$scaling$cumulative_scale,
                     c(X=8, Y=6))))
  metadata <- jsonlite::read_json(paste0(prefix, ".metadata.json"))
  stopifnot(identical(names(metadata$r_scaling$scale), c("X", "Y")))
  # Named scales retain their identity independently of JSON member order.
  metadata$r_scaling$scale <- metadata$r_scaling$scale[c("Y", "X")]
  metadata$r_scaling$cumulative_scale <- metadata$r_scaling$cumulative_scale[c("Y", "X")]
  jsonlite::write_json(metadata, paste0(prefix, ".metadata.json"), auto_unbox=TRUE)
  stopifnot(isTRUE(all.equal(read_moments(prefix)$scaling, z$scaling)))
  # Earlier JSON arrays are interpreted in the documented trait order.
  metadata$r_scaling$scale <- c(2, 3)
  metadata$r_scaling$cumulative_scale <- c(2, 3)
  jsonlite::write_json(metadata, paste0(prefix, ".metadata.json"), auto_unbox=TRUE)
  stopifnot(isTRUE(all.equal(read_moments(prefix)$scaling, z$scaling)))
  metadata$r_scaling$scale <- list(X=2, unknown=3)
  jsonlite::write_json(metadata, paste0(prefix, ".metadata.json"), auto_unbox=TRUE)
  stopifnot(inherits(try(read_moments(prefix), silent=TRUE), "try-error"))
  # Raw metadata can be preserved while incorporating a new R transformation.
  raw <- z
  raw$metadata <- list(raw_json='{"schema_version":"1.0","moment_type":"cumulant","n_tags":123}',
                       r_scaling=z$scaling)
  raw_prefix <- file.path(work, "raw")
  export_gcmsem(raw, raw_prefix)
  raw_restored <- read_moments(raw_prefix)
  stopifnot(raw_restored$metadata$n_tags == 123,
    isTRUE(all.equal(raw_restored$scaling, z$scaling)))
  # Replacing an object must not mix current estimates with previous draws/logs.
  optional_prefix <- file.path(work, "optional")
  previous <- z; previous$log <- "Previous analysis"
  export_gcmsem(previous, optional_prefix)
  replacement <- z; replacement$jackknife <- NULL; replacement$log <- NULL
  replacement$moments$estimate <- replacement$moments$estimate * 2
  export_gcmsem(replacement, optional_prefix, overwrite=TRUE)
  reread <- read_moments(optional_prefix)
  stopifnot(is.null(reread$jackknife), !length(reread$log),
    !file.exists(paste0(optional_prefix, ".jackknife.tsv")),
    !file.exists(paste0(optional_prefix, ".log")),
    isTRUE(all.equal(reread$moments$estimate, replacement$moments$estimate)))
  # An isolated optional member is also an export collision.
  orphan_prefix <- file.path(work, "orphan")
  writeLines("existing", paste0(orphan_prefix, ".log"))
  stopifnot(inherits(try(export_gcmsem(replacement, orphan_prefix), silent=TRUE), "try-error"),
    !file.exists(paste0(orphan_prefix, ".moments.tsv")))
})
cat("Moment transformation checks passed\n")
