# Run from the repository root after installing GCMSEM and its Python engine.
# In R: source("tutorials/height_weight.R"); hw <- height_weight_example()
# In a terminal: Rscript tutorials/height_weight.R tutorial-data

height_weight_example <- function(data_dir = "tutorial-data", moments = NULL,
                                  sumstats = NULL, n_starts = 32L,
                                  seed = 20261105L, verbose = FALSE) {
  if (!requireNamespace("GCMSEM", quietly = TRUE)) stop("Install GCMSEM first")

  # Each row names one phenotype and its SNP-level summary-statistic file.
  # These prepared BETA/SE values already use standardized genotype dosages.
  if (is.null(sumstats)) {
    sumstats <- data.frame(
      trait = c("Height", "Weight"),
      path = file.path(data_dir, "sumstats", c("50_irnt.tsv.gz", "21002_irnt.tsv.gz")),
      effect_scale = "genotype_standardized"
    )
  }
  if (is.null(moments)) {
    moments <- GCMSEM::estimate_moments(
      sumstats, scores = file.path(data_dir, "scores", "eur_hm3.scores.npz"),
      max_order = 4L, n_blocks = 200L, verbose = verbose
    )
  }
  if (!identical(moments$traits, c("Height", "Weight"))) {
    stop("Expected Height and Weight, in that order")
  }
  scaled <- GCMSEM::standardize_moments(moments)

  # The direct loadings of F are equal on the fixed SNP-effect-SD scale.
  # Its variance is freely estimated, including negative Heywood estimates.
  forward_syntax <- '
    Weight ~ forward*Height
    F =~ 1*Height + 1*Weight
    F ~~ shared_variance*F
    k3(F, F, F) = 0
    k4(F, F, F, F) = 0
    Height ~~ variance_x*Height
    Weight ~~ variance_y*Weight
    k3(Height, Height, Height) = free("skew_x")
    k3(Weight, Weight, Weight) = free("skew_y")
    k4(Height, Height, Height, Height) = free("excess_x")
    k4(Weight, Weight, Weight, Weight) = free("excess_y")
  '
  forward_model <- GCMSEM::gcmsem_model(
    forward_syntax, observed = c("Height", "Weight"),
    gaussian_observed = FALSE, gaussian_latent = TRUE
  )
  fit_control <- list(iter.max = 5000L, eval.max = 10000L, rel.tol = 1e-11)
  forward <- GCMSEM::gcmsem(
    forward_model, scaled, orders = 2:4, n_starts = n_starts, seed = seed,
    feasibility = "check", control = fit_control, verbose = verbose
  )

  # Keep the same source restrictions and reverse the directed path.
  reverse_syntax <- sub("Weight ~ forward*Height", "Height ~ reverse*Weight",
                        forward_syntax, fixed = TRUE)
  reverse_model <- GCMSEM::gcmsem_model(
    reverse_syntax, observed = c("Height", "Weight"),
    gaussian_observed = FALSE, gaussian_latent = TRUE
  )
  reverse <- GCMSEM::gcmsem(
    reverse_model, scaled, orders = 2:4, n_starts = n_starts, seed = seed + 1L,
    feasibility = "check", control = fit_control, verbose = verbose
  )

  # A separate sensitivity model estimates both directions with independent
  # trait-specific sources. stable=TRUE chooses the dynamically stable branch.
  reciprocal_model <- GCMSEM::gcmsem_model('
    Weight ~ forward*Height
    Height ~ reverse*Weight
    Height ~~ variance_x*Height
    Weight ~~ variance_y*Weight
    k3(Height, Height, Height) = free("skew_x")
    k3(Weight, Weight, Weight) = free("skew_y")
    k4(Height, Height, Height, Height) = free("excess_x")
    k4(Weight, Weight, Weight, Weight) = free("excess_y")
  ', observed = c("Height", "Weight"), gaussian_observed = FALSE)
  reciprocal <- GCMSEM::gcmsem(
    reciprocal_model, scaled, orders = 2:4, n_starts = n_starts, seed = seed + 2L,
    stable = TRUE, feasibility = "check", control = fit_control, verbose = verbose
  )

  # Transform paths back to the original phenotype units. Scaling constants
  # are held fixed when transforming both the estimate and its formal SE.
  fits <- list(forward = forward, reverse = reverse, reciprocal = reciprocal)
  divisors <- scaled$scaling$scale
  paths <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    p <- GCMSEM::parameter_estimates(fit)
    p <- p[p$parameter %in% c("forward", "reverse"), ]
    from <- ifelse(p$parameter == "forward", "Height", "Weight")
    to <- ifelse(p$parameter == "forward", "Weight", "Height")
    factor <- unname(divisors[to] / divisors[from])
    data.frame(model = name, from = from, to = to,
               effect_sd_estimate = p$estimate, effect_sd_formal_se = p$se,
               original_estimate = p$estimate * factor, original_formal_se = p$se * factor,
               se_status = p$se_status, source_admissible = fit$diagnostics$source_admissible)
  }))
  rownames(paths) <- NULL
  c(list(moments = moments, standardized = scaled), fits, list(paths = paths))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript tutorials/height_weight.R [DATA_DIR]")
  hw <- height_weight_example(if (length(args)) args[1L] else "tutorial-data")
  print(hw$moments)
  print(summary(hw$forward))
  print(hw$paths, digits = 4, row.names = FALSE)
}
