# Run from the repository root after installing the R and Python packages:
# Rscript tutorials/standalone_polygenicity.R tutorial-data tutorial-results/standalone-polygenicity
library(GCMSEM)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1L] else "tutorial-data"
result_dir <- if (length(args) > 1L) args[2L] else "tutorial-results/standalone-polygenicity"

# Start with height: a relatively precise estimate in the released panel.
height <- estimate_polygenicity(
  sumstats = c(Height = file.path(data_dir, "sumstats/50_irnt.tsv.gz")),
  scores = file.path(data_dir, "scores/eur_hm3.scores.npz"),
  effect_scale = "genotype_standardized",
  n_blocks = 200L
)
height
summary(height, trait = "Height")

# Then use the same settings for depression to inspect weak precision.
depression <- estimate_polygenicity(
  sumstats = c(Depression = file.path(data_dir, "sumstats/20544_11.tsv.gz")),
  scores = file.path(data_dir, "scores/eur_hm3.scores.npz"),
  effect_scale = "genotype_standardized",
  n_blocks = 200L
)
depression
summary(depression, trait = "Depression")

comparison <- rbind(as.data.frame(height), as.data.frame(depression))
comparison$relative_se_percent <- 100 * comparison$se / abs(comparison$effective_count)
comparison[, c("trait", "effective_count", "se", "delta_se", "relative_se_percent")]

# Explicit saves of both returned R objects and their deletion plots.
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(height, file.path(result_dir, "height.rds"))
saveRDS(depression, file.path(result_dir, "depression.rds"))
png(file.path(result_dir, "height-deletions.png"), width = 1100, height = 750, res = 140)
plot(height, trait = "Height")
invisible(dev.off())
png(file.path(result_dir, "depression-deletions.png"), width = 1100, height = 750, res = 140)
plot(depression, trait = "Depression")
invisible(dev.off())
