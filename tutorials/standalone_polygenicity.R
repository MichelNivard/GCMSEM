# Run from the repository root after installing the R and Python packages:
# Rscript tutorials/standalone_polygenicity.R tutorial-data tutorial-results/standalone-polygenicity
library(GCMSEM)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1L] else "tutorial-data"
result_dir <- if (length(args) > 1L) args[2L] else "tutorial-results/standalone-polygenicity"

p <- estimate_polygenicity(
  sumstats = c(Depression = file.path(data_dir, "sumstats/20544_11.tsv.gz")),
  scores = file.path(data_dir, "scores/eur_hm3.scores.npz"),
  effect_scale = "genotype_standardized",
  n_blocks = 200L
)
p
p$details$Depression$stages
summary(p, trait = "Depression")

# These are explicit saves of the returned R object and an illustrative plot.
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(p, file.path(result_dir, "depression.rds"))
png(file.path(result_dir, "depression-deletions.png"), width = 1100, height = 750, res = 140)
plot(p, trait = "Depression")
invisible(dev.off())
