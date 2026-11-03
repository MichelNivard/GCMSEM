#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript examples/fit_bundle.R PREFIX")
library(GCMSEM)
bundle <- read_moments(args[1])
if (!identical(bundle$traits, c("A", "B", "C"))) {
  stop("This example expects the A, B, C toy or validation bundle")
}
model <- gcmsem_model('B ~ ab*A; C ~ bc*B',
                      observed = bundle$traits, gaussian_observed = FALSE)
fit <- gcmsem(model, bundle, n_starts = 8, seed = 20261103,
               log_file = paste0(args[1], ".sem.log"))
summary(fit)
write.table(parameter_estimates(fit), paste0(args[1], ".parameters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(fit, paste0(args[1], ".fit.rds"))
