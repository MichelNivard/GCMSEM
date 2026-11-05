# Tutorial figures

`correlation-scatterplot.png` and `correlation-scatterplot.pdf` show the same
105 distinct trait pairs from the archived 15-trait GCMSEM analysis and
standard LDSC results. The PNG is embedded in the tutorial; the PDF is a
vector version for export.

Follow sections 3–5 of
[the correlation tutorial](../Fifteen-trait-correlations.md) to read the
external result files and create the `comparison` R object. Then run:

```r
png("correlation-scatterplot.png", width = 1600, height = 1500, res = 220)
plot(comparison)
dev.off()

pdf("correlation-scatterplot.pdf", width = 1600 / 220, height = 1500 / 220)
plot(comparison)
dev.off()
```

The numerical datasets and full GWAS results are distributed separately.
