## Load data (genes x conditions)
data <- as.matrix(read.csv("C:/Users/geoap/Documents/WUR/msc_bioinformatics_and_systems_biology/period_3/molecular_systems_biology/msb_assignment/datasets/dataset_19.csv", row.names = 1, check.names = FALSE))

## 1) How many genes in how many conditions?
n_genes <- nrow(data)
n_conditions <- ncol(data)

cat("Genes:", n_genes, "\n")
cat("Conditions:", n_conditions, "\n")

## 2) Check for missing data points (NA / NaN / Inf)
n_na   <- sum(is.na(data))
n_nan  <- sum(is.nan(data))
n_inf  <- sum(is.infinite(data))

cat("NA:", n_na, " NaN:", n_nan, " Inf:", n_inf, "\n")

## Which genes have any missing values?
genes_with_missing <- rownames(data)[rowSums(!is.finite(data)) > 0]
cat("Genes with any non-finite values:", length(genes_with_missing), "\n")

## Which conditions have any missing values?
conds_with_missing <- colnames(data)[colSums(!is.finite(data)) > 0]
cat("Conditions with any non-finite values:", length(conds_with_missing), "\n")

## Quick peek (first few problematic entries)
if (length(genes_with_missing) > 0 || length(conds_with_missing) > 0) {
  idx <- which(!is.finite(data), arr.ind = TRUE)
  print(head(data[idx], 20))
}

## 3) Check data normalization (visual diagnostics)
## If values are already DESeq2-normalized counts, they can still be skewed:
## Use log2(x+1) for sanity-check plots.
log_data <- log2(data + 1)

## (a) Boxplots across conditions: similar medians/spreads is a good sign
par(mar=c(8,4,2,1))
boxplot(log_data, las=2, main="log2(expression+1) per condition", ylab="log2(expr+1)")

## (b) Density curves: distributions should broadly overlap
plot(density(log_data[,1]), main="Density of log2(expression+1) across conditions", xlab="log2(expr+1)")
for (j in 2:ncol(log_data)) lines(density(log_data[,j]))

## (c) Total signal per condition (not “library size” anymore, but still useful to spot outliers)
totals <- colSums(data, na.rm = TRUE)
barplot(totals, las=2, main="Column sums (sanity check)", ylab="sum(expression)")

## (d) PCA: do samples cluster sensibly? any extreme outlier?
pca <- prcomp(t(log_data), center = TRUE, scale. = FALSE)
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", pch=19)
text(pca$x[,1], pca$x[,2], labels=colnames(log_data), pos=3, cex=0.8)

## Distance between samples
dist_samples <- dist(t(log_data), method = "euclidean")

## Hierarchical clustering
hc_samples <- hclust(dist_samples, method = "complete")

## Plot dendrogram
plot(hc_samples,
     main = "Hierarchical clustering of conditions",
     xlab = "",
     sub = "")

## Log-transform for exploratory analysis
log_data <- log2(data + 1)

## Sample–sample correlation
cor_mat <- cor(log_data, method = "pearson", use = "pairwise.complete.obs")

## Inspect correlation values
round(cor_mat[1:5, 1:5], 2)

heatmap(cor_mat,
        symm = TRUE,
        margins = c(6, 6),
        main = "Correlation between conditions")

## Gene-wise variance
gene_var <- apply(log_data, 1, var)

## Top 50 most variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:50]

heatmap(log_data[top_genes, ],
        scale = "row",
        margins = c(6, 6),
        main = "Heatmap of top 50 variable genes")

dist_samples <- dist(t(log_data))
hc_samples <- hclust(dist_samples)

heatmap(cor_mat,
        symm = TRUE,
        Colv = as.dendrogram(hc_samples),
        Rowv = as.dendrogram(hc_samples),
        margins = c(6,6),
        main = "Correlation heatmap with hierarchical clustering")

log_data_imp <- log_data

for (i in 1:nrow(log_data_imp)) {
  missing <- !is.finite(log_data_imp[i, ])
  if (any(missing)) {
    log_data_imp[i, missing] <- mean(log_data_imp[i, !missing], na.rm = TRUE)
  }
}

missing_frac <- rowMeans(!is.finite(log_data))
log_data_imp <- log_data_imp[missing_frac <= 0.2, ]

cat("Genes retained after filtering:", nrow(log_data_imp), "\n")

gene_var <- apply(log_data_imp, 1, var)
log_data_imp <- log_data_imp[gene_var > 0, ]
pca <- prcomp(t(log_data_imp), center = TRUE, scale. = FALSE)

plot(pca$x[,1], pca$x[,2],
     pch = 19,
     xlab = paste0("PC1 (", round(100*summary(pca)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca)$importance[2,2],1), "%)"),
     main = "PCA of conditions (imputed data)")

text(pca$x[,1], pca$x[,2],
     labels = colnames(log_data_imp),
     pos = 3, cex = 0.8)

var_explained <- summary(pca)$importance[2, ]

plot(var_explained,
     type = "b",
     xlab = "Principal component",
     ylab = "Proportion of variance explained",
     main = "Scree plot")

cumsum(var_explained)

plot(pca$x[,1], pca$x[,3],
     xlab = paste0("PC1 (", round(100*var_explained[1],1), "%)"),
     ylab = paste0("PC3 (", round(100*var_explained[3],1), "%)"),
     pch = 19)

text(pca$x[,1], pca$x[,3], labels = colnames(log_data_imp), pos = 3, cex = 0.7)

plot(pca$x[,2], pca$x[,3],
     xlab = paste0("PC2 (", round(100*var_explained[2],1), "%)"),
     ylab = paste0("PC3 (", round(100*var_explained[3],1), "%)"),
     pch = 19)

## log_data_imp = your log-transformed + imputed matrix (genes x samples)
## If you only have data: log_data_imp <- log2(data + 1) then impute

focus_patterns <- c("^iron\\.starvation", "^stationary\\.phase", "^biofilm\\.", "^anoxic")

focus_cols <- grepl(paste(focus_patterns, collapse="|"), colnames(log_data_imp))
focus_mat  <- log_data_imp[, focus_cols]

cat("Focus samples:", ncol(focus_mat), "\n")
head(colnames(focus_mat))

cor_mat <- cor(t(focus_mat), method = "pearson", use = "pairwise.complete.obs")

library(igraph)

thresholds <- c(0.75, 0.80, 0.85, 0.90)

summ <- lapply(thresholds, function(thr){
  idx <- which(abs(cor_mat) >= thr & upper.tri(cor_mat), arr.ind = TRUE)
  edges <- data.frame(
    source = rownames(cor_mat)[idx[,1]],
    target = colnames(cor_mat)[idx[,2]],
    weight = cor_mat[idx]
  )
  g <- graph_from_data_frame(edges, directed = FALSE)
  cc <- components(g)
  data.frame(threshold=thr, nodes=vcount(g), edges=ecount(g),
             components=cc$no, largest=max(cc$csize))
})

do.call(rbind, summ)

thr_best <- 0.85  # example

idx <- which(abs(cor_mat) >= thr_best & upper.tri(cor_mat), arr.ind = TRUE)
edges <- data.frame(
  source = rownames(cor_mat)[idx[,1]],
  target = colnames(cor_mat)[idx[,2]],
  weight = cor_mat[idx]
)

write.table(edges, "PA14_focus1_edges.tsv", sep="\t", row.names=FALSE, quote=FALSE)

