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

## 0) Make sure focus_mat exists and has focus samples
stopifnot(exists("focus_mat"))
stopifnot(ncol(focus_mat) > 1)

## 1) Filter zero-variance genes (prevents weird behavior)
gene_var_focus <- apply(focus_mat, 1, var)
focus_mat_filt <- focus_mat[gene_var_focus > 0, , drop = FALSE]

cat("Focus samples:", ncol(focus_mat_filt), "\n")
cat("Genes kept for network:", nrow(focus_mat_filt), "\n")

stopifnot(nrow(focus_mat_filt) > 2)

## 2) Create gene–gene correlation matrix (THIS defines cor_gene)
cor_gene <- cor(
  t(focus_mat_filt),
  method = "pearson",
  use = "pairwise.complete.obs"
)

diag(cor_gene) <- 0

## 3) Confirm it exists (so you never see that error again)
cat("cor_gene created:", exists("cor_gene"), "\n")
cat("cor_gene dim:", paste(dim(cor_gene), collapse = " x "), "\n")
ls()[grepl("cor_gene", ls())]


library(igraph)

## Choose a few candidate sparsity levels (top % of absolute correlations)
top_props <- c(0.001, 0.002, 0.005, 0.01)  # 0.1%, 0.2%, 0.5%, 1%

get_edges_by_top_prop <- function(cor_mat, prop){
  w <- abs(cor_mat[upper.tri(cor_mat)])
  thr <- quantile(w, probs = 1 - prop, na.rm = TRUE)
  idx <- which(abs(cor_mat) >= thr & upper.tri(cor_mat), arr.ind = TRUE)
  edges <- data.frame(
    source = rownames(cor_mat)[idx[,1]],
    target = colnames(cor_mat)[idx[,2]],
    weight = cor_mat[idx],
    abs_weight = abs(cor_mat[idx])
  )
  list(edges = edges, thr = as.numeric(thr))
}

summ_top <- lapply(top_props, function(p){
  res <- get_edges_by_top_prop(cor_gene, p)
  g <- graph_from_data_frame(res$edges, directed = FALSE)
  cc <- components(g)
  data.frame(
    strategy = "top_prop",
    top_prop = p,
    implied_abs_r_threshold = res$thr,
    nodes = vcount(g),
    edges = ecount(g),
    components = cc$no,
    largest_component = max(cc$csize),
    density = edge_density(g, loops = FALSE)
  )
})

summ_top_df <- do.call(rbind, summ_top)
summ_top_df

plot(summ_top_df$implied_abs_r_threshold, summ_top_df$largest_component,
     xlab="Implied |correlation| threshold",
     ylab="Largest component size",
     main="Threshold choice: connectivity vs strictness",
     pch=19)

plot(summ_top_df$implied_abs_r_threshold, summ_top_df$edges,
     xlab="Implied |correlation| threshold",
     ylab="Number of edges",
     main="Threshold choice: edges vs strictness",
     pch=19)

## choose one row (example: top 0.5% strongest edges)
chosen_prop <- 0.005

res <- get_edges_by_top_prop(cor_gene, chosen_prop)
edges <- res$edges
thr_best <- res$thr

cat("Chosen top proportion:", chosen_prop, "\n")
cat("Implied |r| threshold:", thr_best, "\n")
cat("Edges:", nrow(edges), " Nodes:", length(unique(c(edges$source, edges$target))), "\n")

write.table(edges[,c("source","target","weight")],
            "PA14_focus1_edges.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

library(igraph)

## Try less strict sparsity levels
top_props <- c(0.02, 0.05, 0.10)   # 2%, 5%, 10%

summ_top2 <- lapply(top_props, function(p){
  res <- get_edges_by_top_prop(cor_gene, p)
  g <- graph_from_data_frame(res$edges, directed = FALSE)
  cc <- components(g)
  data.frame(
    top_prop = p,
    implied_abs_r_threshold = res$thr,
    nodes = vcount(g),
    edges = ecount(g),
    components = cc$no,
    largest_component = max(cc$csize),
    density = edge_density(g, loops = FALSE)
  )
})

summ_top2_df <- do.call(rbind, summ_top2)
summ_top2_df

chosen_prop <- 0.05

res <- get_edges_by_top_prop(cor_gene, chosen_prop)
edges <- res$edges
thr_best <- res$thr

cat("Chosen top proportion:", chosen_prop, "\n")
cat("Implied |r| threshold:", thr_best, "\n")
cat("Nodes:", length(unique(c(edges$source, edges$target))), 
    "Edges:", nrow(edges), "\n")

write.table(edges[,c("source","target","weight")],
            "PA14_focus1_edges_5pct.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

## Count replicates per condition (strip trailing _number)
cond_name <- sub("_[0-9]+$", "", colnames(log_data_imp))
rep_counts <- sort(table(cond_name), decreasing = TRUE)
rep_counts

conditions_4 <- c("iron.starvation", "stationary.phase", "biofilm.48h", "transitional.growth")

library(igraph)

build_condition_network <- function(expr_mat, condition_prefix,
                                    top_prop = 0.05,
                                    cor_method = "pearson",
                                    out_dir = ".",
                                    min_reps_warn = 3) {
  
  ## Select replicate columns for that condition
  pattern <- paste0("^", gsub("\\.", "\\\\.", condition_prefix), "_[0-9]+$")
  cols <- grepl(pattern, colnames(expr_mat))
  sub <- expr_mat[, cols, drop = FALSE]
  
  cat("\n=== ", condition_prefix, " ===\n", sep="")
  cat("Replicates:", ncol(sub), "\n")
  
  if (ncol(sub) < min_reps_warn) {
    cat("NOTE: < ", min_reps_warn, " replicates. Correlations may be unreliable; report this.\n", sep="")
  }
  if (ncol(sub) < 2) stop("Not enough replicates to compute correlation.")
  
  ## Filter zero-variance genes
  gv <- apply(sub, 1, var)
  sub <- sub[gv > 0, , drop = FALSE]
  cat("Genes kept (var>0):", nrow(sub), "\n")
  if (nrow(sub) < 10) stop("Too few genes after filtering; check condition name or data.")
  
  ## Gene–gene correlation
  cor_gene <- cor(t(sub), method = cor_method, use = "pairwise.complete.obs")
  diag(cor_gene) <- 0
  
  ## Threshold by top_prop strongest |r|
  w <- abs(cor_gene[upper.tri(cor_gene)])
  thr <- as.numeric(quantile(w, probs = 1 - top_prop, na.rm = TRUE))
  
  idx <- which(abs(cor_gene) >= thr & upper.tri(cor_gene), arr.ind = TRUE)
  
  ## IMPORTANT: avoid naming column "weight" so igraph doesn't treat signed values as weights
  edges <- data.frame(
    source = rownames(cor_gene)[idx[,1]],
    target = colnames(cor_gene)[idx[,2]],
    corr   = cor_gene[idx],
    abs_corr = abs(cor_gene[idx])
  )
  
  ## Build graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  ## Connected components
  cc <- components(g)
  
  ## Centralities UNWEIGHTED (safe for signed correlations)
  deg <- degree(g)
  btw <- betweenness(g, normalized = TRUE)
  clo <- closeness(g, normalized = TRUE)
  clu <- transitivity(g, type = "local", isolates = "zero")
  
  ## Modules: use positive weights (abs correlation)
  E(g)$w_abs <- E(g)$abs_corr
  comm <- cluster_louvain(g, weights = E(g)$w_abs)
  
  node_table <- data.frame(
    gene = names(deg),
    degree = as.numeric(deg),
    betweenness = as.numeric(btw[names(deg)]),
    closeness = as.numeric(clo[names(deg)]),
    clustering_coeff = as.numeric(clu[names(deg)]),
    module = as.integer(membership(comm)[names(deg)]),
    component = as.integer(cc$membership[match(names(deg), V(g)$name)])
  )
  
  ## Export
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  edge_file <- file.path(out_dir, paste0("PA14_", condition_prefix, "_edges_top", top_prop*100, "pct.tsv"))
  node_file <- file.path(out_dir, paste0("PA14_", condition_prefix, "_nodes_top", top_prop*100, "pct.tsv"))
  
  write.table(edges[, c("source","target","corr")],
              edge_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(node_table,
              node_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## Print summary for your report
  cat("Implied |r| threshold:", round(thr, 4), "\n")
  cat("Nodes:", vcount(g), "Edges:", ecount(g), "\n")
  cat("Components:", cc$no, "Largest component:", max(cc$csize), "\n")
  cat("Saved:\n  ", edge_file, "\n  ", node_file, "\n", sep="")
  
  invisible(list(graph=g, edges=edges, nodes=node_table, thr=thr))
}

out_dir <- "cytoscape_condition_networks"
top_prop <- 0.05
cor_method <- "pearson"

networks <- lapply(conditions_4, function(cond) {
  build_condition_network(
    expr_mat = log_data_imp,
    condition_prefix = cond,
    top_prop = top_prop,
    cor_method = cor_method,
    out_dir = out_dir
  )
})
