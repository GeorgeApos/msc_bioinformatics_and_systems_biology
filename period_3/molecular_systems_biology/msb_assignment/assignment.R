############################################################
# Molecular Systems Biology Assignment – Dataset 19
# Author: George David Apostolidis
# Date: 2026-01-22
#
# Purpose:
# 1) Preliminary analysis / QC:
#    - dimensions, missingness, log-transform, normalization diagnostics
#    - sample clustering (dendrogram), correlation heatmaps
#    - PCA on imputed data + scree plot
# 2) Network inference:
#    - Build ONE global gene–gene co-expression "backbone" network using ALL samples
#    - Export global edge list + node topology for Cytoscape
# 3) Condition-specific subnetworks (in Cytoscape):
#    - Export per-condition node activity tables (mean expression, mean z-score, rank)
#    - In Cytoscape: filter active genes per condition and create subnetworks via selection
#
# Input:
# - CSV matrix: genes (rows) x samples (columns)
# - Sample names expected to be like: condition_1, condition_2, ...
#
# Output (for Cytoscape):
# - cytoscape_global_backbone/PA14_GLOBAL_edges.tsv
# - cytoscape_global_backbone/PA14_GLOBAL_nodes_topology.tsv
# - cytoscape_condition_activity_tables/PA14_<condition>_node_activity.tsv
# - cytoscape_condition_activity_tables/PA14_core_shared_genes_topN.tsv
############################################################

############################
# 0) Setup
############################
suppressPackageStartupMessages({
  library(igraph)
})

cat("\n============================================================\n")
cat("Molecular Systems Biology Assignment – Dataset 19\n")
cat("Author: George David Apostolidis\n")
cat("============================================================\n\n")

############################
# 1) Load data (genes x samples)
############################
infile <- "C:/Users/geoap/Documents/WUR/msc_bioinformatics_and_systems_biology/period_3/molecular_systems_biology/msb_assignment/datasets/dataset_19.csv"
cat("[1] Loading data from:\n", infile, "\n\n")

data <- as.matrix(read.csv(infile, row.names = 1, check.names = FALSE))

cat("Dimensions (genes x samples): ", nrow(data), " x ", ncol(data), "\n", sep = "")
cat("First 5 genes: ", paste(head(rownames(data), 5), collapse = ", "), "\n", sep = "")
cat("First 5 samples: ", paste(head(colnames(data), 5), collapse = ", "), "\n\n", sep = "")

############################
# 2) Basic counts + missingness
############################
cat("[2] Basic counts\n")
cat("Genes: ", nrow(data), "\n", sep = "")
cat("Samples: ", ncol(data), "\n\n", sep = "")

cat("[3] Checking missing/non-finite values (NA/NaN/Inf)\n")
n_na  <- sum(is.na(data))
n_nan <- sum(is.nan(data))
n_inf <- sum(is.infinite(data))
cat("NA:  ", n_na, "\n", sep = "")
cat("NaN: ", n_nan, "\n", sep = "")
cat("Inf: ", n_inf, "\n\n", sep = "")

genes_with_missing <- rownames(data)[rowSums(!is.finite(data)) > 0]
samples_with_missing <- colnames(data)[colSums(!is.finite(data)) > 0]
cat("Genes with any non-finite values: ", length(genes_with_missing), "\n", sep = "")
cat("Samples with any non-finite values: ", length(samples_with_missing), "\n\n", sep = "")

if (length(genes_with_missing) > 0 || length(samples_with_missing) > 0) {
  cat("Non-finite entries (row, col) preview:\n")
  idx <- which(!is.finite(data), arr.ind = TRUE)
  print(head(idx, 20))
  cat("\n")
}

############################
# 3) Log transform
############################
cat("[4] Log2 transform: log2(x + 1)\n\n")
log_data <- log2(data + 1)

############################
# 4) Impute missing values early (so density/PCA do not fail)
############################
cat("[5] Imputation (gene-wise mean) on log2 data\n")
log_data_imp <- log_data

for (i in 1:nrow(log_data_imp)) {
  missing <- !is.finite(log_data_imp[i, ])
  if (any(missing)) {
    log_data_imp[i, missing] <- mean(log_data_imp[i, !missing], na.rm = TRUE)
  }
}

missing_frac <- rowMeans(!is.finite(log_data))
log_data_imp <- log_data_imp[missing_frac <= 0.2, , drop = FALSE]

gene_var_imp <- apply(log_data_imp, 1, var)
log_data_imp <- log_data_imp[gene_var_imp > 0, , drop = FALSE]

cat("Genes retained after missingness filter (<=20%): ", nrow(log_data_imp), "\n", sep = "")
cat("Genes retained after removing zero-variance: ", nrow(log_data_imp), "\n\n", sep = "")

############################
# 5) Normalization sanity checks (visual)
############################
cat("[6] Plotting normalization diagnostics\n")
cat("    - Boxplots, density curves (NA-safe), column sums\n\n")

par(mar=c(8,4,2,1))
boxplot(log_data_imp, las=2, main="log2(expression+1) per sample (imputed)", ylab="log2(expr+1)")

plot(density(log_data_imp[,1], na.rm = TRUE),
     main="Density: log2(expression+1) across samples (imputed)",
     xlab="log2(expr+1)")
if (ncol(log_data_imp) > 1) {
  for (j in 2:ncol(log_data_imp)) lines(density(log_data_imp[,j], na.rm = TRUE))
}

totals <- colSums(data, na.rm = TRUE)
barplot(totals, las=2, main="Column sums (sanity check)", ylab="sum(expression)")

############################
# 6) Sample clustering + correlation heatmaps (preliminary analysis)
############################
cat("[7] Sample clustering (dendrogram) and correlation heatmaps\n\n")

dist_samples <- dist(t(log_data_imp), method = "euclidean")
hc_samples <- hclust(dist_samples, method = "complete")
plot(hc_samples, main="Hierarchical clustering of samples (imputed)", xlab="", sub="")

cor_samples <- cor(log_data_imp, method = "pearson", use = "pairwise.complete.obs")
cat("Top-left 5x5 of sample correlation matrix:\n")
print(round(cor_samples[1:5, 1:5], 2))
cat("\n")

heatmap(cor_samples,
        symm = TRUE,
        margins = c(6, 6),
        main = "Correlation between samples (Pearson)")

hc_samples2 <- hclust(dist(t(log_data_imp)))
heatmap(cor_samples,
        symm = TRUE,
        Colv = as.dendrogram(hc_samples2),
        Rowv = as.dendrogram(hc_samples2),
        margins = c(6,6),
        main = "Correlation heatmap with hierarchical clustering")

############################
# 7) Heatmap of top variable genes
############################
cat("[8] Heatmap of top variable genes\n\n")

gene_var <- apply(log_data_imp, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:50]

heatmap(log_data_imp[top_genes, ],
        scale = "row",
        margins = c(6, 6),
        main = "Top 50 most variable genes (row-scaled, imputed)")

############################
# 8) PCA on imputed data + scree plots
############################
cat("[9] PCA on imputed data + scree plots\n\n")

pca <- prcomp(t(log_data_imp), center = TRUE, scale. = FALSE)
var_explained <- summary(pca)$importance[2, ]

plot(pca$x[,1], pca$x[,2],
     pch = 19,
     xlab = paste0("PC1 (", round(100*var_explained[1],1), "%)"),
     ylab = paste0("PC2 (", round(100*var_explained[2],1), "%)"),
     main = "PCA of samples (imputed data)")
text(pca$x[,1], pca$x[,2], labels = colnames(log_data_imp), pos = 3, cex = 0.7)

plot(var_explained,
     type = "b",
     xlab = "Principal component",
     ylab = "Proportion of variance explained",
     main = "Scree plot (variance explained)")

cat("Cumulative variance explained (first 10 PCs):\n")
print(round(cumsum(var_explained)[1:min(10, length(var_explained))], 3))
cat("\n")

plot(pca$x[,1], pca$x[,3],
     xlab = paste0("PC1 (", round(100*var_explained[1],1), "%)"),
     ylab = paste0("PC3 (", round(100*var_explained[3],1), "%)"),
     pch = 19,
     main = "PCA: PC1 vs PC3")
text(pca$x[,1], pca$x[,3], labels = colnames(log_data_imp), pos = 3, cex = 0.7)

plot(pca$x[,2], pca$x[,3],
     xlab = paste0("PC2 (", round(100*var_explained[2],1), "%)"),
     ylab = paste0("PC3 (", round(100*var_explained[3],1), "%)"),
     pch = 19,
     main = "PCA: PC2 vs PC3")
text(pca$x[,2], pca$x[,3], labels = colnames(log_data_imp), pos = 3, cex = 0.7)

############################################################
# Part 3–4: Global backbone network + subnetwork extraction
############################################################

############################
# 9) GLOBAL BACKBONE NETWORK (all samples)
############################
cat("\n[10] Global backbone network inference (all samples)\n")
cat("     - Gene–gene Pearson correlation across all 47 samples\n")
cat("     - Threshold edges by top proportion of absolute correlations\n")
cat("     - Export global edge list + node topology for Cytoscape\n\n")

# Filter genes with zero variance across ALL samples to avoid undefined correlations
gene_var_all <- apply(log_data_imp, 1, var)
expr_all <- log_data_imp[gene_var_all > 0, , drop = FALSE]
cat("Genes kept for global network (var>0): ", nrow(expr_all), "\n", sep = "")
cat("Samples used: ", ncol(expr_all), "\n\n", sep = "")

# Correlation across all samples
cor_global <- cor(t(expr_all), method = "pearson", use = "pairwise.complete.obs")
diag(cor_global) <- 0

# Helper: keep top proportion of |correlation|
get_edges_by_top_prop <- function(cor_mat, prop){
  w <- abs(cor_mat[upper.tri(cor_mat)])
  thr <- as.numeric(quantile(w, probs = 1 - prop, na.rm = TRUE))
  idx <- which(abs(cor_mat) >= thr & upper.tri(cor_mat), arr.ind = TRUE)
  edges <- data.frame(
    source = rownames(cor_mat)[idx[,1]],
    target = colnames(cor_mat)[idx[,2]],
    corr   = cor_mat[idx],
    abs_corr = abs(cor_mat[idx])
  )
  list(edges = edges, thr = thr)
}

# Evaluate a few thresholds for justification
top_props <- c(0.02, 0.05, 0.10)

summ_global <- lapply(top_props, function(p){
  res <- get_edges_by_top_prop(cor_global, p)
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
summ_global_df <- do.call(rbind, summ_global)

cat("Global threshold evaluation table:\n")
print(summ_global_df)
cat("\n")

# Choose threshold (default: 0.05)
chosen_prop <- 0.05
res_global <- get_edges_by_top_prop(cor_global, chosen_prop)
edges_global <- res_global$edges
thr_global <- res_global$thr

cat("Chosen global top_prop: ", chosen_prop, "\n", sep = "")
cat("Implied global |r| threshold: ", round(thr_global, 6), "\n", sep = "")
cat("Global edges: ", nrow(edges_global), "\n", sep = "")
cat("Global nodes: ", length(unique(c(edges_global$source, edges_global$target))), "\n\n", sep = "")

# Build global graph + topology
g_global <- graph_from_data_frame(edges_global, directed = FALSE)

deg <- degree(g_global)
btw <- betweenness(g_global, normalized = TRUE)
clo <- closeness(g_global, normalized = TRUE)
clu <- transitivity(g_global, type = "local", isolates = "zero")

# Community detection using abs(correlation) as positive weights
E(g_global)$w_abs <- E(g_global)$abs_corr
comm <- cluster_louvain(g_global, weights = E(g_global)$w_abs)
cc <- components(g_global)

nodes_global <- data.frame(
  gene = V(g_global)$name,
  degree = as.numeric(deg[V(g_global)$name]),
  betweenness = as.numeric(btw[V(g_global)$name]),
  closeness = as.numeric(clo[V(g_global)$name]),
  clustering_coeff = as.numeric(clu[V(g_global)$name]),
  module = as.integer(membership(comm)[V(g_global)$name]),
  component = as.integer(cc$membership[match(V(g_global)$name, V(g_global)$name)])
)

cat("Global network topology summary:\n")
cat("Nodes: ", vcount(g_global), " | Edges: ", ecount(g_global), "\n", sep = "")
cat("Connected components: ", cc$no, " | Largest component size: ", max(cc$csize), "\n\n", sep = "")

# Export for Cytoscape
out_dir <- "cytoscape_global_backbone"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

global_edges_file <- file.path(out_dir, "PA14_GLOBAL_edges.tsv")
global_nodes_file <- file.path(out_dir, "PA14_GLOBAL_nodes_topology.tsv")

write.table(edges_global[, c("source","target","corr")],
            global_edges_file, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(nodes_global,
            global_nodes_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Saved Cytoscape global backbone files:\n")
cat("  ", global_edges_file, "\n", sep = "")
cat("  ", global_nodes_file, "\n\n", sep = "")

############################
# 10) CONDITION ACTIVITY TABLES (for Cytoscape subnetwork extraction)
############################
cat("[11] Exporting condition-specific node activity tables\n")

# Parse condition name (strip trailing _number)
cond_name <- sub("_[0-9]+$", "", colnames(log_data_imp))

# Select 4 conditions to compare (edit here)
conditions_4 <- c("iron.starvation", "stationary.phase", "biofilm.48h", "transitional.growth")
cat("Selected conditions for subnetwork extraction:\n")
print(conditions_4)
cat("\n")

# Gene-wise z-scores across all samples (relative activity)
z_all <- t(scale(t(log_data_imp)))
z_all[!is.finite(z_all)] <- 0

export_condition_activity <- function(cond, expr_mat, z_mat, cond_name_vec, out_dir){
  cols <- which(cond_name_vec == cond)
  if (length(cols) == 0) stop(paste("No samples found for condition:", cond))
  
  expr_sub <- expr_mat[, cols, drop = FALSE]
  z_sub <- z_mat[, cols, drop = FALSE]
  
  mean_expr <- rowMeans(expr_sub, na.rm = TRUE)
  mean_z <- rowMeans(z_sub, na.rm = TRUE)
  
  # Flags and rank to support multiple ways of selecting active genes
  active_ge_1 <- mean_z >= 1
  active_ge_0_5 <- mean_z >= 0.5
  rank_mean_z <- rank(-mean_z, ties.method = "min")  # 1 = highest mean_z
  
  tbl <- data.frame(
    gene = rownames(expr_mat),
    condition = cond,
    n_replicates = length(cols),
    mean_expr = as.numeric(mean_expr),
    mean_z = as.numeric(mean_z),
    rank_mean_z = as.integer(rank_mean_z),
    active_z_ge_1 = active_ge_1,
    active_z_ge_0_5 = active_ge_0_5
  )
  
  f <- file.path(out_dir, paste0("PA14_", cond, "_node_activity.tsv"))
  write.table(tbl, f, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Condition: ", cond,
      " | Replicates: ", length(cols),
      " | Saved: ", f, "\n", sep = "")
  invisible(tbl)
}

out_dir2 <- "cytoscape_condition_activity_tables"
dir.create(out_dir2, showWarnings = FALSE, recursive = TRUE)

activity_tables <- lapply(conditions_4, function(cond){
  export_condition_activity(cond, log_data_imp, z_all, cond_name, out_dir2)
})

# Core shared genes by "top N per condition" (more robust than strict z-thresholding)
topN <- 30
active_sets_topN <- lapply(activity_tables, function(tbl){
  tbl$gene[tbl$rank_mean_z <= topN]
})
core_shared_genes <- Reduce(intersect, active_sets_topN)

core_file <- file.path(out_dir2, "PA14_core_shared_genes_topN.tsv")
write.table(data.frame(gene = core_shared_genes),
            core_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nCore shared genes (intersection of top ", topN, " active genes per condition): ",
    length(core_shared_genes), "\n", sep = "")
cat("Saved core gene list: ", core_file, "\n\n", sep = "")

cat("============================================================\n")
cat("Script completed successfully.\n\n")
cat("============================================================\n")

