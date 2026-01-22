if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("remotes")
BiocManager::install("pachterlab/sleuth")

## ------------------------------------------------------------------------
library(DESeq2)
source("http://www.bioinformatics.nl/courses/RNAseq/DEseq2Exercise.R")

## ------------------------------------------------------------------------
expression_data=read.table(
  "http://www.bioinformatics.nl/courses/RNAseq/ST_vs_HT.csv", row.names=1, header=TRUE, sep =",", stringsAsFactors=FALSE) 

## ------------------------------------------------------------------------
dim(expression_data) 

## ------------------------------------------------------------------------
colnames(expression_data) 

## ------------------------------------------------------------------------
rownames(expression_data) 

## ------------------------------------------------------------------------
summary(expression_data)

## ------------------------------------------------------------------------
apply(expression_data, 2, sum)

## ------------------------------------------------------------------------
mx = apply( expression_data, 1, max )

## ------------------------------------------------------------------------
expression_data = expression_data[ mx > 10, ] 

## ------------------------------------------------------------------------
condition = factor(c("ST","ST","ST","HT","HT","HT"),c("ST","HT"))
col_data = data.frame(condition)
dds = DESeqDataSetFromMatrix(expression_data, col_data, ~condition)

## ------------------------------------------------------------------------
dds = estimateSizeFactors(dds) 

## ------------------------------------------------------------------------
sizeFactors(dds)

## ------------------------------------------------------------------------
norm_versus_non_norm( dds, 2, 4, left = 2, right = 8 )

## ------------------------------------------------------------------------
rld = rlog(dds)

## ------------------------------------------------------------------------
plot(density(assay(dds)[,1]), main="counts")
plot(density(assay(rld)[,1]), main="log counts")

## ------------------------------------------------------------------------
dists = dist(t(assay(rld)))

## ------------------------------------------------------------------------
plot(hclust(dists))  

## ------------------------------------------------------------------------
dds = estimateDispersions(dds)

## ------------------------------------------------------------------------
plotDispEsts(dds) 

## ------------------------------------------------------------------------
dds = nbinomWaldTest(dds)

## ------------------------------------------------------------------------
res = results(dds)

## ------------------------------------------------------------------------
head(res)

## ------------------------------------------------------------------------
res$padj = ifelse(is.na(res$padj), 1, res$padj)

write.table(res, col.names=NA, row.names=T, file ="expressions.tsv", sep ="\t")

## ------------------------------------------------------------------------
plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)

## ------------------------------------------------------------------------
library(sleuth)

## ------------------------------------------------------------------------
setwd("~/kallisto") # or where you put the downloaded kallisto folder

## ------------------------------------------------------------------------
base_dir = "."
sample_ids = dir(base_dir,"_")
kal_dirs = sapply(sample_ids, function(id) file.path(base_dir, id))

## ------------------------------------------------------------------------
s2c = read.table(file.path(base_dir,"design.txt"), header=TRUE, stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
s2c$path = kal_dirs

## ------------------------------------------------------------------------
print(s2c)

## ------------------------------------------------------------------------
so = sleuth_prep(s2c, ~condition,extra_bootstrap_summary=TRUE)

## ------------------------------------------------------------------------
so = sleuth_fit(so, ~condition, 'full')

## ------------------------------------------------------------------------
so = sleuth_fit(so, ~1, 'reduced')

## ------------------------------------------------------------------------
so = sleuth_lrt(so, 'reduced', 'full')

## ------------------------------------------------------------------------
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE) 

## ------------------------------------------------------------------------
sleuth_significant <- sleuth_table[sleuth_table$qval<=0.05,]
head(sleuth_significant)

## ------------------------------------------------------------------------
plot_bootstrap(so, "AT1G56600.1")

## ------------------------------------------------------------------------
plot_bootstrap(so,"AT1G01810.1")

## ------------------------------------------------------------------------
so = sleuth_wt(so,which_beta="conditionST", which_model="full")

## ------------------------------------------------------------------------
sleuth_table <- sleuth_results(so, 'conditionST', 'wt', show_all=FALSE)

## ------------------------------------------------------------------------
sleuth_live(so)

## ------------------------------------------------------------------------
sleuth_matrix = sleuth_to_matrix(so, 'obs_norm', 'tpm')
transcripts = rownames(sleuth_matrix)
t2g = data.frame(target_id=transcripts, gene_id=substring(transcripts,1,9),stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
so = sleuth_prep(s2c, ~ condition,target_mapping = t2g, aggregation_column = 'gene_id')
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_wt(so,which_beta="conditionST", which_model="full")
sleuth_table <- sleuth_results(so, 'conditionST', 'wt', show_all=FALSE)
sleuth_significant <- sleuth_table[sleuth_table$qval<=0.05,]

## ------------------------------------------------------------------------
deseq2_significant = data.frame(res[res$padj<=0.05,])

## ------------------------------------------------------------------------
both_significant = merge(deseq2_significant,sleuth_significant,by.x=0,by.y="target_id",all=FALSE)

## ------------------------------------------------------------------------
nrow(deseq2_significant)
nrow(sleuth_significant)
nrow(both_significant)
